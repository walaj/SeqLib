#include "SeqLib/BWAAligner.h"
#include <htslib/sam.h>

extern "C" {
  // signature must match exactly how theyre defined in bwamem_extra.c / bwamem.c
  mem_alnreg_v mem_align1_core(const mem_opt_t *opt,
                                const bwt_t *bwt,
                                const bntseq_t *bns,
                                const uint8_t *pac,
                                int l_seq,
                                char *seq,
                                void *buf);

  void mem_mark_primary_se(const mem_opt_t *opt,
                           int n,
                           mem_alnreg_t *a,
                           int64_t id);
}

namespace SeqLib {

  namespace {
    // sort primary by descending MAPQ, then by (rid,pos)
    bool aln_sort(const mem_aln_t& a, const mem_aln_t& b) {
      if (a.mapq != b.mapq)    return a.mapq > b.mapq;
      if (a.rid  != b.rid)     return a.rid  < b.rid;
      return a.pos    < b.pos;
    }
  }
    
  void BWAAligner::SetGapOpen(int gap_open) {
    if (gap_open < 0) {
      throw std::invalid_argument{"SetGapOpen: gap_open must be >= 0"};
    }
    memopt_->o_ins = memopt_->o_del = gap_open;
  }
  
  void BWAAligner::SetGapExtension(int gap_ext) {
    if (gap_ext < 0) {
      throw std::invalid_argument{"SetGapExtension: gap_ext must be >= 0"};
    }
    memopt_->e_ins = memopt_->e_del = gap_ext;
  }
  
  void BWAAligner::SetMismatchPenalty(int mismatch) {
    if (mismatch < 0) {
      throw std::invalid_argument{"SetMismatchPenalty: mismatch must be >= 0"};
    }
    memopt_->b = mismatch;
    bwa_fill_scmat(memopt_->a, memopt_->b, memopt_->mat);
  }
  
  void BWAAligner::SetZDropoff(int zdrop) {
    if (zdrop < 0) {
      throw std::invalid_argument{"SetZDropoff: zdrop must be >= 0"};
    }
    memopt_->zdrop = zdrop;
  }
  
  void BWAAligner::SetAScore(int a) {
    if (a < 0) {
      throw std::invalid_argument{"SetAScore: a must be >= 0"};
    }
    // scale all related penalties by a
    memopt_->a            = a;
    memopt_->b           *= a;
    memopt_->T           *= a;
    memopt_->o_ins       *= a;
    memopt_->o_del       *= a;
    memopt_->e_ins       *= a;
    memopt_->e_del       *= a;
    memopt_->zdrop       *= a;
    memopt_->pen_clip5   *= a;
    memopt_->pen_clip3   *= a;
    memopt_->pen_unpaired*= a;
  }
  
  void BWAAligner::Set3primeClippingPenalty(int penalty) {
    if (penalty < 0) {
      throw std::invalid_argument{"Set3primeClippingPenalty: penalty must be >= 0"};
    }
    memopt_->pen_clip3 = penalty;
  }
  
  void BWAAligner::Set5primeClippingPenalty(int penalty) {
    if (penalty < 0) {
      throw std::invalid_argument{"Set5primeClippingPenalty: penalty must be >= 0"};
    }
    memopt_->pen_clip5 = penalty;
  }
  
  void BWAAligner::SetBandwidth(int bw) {
    if (bw < 0) {
      throw std::invalid_argument{"SetBandwidth: bandwidth must be >= 0"};
    }
    memopt_->w = bw;
  }
  
  void BWAAligner::SetReseedTrigger(float trigger) {
    if (trigger < 0.0f) {
      throw std::invalid_argument{"SetReseedTrigger: trigger must be >= 0"};
    }
    memopt_->split_factor = trigger;
  }

  void BWAAligner::allocBuffer(size_t read_len) {
    read_len_ = read_len;
    seq_buf = (char*) malloc(read_len_);
    core_buf = calloc(1, 1 << 16);  // 64KB internal working space (can tune)
  }

  mem_alnreg_v BWAAligner::mem_align1_reuse(const mem_opt_t* opt,
				const bwt_t* bwt,
				const bntseq_t* bns,
				const uint8_t* pac,
				int l_seq,
				const char* seq_in,
				char* seq_buf,
				void* core_buf) {
    memcpy(seq_buf, seq_in, l_seq);  // cheaper since you reuse the same buffer
    mem_alnreg_v ar = mem_align1_core(opt, bwt, bns, pac, l_seq, seq_buf, core_buf);
    mem_mark_primary_se(opt, ar.n, ar.a, 0); //lrand48());
    return ar;
  }

  bseq1_t BWAAligner::convertToBseq(const UnalignedSequence& us) {
    bseq1_t bs;
    // 1) qname
    bs.name    = strdup(us.Name.c_str());
    bs.comment = nullptr;
    
    // 2) seq: 4-bit encoded using seq_nt16_table[]
    bs.l_seq = us.Seq.size();
    bs.seq   = (char*)malloc(bs.l_seq);
    memcpy(bs.seq, us.Seq.data(), bs.l_seq);    
    
    // 3) qual: Phred-scaled, robust to empty
    bs.qual   = nullptr;
    
    bs.sam = nullptr;
      
    return bs;
  }

std::vector<Alignment> BWAAligner::alignBatch(
			 const UnalignedSequenceVector& inputs) {
  
  int n = inputs.size();
  // 1) Build array of bseq1_t
  std::vector<bseq1_t> bseqs(n);
  for(int i = 0; i < n; ++i)
    bseqs[i] = convertToBseq(inputs[i]);

  // 2) Call the batch aligner
  mem_process_seqs(
    memopt_,
    index_->idx_->bwt,
    index_->idx_->bns,
    index_->idx_->pac,
    processedCount_,  // number of reads already done
    n,
    bseqs.data(),
    /* pes0 = */ nullptr   // let BWA infer fragment statistics
		   );
  
    // 3) Parse the SAM text in bseqs[i].sam into Alignment structs
  std::vector<Alignment> out;
  out.reserve(n * 2);

      for (int i = 0; i < n; ++i) {
      char* ptr = bseqs[i].sam;
      if (!ptr) continue;

      while (*ptr) {
        // grab one line
        char* nl = strchr(ptr, '\n');
        size_t len = nl ? (nl - ptr) : strlen(ptr);
        std::string line(ptr, len);
        ptr = nl ? nl + 1 : ptr + len;
        if (line.empty()) continue;

        // split on tabs
        std::istringstream iss(line);
        std::vector<std::string> f;
        std::string fld;
        while (std::getline(iss, fld, '\t')) f.push_back(fld);

        // 1) flag is always numeric
        int flag = std::stoi(f[1]);
        // 2) skip unmapped reads
        if (flag & 0x4) continue;

        // 3) now its safe to parse POS and MAPQ
        int64_t pos = std::stoll(f[3]) - 1;   // SAM is 1-based
        int     mapq = std::stoi(f[4]);
        std::string cigar = f[5];

        // 4) map contig name (f[2]) to an integer rid if you really need one
        //    or store it as string if thats sufficient.
        int rid = -1;
        std::string rname = f[2];
        for (int j = 0; j < index_->idx_->bns->n_seqs; ++j) {
          if (rname == index_->idx_->bns->anns[j].name) { rid = j; break; }
        }

        // 5) strand
        bool is_rev = (flag & 0x10) != 0;

        // 6) push your Alignment
        Alignment A;
        A.qname  = f[0];
        A.rid    = rid;      // or store rname if you change your struct
        A.pos    = pos;
        A.mapq   = mapq;
        A.cigar  = std::move(cigar);
        A.is_rev = is_rev;
        out.push_back(std::move(A));
      }
    
    // 4) clean up this bseq1_t
    free(bseqs[i].name);
    free(bseqs[i].seq);
    if (bseqs[i].qual) free(bseqs[i].qual);
    free(bseqs[i].sam);
  }
  
  // 5) advance counter
  processedCount_ += n;
  return out;
}
  
  // void BWAAligner::alignBatch(const UnalignedSequenceVector& inputs,
  // 			      BamRecordPtrVector& out,
  // 			      bool hardclip, double keepSecFrac, int maxSecondary) const
  // {
  //   int batchSize = inputs.size();
  //   std::vector<bam1_t*> bamBatch(batchSize);
  //   // 1) Convert your UnalignedSequence bam1_t* for each read
  //   for(int i = 0; i < batchSize; ++i) {
  //     bamBatch[i] = convertToBam1(inputs[i]); 
  //   }
    
  //   // 2) Call the batched aligner
  //   int64_t id = 0; // ok this is a seed for breaking "ties"
  //   // for primary vesrus secondary. I don't care that much about
  //   // these reads to will set as a fixed number (zero) rather than
  //   // allowing for randomness
    
  //   mem_alnreg_v* regs = mem_process_seqs(
  //       memopt_,
  //       index_->idx_->bwt,
  //       index_->idx_->bns,
  //       index_->idx_->pac,
  //       processedCount_,
  //       batchSize,
  //       bamBatch.data(),
  //       id
  //   );

  //   // 3) For each read, there may be multiple regions (regs[i].n)
  //   for(size_t i = 0; i < batchSize; ++i) {
  //       auto& vr = regs[i];

  //       // The original unmapped record (we'll clone it per region)
  //       bam1_t* original = bamBatch[i];

  //       for(int r = 0; r < vr.n; ++r) {
  //           // 4) Clone the template record
  //           bam1_t* aligned = bam_dup1(original);

  //           // 5) Populate it for *only that one* region
  //           //    - n = 1, a = &vr.a[r]
  //           mem_reg2aln(
  //               memopt_,
  //               1,
  //               &vr.a[r],
  //               aligned,
  //               hardclip ? 1 : 0,
  //               keepSecFrac
  //           );
  //           // note: mem_reg2aln does not write SA tags; it's one region only

  //           // 6) Wrap in your BamRecordPtr and store
  //           out.emplace_back(std::make_shared<SeqLib::BamRecord>(aligned));

  //           // free the bam1_t inside the shared_ptr
  //           // (SeqLib::BamRecords destructor should do bam_destroy1 for you)
  //       }

  //       // clean up this reads regions
  //       free(vr.a);
  //       // and the original unmapped record
  //       bam_destroy1(original);
  //   }

  //   free(regs);
  //   processedCount_ += batchSize;
  // }
  
  
  void BWAAligner::alignSequence(const std::string& seq,
				 const std::string& name,
				 BamRecordPtrVector& out,
				 bool hardclip,
				 double keepSecFrac,
				 int maxSecondary)
  {

    //assert(out.empty());
    //out.clear();
    
    // nothing to do if no index
    if (index_->IsEmpty()) return;
    
    // run BWA-MEM core
    // auto t0 = std::chrono::high_resolution_clock::now();
    
    // mem_alnreg_v regs = mem_align1_reuse(memopt_,
    // 					 index_->idx_->bwt,
    // 					 index_->idx_->bns,
    // 					 index_->idx_->pac,
    // 					 seq.size(),
    // 					 seq.data(),
    // 					 seq_buf,
    // 					 core_buf);

  // mem_alnreg_t:
  //   int64_t rb, re; // [rb,re): reference sequence in the alignment
  //   int qb, qe;     // [qb,qe): query sequence in the alignment
  //   int rid;        // reference seq ID
  //   int score;      // best local SW score
  //   int truesc;     // actual score corresponding to the aligned region; possibly smaller than $score
  //   int sub;        // 2nd best SW score
  //   int alt_sc;
  //   int csub;       // SW score of a tandem hit
  //   int sub_n;      // approximate number of suboptimal hits
  //   int w;          // actual band width used in extension
  //   int seedcov;    // length of regions coverged by seeds
  //   int secondary;  // index of the parent hit shadowing the current hit; <0 if primary
  //   int secondary_all;
  //   int seedlen0;   // length of the starting seed
  //   int n_comp:30, is_alt:2; // number of sub-alignments chained together
  //   float frac_rep;
  //   uint64_t hash;
    mem_alnreg_v regs = mem_align1(memopt_,
				   index_->idx_->bwt,
				   index_->idx_->bns,
				   index_->idx_->pac,
				   seq.size(),
				   seq.data());
    // auto t1 = std::chrono::high_resolution_clock::now();
    // std::cout << "Part A: "
    //           << std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count()
    //           << " ns\n";
    
    double primaryScore = 0;
    int secondaryCount = 0;

//     typedef struct { // This struct is only used for the convenience of API.
// 	int64_t pos;     // forward strand 5'-end mapping position
// 	int rid;         // reference sequence index in bntseq_t; <0 for unmapped
// 	int flag;        // extra flag
// 	uint32_t is_rev:1, is_alt:1, mapq:8, NM:22; // is_rev: whether on the reverse strand; mapq: mapping quality; NM: edit distance
// 	int n_cigar;     // number of CIGAR operations
// 	uint32_t *cigar; // CIGAR in the BAM encoding: opLen<<4|op; op to integer mapping: MIDSH=>01234
// 	char *XA;        // alternative mappings

// 	int score, sub, alt_sc;
// } mem_aln_t;
    std::vector<mem_aln_t> hits;
    hits.reserve(regs.n);

    for (int i = 0; i < regs.n; ++i) {
      auto& r = regs.a[i];

      // don't overdo seconarry
      if (r.secondary >= 0) { // <0 is primay
	++secondaryCount;
	if (secondaryCount > maxSecondary)
	  continue;
      }
	
      hits.push_back(mem_reg2aln(memopt_,
				 index_->idx_->bns,
				 index_->idx_->pac,
				 seq.size(),
				 seq.data(),
				 &r));
    }

    std::free(regs.a);
    std::sort(hits.begin(), hits.end(), aln_sort);

    // emit a BamRecord for each hit
      
    for (size_t i = 0; i < hits.size(); ++i) {
      auto& h = hits[i];
      
      // set the score for the primary
      bool isSecondary   = (h.flag & BAM_FSECONDARY);      
      if (!isSecondary)
	primaryScore = h.score;
      
      // don't add if the alignment score a secondary is too low
      if (isSecondary && (h.score < primaryScore * keepSecFrac)) {
	std::free(h.cigar);
	std::free(h.XA);
	continue;
      }
      
      auto b = std::make_shared<BamRecord>();
      
      b->b->core.tid     = h.rid;
      b->b->core.pos     = h.pos;
      b->b->core.qual    = h.mapq;
      b->b->core.flag    = h.flag;
      b->b->core.n_cigar = h.n_cigar;
      b->b->core.mtid    = -1;
      b->b->core.mpos    = -1;
      b->b->core.isize   = 0;
      if (h.is_rev)
	b->b->core.flag |= BAM_FREVERSE;

      // optionally hard clip out
      std::string clipped = seq;
      if (hardclip) {
	size_t tstart = 0, clen = 0;
	for (int c = 0; c < h.n_cigar; ++c) {
	  auto op = bam_cigar_op(h.cigar[c]);
	  if (c == 0 && op == BAM_CREF_SKIP)
	    tstart = bam_cigar_oplen(h.cigar[c]);
	  else if (bam_cigar_type(op) & 1)
	    clen += bam_cigar_oplen(h.cigar[c]);
	}
	assert(clen && tstart + clen <= seq.size());
	clipped = seq.substr(tstart, clen);
      }
      
      b->b->core.l_qname = name.size() + 1;
      b->b->core.l_qseq  = clipped.size();
      b->b->l_data       = b->b->core.l_qname +
	(h.n_cigar << 2) +
	((b->b->core.l_qseq + 1) >> 1) +
	b->b->core.l_qseq;
      
      b->b->data = static_cast<uint8_t*>(std::malloc(b->b->l_data));
      if (!b->b->data)
	throw std::bad_alloc();
      
      std::memcpy(b->b->data, name.c_str(), name.size() + 1); //qname
      std::memcpy(b->b->data + b->b->core.l_qname, h.cigar, h.n_cigar << 2); //cigar
      
      int newOp = hardclip ? BAM_CHARD_CLIP : BAM_CSOFT_CLIP;
      uint8_t* raw = reinterpret_cast<uint8_t*>(bam_get_cigar(b->b));
      for (size_t k = 0; k < b->b->core.n_cigar; ++k) {
	uint32_t cigarOp;
	std::memcpy(&cigarOp, raw + k * sizeof(uint32_t), sizeof(cigarOp));
	if ((cigarOp & BAM_CIGAR_MASK) == BAM_CREF_SKIP) {
	  cigarOp = (cigarOp & ~BAM_CIGAR_MASK) | static_cast<uint32_t>(newOp);
	  std::memcpy(raw + k * sizeof(uint32_t), &cigarOp, sizeof(cigarOp));
	}
      }
      
      std::free(h.cigar);
      
      auto* seqbuf = b->b->data + b->b->core.l_qname + (b->b->core.n_cigar << 2);
      int sl = clipped.size();
      if (h.is_rev) {
	int j = 0;
	for (int p = sl - 1; p >= 0; --p, ++j) {
	  uint8_t v = 15; // N
	  switch (clipped[p]) {
	  case 'A': v = 8; break; // T
	  case 'C': v = 4; break; // G
	  case 'G': v = 2; break; // C
	  case 'T': v = 1; break; // A
	  case 'N': v = 15; break;
	  }
	  seqbuf[j >> 1] &= ~(0xF << ((~j & 1) << 2));
	  seqbuf[j >> 1] |= v << ((~j & 1) << 2);
	}
      } else {
	for (int p = 0; p < sl; ++p) {
	  uint8_t v = 15; // N
	  switch (clipped[p]) {
	  case 'A': v = 1; break;
	  case 'C': v = 2; break;
	  case 'G': v = 4; break;
	  case 'T': v = 8; break;
	  case 'N': v = 15; break;
	  }
	  seqbuf[p >> 1] &= ~(0xF << ((~p & 1) << 2));
	  seqbuf[p >> 1] |= v << ((~p & 1) << 2);
	}
      }      
      
      auto* q = bam_get_qual(b->b);
      if (h.is_rev) {
	for (int j = 0, p = sl - 1; p >= 0; --p, ++j) {
	  q[j] = 0x23; // or your real qualities; 0xFF is special ("missing")
	}
      } else {
	for (int p = 0; p < sl; ++p) {
	  q[p] = 0x23; // placeholder; set real per-base Q if you have it
	}
      }
      
      //q[0] = 0xff; //old
      
      // add score tags
      b->AddIntTag("NM", h.NM);
      if (h.XA)
	b->AddZTag("XA", std::string(h.XA));
      b->AddIntTag("AS", h.score);
      b->AddIntTag("XS", h.sub);
      
      std::free(h.XA);
      
      out.push_back(b);
    }
    // auto t2 = std::chrono::high_resolution_clock::now();
    // std::cout << "Part B: "
    //           << std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count()
    //           << " ns\n";
    
    
  }
  
  void BWAAligner::alignSequence(const UnalignedSequence& us,
				 BamRecordPtrVector&           out,
				 bool                       hardclip,
				 double                     keepSecFrac,
				 int                        maxSecondary)
  {
    // delegate then optionally append BC tag
    alignSequence(us.Seq, us.Name, out, hardclip, keepSecFrac, maxSecondary);
    if (!copyComment_) return;
    for (auto& rec : out)
      rec->AddZTag("BC", us.Com);
  }
  
}
