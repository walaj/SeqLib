#include "SeqLib/BWAAligner.h"

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

void BWAAligner::alignSequence(const std::string& seq,
                               const std::string& name,
                               BamRecordVector& out,
                               bool hardclip,
                               double keepSecFrac,
                               int maxSecondary) const
{
  // nothing to do if no index
  if (index_->IsEmpty()) return;

  // run BWA-MEM core
  mem_alnreg_v regs = mem_align1(memopt_, index_->idx_->bwt, index_->idx_->bns,
				 index_->idx_->pac,
                                 seq.size(), seq.data());

  double primaryScore = 0;
  int   secondaryCount = 0;

  // filter & convert to mem_aln_t
  std::vector<mem_aln_t> hits;
  hits.reserve(regs.n);
  for (int i = 0; i < regs.n; ++i) {
    auto& r = regs.a[i];
    
    // skip secondaries if frac param invalid
    if ((r.secondary) &&
        (keepSecFrac < 0.0 || keepSecFrac > 1.0))
      continue;
    hits.push_back(mem_reg2aln(memopt_,
                               index_->idx_->bns,
                               index_->idx_->pac,
                               seq.size(),
                               seq.data(),
                               &r));
  }
  std::free(regs.a);

  // sort
  std::sort(hits.begin(), hits.end(), aln_sort);

  // emit BamRecord for each hit
  for (size_t i = 0; i < hits.size(); ++i) {
    auto& h = hits[i];
    bool isSec   = (h.flag & BAM_FSECONDARY);
    bool tooLow  = isSec && (primaryScore * keepSecFrac > h.score);
    bool tooMany = isSec && (int(i) > maxSecondary);
    if (tooLow || tooMany) {
      std::free(h.cigar);
      continue;
    }
    if (!isSec)
      primaryScore = h.score;

    BamRecord b;
    b.init();
    b.b->core.tid       = h.rid;
    b.b->core.pos       = h.pos;
    b.b->core.qual      = h.mapq;
    b.b->core.flag      = h.flag;
    b.b->core.n_cigar   = h.n_cigar;
    b.b->core.mtid      = -1;
    b.b->core.mpos      = -1;
    b.b->core.isize     = 0;
    if (h.is_rev) b.b->core.flag |= BAM_FREVERSE;

    // optionally hardclip out
    std::string clipped = seq;
    if (hardclip) {
      size_t tstart = 0, clen = 0;
      for (int c = 0; c < h.n_cigar; ++c) {
        auto op = bam_cigar_op(h.cigar[c]);
        if (c == 0 && op == BAM_CREF_SKIP)
          tstart = bam_cigar_oplen(h.cigar[c]);
        else if (bam_cigar_type(op)&1)
          clen += bam_cigar_oplen(h.cigar[c]);
      }
      assert(clen && tstart+clen <= seq.size());
      clipped = seq.substr(tstart, clen);
    }

    // allocate & fill
    b.b->core.l_qname = name.size()+1;
    b.b->core.l_qseq  = clipped.size();
    b.b->l_data       = b.b->core.l_qname
                      + (h.n_cigar<<2)
                      + ((b.b->core.l_qseq+1)>>1)
                      + b.b->core.l_qseq;
    b.b->data  = static_cast<uint8_t*>(std::malloc(b.b->l_data));
    if (!b.b->data) {
      throw std::bad_alloc();
    }
    //b.b->data = uint8_t(malloc(b.b->l_data));

    // qname
    std::memcpy(b.b->data, name.c_str(), name.size()+1);
    // CIGAR
    std::memcpy(b.b->data + b.b->core.l_qname,
                h.cigar, h.n_cigar<<2);
    // convert N to S/H
    int newOp = hardclip ? BAM_CHARD_CLIP : BAM_CSOFT_CLIP;
    auto* cig = bam_get_cigar(b.b);
    for (int k = 0; k < b.b->core.n_cigar; ++k) {
      if ((cig[k]&BAM_CIGAR_MASK)==BAM_CREF_SKIP) {
        cig[k] &= ~BAM_CIGAR_MASK;
        cig[k] |= newOp;
      }
    }
    std::free(h.cigar);

    // sequence
    auto* seqbuf = b.b->data
                 + b.b->core.l_qname
                 + (b.b->core.n_cigar<<2);
    int sl = clipped.size();
    if (h.is_rev) {
      int j=0;
      for (int p=sl-1; p>=0; --p, ++j) {
        uint8_t v=15;
        switch (clipped[p]) {
          case 'A': v=8; break;
          case 'C': v=2; break;
          case 'G': v=4; break;
          case 'T': v=1; break;
        }
        seqbuf[j>>1] &= ~(0xF<<((~j&1)<<2));
        seqbuf[j>>1] |= v<<((~j&1)<<2);
      }
    } else {
      for (int p=0; p<sl; ++p) {
        uint8_t v=15;
        switch (clipped[p]) {
          case 'A': v=1; break;
          case 'C': v=2; break;
          case 'G': v=4; break;
          case 'T': v=8; break;
        }
        seqbuf[p>>1] &= ~(0xF<<((~p&1)<<2));
        seqbuf[p>>1] |= v<<((~p&1)<<2);
      }
    }

    // qual = NULL
    auto* q = bam_get_qual(b.b);
    q[0] = 0xff;

    // tags
    b.AddIntTag("NA", regs.n);
    b.AddIntTag("NM", h.NM);
    if (h.XA) b.AddZTag("XA", std::string(h.XA));
    b.AddIntTag("AS", h.score);
    if (b.SecondaryFlag()) ++secondaryCount;

    out.push_back(std::move(b));
  }

  // add secondary counts
  for (auto& rec : out)
    rec.AddIntTag("SQ", secondaryCount);
}

void BWAAligner::alignSequence(const UnalignedSequence& us,
                               BamRecordVector&           out,
                               bool                       hardclip,
                               double                     keepSecFrac,
                               int                        maxSecondary) const
{
  // delegate then optionally append BC tag
  alignSequence(us.Seq, us.Name, out, hardclip, keepSecFrac, maxSecondary);
  if (!copyComment_) return;
  for (auto& rec : out)
    rec.AddZTag("BC", us.Com);
}
  
}
