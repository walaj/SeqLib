#include "SeqLib/BWAIndex.h"

#include <memory>
#include <vector>
#include <numeric>
#include <cstdlib>     // for std::calloc, std::free
#include <stdexcept>
#include <filesystem>

#include "SeqLib/UnalignedSequence.h"

#define _set_pac(pac, l, c) ((pac)[(l)>>2] |= (c)<<((~(l)&3)<<1))
#define _get_pac(pac, l) ((pac)[(l)>>2]>>((~(l)&3)<<1)&3)

namespace SeqLib {
  
  BWAIndex::~BWAIndex() {
    if (idx_) {
      bwa_idx_destroy(idx_);
      idx_ = nullptr;
    }
  }

  bool BWAIndex::IsEmpty() const noexcept {
    return idx_ == nullptr;
  }  
  
  void BWAIndex::LoadIndex(const std::string& fastaPath) {
    auto newIdx = bwa_idx_load(fastaPath.c_str(), BWA_IDX_ALL);
    if (!newIdx) throw std::runtime_error("Failed to load BWA index");
    if (idx_) bwa_idx_destroy(idx_);
    idx_ = newIdx;
  }
  
  BamHeader BWAIndex::HeaderFromIndex() const 
  {
    std::string my_hdr = printSamHeader();
    
    BamHeader hdr(my_hdr);
    
    return hdr;
  }
  
  int BWAIndex::NumSequences() const {
    
    if (!idx_)
      return 0;
    
    return idx_->bns->n_seqs;
    
  }
  
  std::string BWAIndex::ChrIDToName(int id) const {
    
    if (!idx_)
      throw std::runtime_error("Index has not be loaded / constructed");
    if (id < 0 || id >= idx_->bns->n_seqs) 
      throw std::out_of_range("BWAIndex::ChrIDToName - id out of bounds of refs in index for id of " + tostring(id) + " on IDX of size " + tostring(idx_->bns->n_seqs));
    
    return std::string(idx_->bns->anns[id].name);
  }
  
  std::string BWAIndex::printSamHeader() const {
    if (!idx_ || !idx_->bns)
      return "";
    
    std::ostringstream out;
    const bntseq_t* bns = idx_->bns;
    
    // one @SQ line per contig
    for (int i = 0; i < bns->n_seqs; ++i) {
      out << "@SQ\tSN:" << bns->anns[i].name
	  << "\tLN:" << bns->anns[i].len
	  << "\n";
    }
    
    return out.str();
  }
  
  /// Build an in-memory BWA index from the given reference sequences.
  /// Throws std::invalid_argument if any sequence name or seq is empty,
  /// std::bad_alloc on allocation failure.
  void BWAIndex::ConstructIndex(const UnalignedSequenceVector& refs) {
    if (refs.empty()) return;
    
    // integrity check
    for (auto const& r : refs) {
      if (r.Name.empty() || r.Seq.empty())
        throw std::invalid_argument(
				    "BWAIndex::Construct each reference must have non-empty Name and Seq"
				    );
    }
    
    // clear old index if present
    if (idx_) {
      bwa_idx_destroy(idx_);
      idx_ = nullptr;
    }
    
    // allocate the wrapper struct
    idx_ = static_cast<bwaidx_t*>(
				  std::calloc(1, sizeof(bwaidx_t))
				  );
    if (!idx_) throw std::bad_alloc();
    
    // build the forward only PAC
    auto* fwd_pac = seqlib_make_pac(refs, /*for_only=*/true);
    if (!fwd_pac) {
      throw std::runtime_error("BWAIndex::Construct failed to build forward PAC");
    }
    
    // build the full PAC (forward+reverse) for BWT construction
    auto* pac = seqlib_make_pac(refs, /*for_only=*/false);
    if (!pac) {
      std::free(fwd_pac);
      throw std::runtime_error("BWAIndex::Construct failed to build full PAC");
    }
    
    // total reference length
    const size_t total_len = std::accumulate(
					     refs.begin(), refs.end(),
					     size_t{0},
					     [](size_t sum, auto const& r){ return sum + r.Seq.size(); }
					     );
    
    // build the BWT
    bwt_t* bwt = seqlib_bwt_pac2bwt(pac, total_len * 2);
    if (!bwt) {
      std::free(pac);
      std::free(fwd_pac);
      throw std::runtime_error("BWAIndex::Construct BWT construction failed");
    }
    bwt_bwtupdate_core(bwt);
    std::free(pac);
    
    // finalize suffix array and count table
    bwt_cal_sa(bwt, /*sa_interval=*/32);
    bwt_gen_cnt_table(bwt);
    
    // build the name/length annotations (bns)
    auto* bns = static_cast<bntseq_t*>(
				       std::calloc(1, sizeof(bntseq_t))
				       );
    if (!bns) {
      // clean up everything we allocated
      free(bwt);
      std::free(fwd_pac);
      throw std::bad_alloc();
    }
    bns->l_pac    = total_len;
    bns->n_seqs   = refs.size();
    bns->seed     = 11;   // fixed seed
    bns->n_holes  = 0;
    bns->anns     = static_cast<bntann1_t*>(
					    std::calloc(refs.size(), sizeof(bntann1_t))
					    );
    if (!bns->anns) {
      std::free(bns);
      std::free(fwd_pac);
      throw std::bad_alloc();
    }
    
    // fill the annotations
    size_t offset = 0;
    for (size_t i = 0; i < refs.size(); ++i) {
      seqlib_add_to_anns(
			 refs[i].Name,
			 refs[i].Seq,
			 &bns->anns[i],
			 offset
			 );
      offset += refs[i].Seq.size();
    }
    bns->ambs = nullptr;  // no 'holes' array at this point
    
    // hand off to idx
    idx_->bwt = bwt;
    idx_->bns = bns;
    idx_->pac = fwd_pac;
  }
  
  // modified from bwa (heng li)
  uint8_t* BWAIndex::seqlib_add1(const kseq_t *seq, bntseq_t *bns, uint8_t *pac, int64_t *m_pac, int *m_seqs, int *m_holes, bntamb1_t **q) const
  {
    bntann1_t *p;
    int lasts;
    if (bns->n_seqs == *m_seqs) {
      *m_seqs <<= 1;
      bns->anns = (bntann1_t*)realloc(bns->anns, *m_seqs * sizeof(bntann1_t));
    }
    p = bns->anns + bns->n_seqs;
    p->name = strdup((char*)seq->name.s);
    p->anno = seq->comment.l > 0? strdup((char*)seq->comment.s) : strdup("(null)");
    p->gi = 0; p->len = seq->seq.l;
    p->offset = (bns->n_seqs == 0)? 0 : (p-1)->offset + (p-1)->len;
    p->n_ambs = 0;
    for (size_t i = lasts = 0; i < seq->seq.l; ++i) {
      int c = nst_nt4_table[(int)seq->seq.s[i]];
      if (c >= 4) { // N
	if (lasts == seq->seq.s[i]) { // contiguous N
	  ++(*q)->len;
	} else {
	  if (bns->n_holes == *m_holes) {
	    (*m_holes) <<= 1;
	    bns->ambs = (bntamb1_t*)realloc(bns->ambs, (*m_holes) * sizeof(bntamb1_t));
	  }
	  *q = bns->ambs + bns->n_holes;
	  (*q)->len = 1;
	  (*q)->offset = p->offset + i;
	  (*q)->amb = seq->seq.s[i];
	  ++p->n_ambs;
	  ++bns->n_holes;
	}
      }
      lasts = seq->seq.s[i];
      { // fill buffer
	if (c >= 4) c = lrand48()&3;
	if (bns->l_pac == *m_pac) { // double the pac size
	  *m_pac <<= 1;
	  pac = (uint8_t*)realloc(pac, *m_pac/4);
	  memset(pac + bns->l_pac/4, 0, (*m_pac - bns->l_pac)/4);
	}
	_set_pac(pac, bns->l_pac, c);
	++bns->l_pac;
      }
    }
    ++bns->n_seqs;
    
    return pac;
  }
  
  // modified from bwa (heng li)
  uint8_t* BWAIndex::seqlib_make_pac(const UnalignedSequenceVector& v, bool for_only) const
  {
    
    bntseq_t * bns = (bntseq_t*)calloc(1, sizeof(bntseq_t));
    uint8_t *pac = 0; 
    int32_t m_seqs, m_holes;
    int64_t m_pac, l;
    bntamb1_t *q;
    
    bns->seed = 11; // fixed seed for random generator
    m_seqs = m_holes = 8; m_pac = 0x10000;
    bns->anns = (bntann1_t*)calloc(m_seqs, sizeof(bntann1_t));
    bns->ambs = (bntamb1_t*)calloc(m_holes, sizeof(bntamb1_t));
    pac = (uint8_t*) calloc(m_pac/4, 1);
    q = bns->ambs;
    
    // move through the unaligned sequences
    for (size_t k = 0; k < v.size(); ++k) {
      
      // make the ref name kstring
      kstring_t * name = (kstring_t*)malloc(1 * sizeof(kstring_t));
      name->l = v[k].Name.length() + 1;
      name->m = v[k].Name.length() + 3;
      name->s = (char*)calloc(name->m, sizeof(char));
      memcpy(name->s, v[k].Name.c_str(), v[k].Name.length()+1);
      
      // make the sequence kstring
      kstring_t * t = (kstring_t*)malloc(sizeof(kstring_t));
      t->l = v[k].Seq.length();
      t->m = v[k].Seq.length() + 2;
      //t->s = (char*)calloc(v[k].Seq.length(), sizeof(char));
      t->s = (char*)malloc(t->m);
      memcpy(t->s, v[k].Seq.c_str(), v[k].Seq.length());
      
      // put into a kstring
      kseq_t *ks = (kseq_t*)calloc(1, sizeof(kseq_t));  
      ks->seq = *t;
      ks->name = *name;
      
      // make the forward only pac
      pac = seqlib_add1(ks, bns, pac, &m_pac, &m_seqs, &m_holes, &q);
      
      // clear it out
      free(name->s);
      free(name);
      free(t->s);
      free(t);
      //free(ks->name.s); 
      //free(ks->seq.s);
      //free(ks->f->buf);
      //free(
      free(ks);
      // NOTE free kstring_t?
      //kseq_destroy(s);
    }
    
    if (!for_only) 
      {
	// add the reverse complemented sequence
	m_pac = (bns->l_pac * 2 + 3) / 4 * 4;
	pac = (uint8_t*)realloc(pac, m_pac/4);
	memset(pac + (bns->l_pac+3)/4, 0, (m_pac - (bns->l_pac+3)/4*4) / 4);
	for (l = bns->l_pac - 1; l >= 0; --l, ++bns->l_pac)
	  _set_pac(pac, bns->l_pac, 3-_get_pac(pac, l));
      }
    
    bns_destroy(bns);
    
    return pac;
  }
  
  // modified from bwa (heng li)
  bwt_t *BWAIndex::seqlib_bwt_pac2bwt(const uint8_t *pac, int bwt_seq_lenr) const
  {
    
    bwt_t *bwt;
    ubyte_t *buf;
    int i;
    //FILE *fp;
    
    // initialization
    bwt = (bwt_t*)calloc(1, sizeof(bwt_t));
    bwt->seq_len = bwt_seq_lenr; //bwa_seq_len(fn_pac); //dummy
    bwt->bwt_size = (bwt->seq_len + 15) >> 4;
    //fp = xopen(fn_pac, "rb");
    
    // prepare sequence
    //pac_size = (bwt->seq_len>>2) + ((bwt->seq_len&3) == 0? 0 : 1);
    //buf2 = (ubyte_t*)calloc(pac_size, 1);
    //err_fread_noeof(buf2, 1, pac_size, fp);
    //err_fclose(fp);
    memset(bwt->L2, 0, 5 * 4);
    buf = (ubyte_t*)calloc(bwt->seq_len + 1, 1);
    for (i = 0; i < (int)bwt->seq_len; ++i) {
      buf[i] = pac[i>>2] >> ((3 - (i&3)) << 1) & 3;
      ++bwt->L2[1+buf[i]];
    }
    for (i = 2; i <= 4; ++i) 
      bwt->L2[i] += bwt->L2[i-1];
    //free(buf2);
    
    // Burrows-Wheeler Transform
    bwt->primary = is_bwt(buf, bwt->seq_len);
    bwt->bwt = (u_int32_t*)calloc(bwt->bwt_size, 4);
    for (i = 0; i < (int)bwt->seq_len; ++i)
      bwt->bwt[i>>4] |= buf[i] << ((15 - (i&15)) << 1);
    free(buf);
    return bwt;
  }
  
  // modified from bwa (heng li)
  bntann1_t* BWAIndex::seqlib_add_to_anns(const std::string& name, const std::string& seq, bntann1_t* ann, size_t offset) const
  {
    
    ann->offset = offset;
    ann->name = (char*)malloc(name.length()+1); // +1 for \0
    strncpy(ann->name, name.c_str(), name.length()+1);
    ann->anno = (char*)malloc(7);
    strcpy(ann->anno, "(null)\0");
    ann->len = seq.length();
    ann->n_ambs = 0; // number of "holes"
    ann->gi = 0; // gi?
    ann->is_alt = 0;
    
    return ann;
  }
  
  void BWAIndex::seqlib_write_pac_to_file(const std::string& file) const
  {
    // finalize .pac file
    FILE *fp;
    std::string nm = file + ".pac";
    fp = xopen(nm.c_str(), "wb");
    ubyte_t ct;
    err_fwrite(idx_->pac, 1, (idx_->bns->l_pac>>2) + ((idx_->bns->l_pac&3) == 0? 0 : 1), fp);
    
    // the following codes make the pac file size always (l_pac/4+1+1)
    if (idx_->bns->l_pac % 4 == 0) {
      ct = 0;
      err_fwrite(&ct, 1, 1, fp);
    }
    ct = idx_->bns->l_pac % 4;
    err_fwrite(&ct, 1, 1, fp);
    
    // close .pac file
    err_fflush(fp);
    err_fclose(fp);
  }
  
  void BWAIndex::WriteIndex(const std::string& prefix) const {
    if (!idx_) {
      throw std::runtime_error("BWAIndex::writeIndex: no index loaded");
    }
    
    // Prepare filenames
    const std::string bwtFile = prefix + ".bwt";
    const std::string saFile  = prefix + ".sa";
    const std::string bnsFile = prefix + ".bns";
    
    // Dump BWT and SA
    if (std::filesystem::is_directory(std::filesystem::path(prefix))) {
      throw std::runtime_error("BWAIndex::writeIndex: prefix refers to a directory");
    }
    bwt_dump_bwt(bwtFile.c_str(), idx_->bwt);
    bwt_dump_sa (saFile.c_str(),  idx_->bwt);
    
    // Dump the BNS file
    bns_dump(idx_->bns, prefix.c_str());
    
    // Dump the .pac file
    seqlib_write_pac_to_file(prefix);
    
    return;
  }
  
  std::ostream& operator<<(std::ostream& os, const BWAIndex& idx) {
    if (!idx.idx_) {
      os << "[BWAIndex] <no index loaded>";
    } else {
      os << "[BWAIndex] #seqs=" << idx.idx_->bns->n_seqs
	 << " pac_len=" << idx.idx_->bns->l_pac
	 << " holes=" << idx.idx_->bns->n_holes;
    }
    return os;
  }
  
}
