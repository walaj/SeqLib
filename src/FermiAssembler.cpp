#include "SeqLib/FermiAssembler.h"

namespace SeqLib {

  FermiAssembler::FermiAssembler() {
    fml_opt_init(&opt);
  }
  
  FermiAssembler::~FermiAssembler() {
    ClearReads();
    ClearContigs();
    //if (opt.mag_opt)
    //  free(opt.mag_opt);
  }
  
  void FermiAssembler::AddReads(const BamRecordVector& brv) {

    // alloc the memory
    m_seqs = (fseq1_t*)realloc(m_seqs, (n_seqs + brv.size()) * sizeof(fseq1_t));

    int m = 0;
    uint64_t size = 0;
    for (auto& r : brv) {
      m_names.push_back(r.Qname());
      fseq1_t *s;

      s = &m_seqs[n_seqs];

      s->seq   = strdup(r.Sequence().c_str());
      s->qual  = strdup(r.Qualities().c_str());

      s->l_seq = r.Sequence().length();
      size += m_seqs[n_seqs++].l_seq;
    }
    
  }

  void FermiAssembler::ClearContigs() {
    fml_utg_destroy(n_utgs, m_utgs);  
    m_utgs = 0;
    n_utgs = 0;
  }

  void FermiAssembler::ClearReads() {  
    if (!m_seqs)
      return; //already cleared

    for (size_t i = 0; i < n_seqs; ++i) {
      fseq1_t * s = &m_seqs[i];
      if (s->qual)
       free(s->qual); 
      s->qual = nullptr;
      if (s->seq)
	free(s->seq);
      s->seq = nullptr;
    }
    free(m_seqs);
    m_seqs = nullptr;
      
  }

  void FermiAssembler::CorrectReads() {  
    fml_correct(&opt, n_seqs, m_seqs);
  }

  void FermiAssembler::CorrectAndFilterReads() {  
    fml_fltuniq(&opt, n_seqs, m_seqs);
  }

  void FermiAssembler::PerformAssembly() {
    m_utgs = fml_assemble(&opt, n_seqs, m_seqs, &n_utgs); // assemble!
  }
  
  std::vector<std::string> FermiAssembler::GetContigs() const {
    std::vector<std::string> c;
    for (size_t i = 0; i < n_utgs; ++i)
      c.push_back(std::string(m_utgs[i].seq));
    return c;
  }

  /*void FermiAssembler::count() {

    // initialize BFC options
    uint64_t tot_len = 0;
    for (int i = 0; i < n_seqs; ++i) 
      tot_len += m_seqs[i].l_seq; // compute total length
    int l_pre = tot_len - 8 < 20? tot_len - 8 : 20;

    //bfc_ch_t *ch = fml_count(n_seqs, m_seqs, opt.ec_k, 20, l_pre, opt.n_threads);
    //std::cerr << " ch->k " << ch->k << " ch->l_pre " << ch->l_pre << std::endl;

    // directly from fml count
    cnt_step_t cs;
    cs.n_seqs = n_seqs, cs.seqs = m_seqs, cs.k = opt.ec_k, cs.q = 20;
    cs.ch = bfc_ch_init(cs.k, l_pre);
    }*/

  UnalignedSequenceVector FermiAssembler::GetSequences() const {
    
    UnalignedSequenceVector r;
    for (size_t i = 0; i < n_seqs; ++i) {
      fseq1_t * s = &m_seqs[i];
      UnalignedSequence read;
      if (s->seq)
	read.Seq = (std::string(s->seq));
      read.Name = m_names[i];
      r.push_back(read);
    }
    return r;
  }
  
}
