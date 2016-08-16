#include "SeqLib/FermiAssembler.h"
#include "fermi-lite/bfc.h"

namespace SeqLib {

  FermiAssembler::FermiAssembler() {
    fml_opt_init(&opt);
  }
  
  FermiAssembler::~FermiAssembler() {
    ClearReads();
    
    //if (opt.mag_opt)
    //  free(opt.mag_opt);
  }

  void FermiAssembler::AddReads(const BamRecordVector& brv) {
    //m_brv.insert(m_brv.end(), brv.begin(), brv.end());

    // alloc the memory
    m_seqs = (fseq1_t*)realloc(m_seqs, (n_seqs + brv.size()) * sizeof(fseq1_t));

    int m = 0;
    uint64_t size = 0;
    for (auto& r : brv) {
      fseq1_t *s;
      //if (n_seqs >= m) {
      //m = m? m<<1 : 256;
      //m_seqs = (fseq1_t*)realloc(m_seqs, m * sizeof(fseq1_t));
      //}

      s = &m_seqs[n_seqs];

      s->seq   = strdup(r.Sequence().c_str());
      s->qual  = strdup(r.Qualities().c_str());

      s->l_seq = r.Sequence().length();
      size += m_seqs[n_seqs++].l_seq;
    }
    
  }

  void FermiAssembler::ClearReads() {  
    for (size_t i = 0; i < n_seqs; ++i) {
      fseq1_t * s = &m_seqs[i];
      if (s->qual)
       free(s->qual); 
      if (s->seq)
	free(s->seq);
    }
    free(m_seqs);
    m_seqs = nullptr;
      
  }

  void FermiAssembler::CorrectReads() {  
    fml_correct(&opt, n_seqs, m_seqs);
  }
  
}
