#include "SeqLib/BFC.h"

namespace SeqLib {

  std::vector<std::pair<std::string,std::string>> BFC::ErrorCorrect(BamRecordVector& brv) {

    // reads to assemble
    fseq1_t *m_seqs = 0;
    
    // number of reads
    size_t n_seqs = 0;

    // alloc the memory
    m_seqs = (fseq1_t*)realloc(m_seqs, (n_seqs + brv.size()) * sizeof(fseq1_t));

    // options
    fml_opt_t opt;
    fml_opt_init(&opt);

    std::vector<std::string> m_names;

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
  
  
    // 
    int flt_uniq = 0; // from fml_correct call

  // initialize BFC options
  for (int i = 0; i < n_seqs; ++i) 
    tot_len += m_seqs[i].l_seq; // compute total length
  int l_pre = tot_len - 8 < 20? tot_len - 8 : 20;
  
  //  setup the counting of kmers
  ec_step_t es;
  memset(&es, 0, sizeof(ec_step_t));
  fprintf(stderr, "N: %d  K: %d  q: %d  L: %d NT: %d\n", n_seqs, m_opt.k, m_opt.q, l_pre, m_opt.n_threads);
  es.opt = &m_opt, es.n_seqs = n_seqs, es.seqs = m_seqs, es.flt_uniq = flt_uniq;

  bfc_ch_t *ch;

  // do the counting
  es.ch = ch = fml_count(n_seqs, m_seqs, 33, m_opt.q, /*m_opt.*/l_pre, m_opt.n_threads);

  // make the histogram?
  int mode = bfc_ch_hist(es.ch, hist, hist_high);

  for (int i = opt.min_cnt; i < 256; ++i) {
    sum_k += hist[i], tot_k += i * hist[i];
    //std::cerr << "hist[" << i << "]: " << hist[i] << std::endl;
  }
  //for (int i = opt.min_cnt; i < 64; ++i) {
  //  std::cerr << "hist_high[" << i << "]: " << hist_high[i] << std::endl;
  // }

  //std::cerr << " N_SEQS " << n_seqs << " sum_k " << sum_k << " MIN CNT " << opt.min_cnt << std::endl;

  kcov = (float)tot_k / sum_k;
  m_opt.min_cov = (int)(BFC_EC_MIN_COV_COEF * kcov + .499);
  m_opt.min_cov = m_opt.min_cov < opt.max_cnt? m_opt.min_cov : opt.max_cnt;
  m_opt.min_cov = m_opt.min_cov > opt.min_cnt? m_opt.min_cov : opt.min_cnt;
  
  std::cerr << " M OPT MIN COV " << m_opt.min_cov << " KCOV " << kcov << std::endl;

  // do the actual error correction
  kmer_correct(&es, mode,ch);

  bfc_ch_destroy(ch);

  std::vector<std::pair<std::string, std::string>> corrected;
  for (int i = 0; i < n_seqs; ++i) {
    corrected.push_back(std::pair<std::string, std::string>(brv[i].Qname(), std::string(m_seqs[i].seq)));
  }

  return corrected;
}

}
