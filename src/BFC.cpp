/*
A significant portion of this code is derived from Heng Li's BFC
repository: https://github.com/lh3/bfc

BFC is copyrighted by Heng Li with the following license:

The MIT License
 
Copyright (c) 2015 Broad Institute
 
Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:
 
The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.
 
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/

#include "SeqLib/BFC.h"

#include <stdexcept>
#include <algorithm>

namespace SeqLib {

  BFC::BFC() {
    m_idx     = 0; //< iterator for getting sequences back when done
    bfc_opt_init(&bfc_opt);
    fml_opt_init(&fml_opt);
    ch        = nullptr;
    kmer      = 0;
    flt_uniq  = 0;
    kcov      = 0;
  }

  BFC::~BFC() {
    ClearReads();
    // clear the old if there
    if (ch) {
      bfc_ch_destroy(ch);
      ch = nullptr;
    }
  }
  
  bool BFC::AddSequence(std::string_view seq,
			std::string_view qual,
			std::string_view name) 
  {
    // refuse empty seq or mismatched qual
    if (seq.empty() || (!qual.empty() && qual.size() != seq.size()))
      return false;
    
    // build a new fseq1_t
    fseq1_t s{};
    s.l_seq = static_cast<uint32_t>(seq.size());
    s.seq   = strndup(seq.data(), seq.size());
    s.qual  = qual.empty() 
      ? nullptr 
      : strndup(qual.data(), qual.size());
    
    // append it
    m_seqs.push_back(s);
    
    // store the name
    m_names.emplace_back(name);
    
    // sanity check
    assert(m_names.size() == m_seqs.size());
    
    return true;
  }
  
  
  // bool BFC::AddSequence(std::string_view seq,
  // 			std::string_view qual,
  // 			std::string_view name) {
    
  //   // do the intial allocation
  //   if (n_seqs == 0 && !m_seqs) {
  //     m_seqs_size = 32;
  //     m_seqs = (fseq1_t*)malloc(m_seqs_size * sizeof(fseq1_t));
  //   }
  //   // realloc if not enough space
  //   else if (n_seqs >= m_seqs_size) {
  //     m_seqs_size = 2 * m_seqs_size;
  //     m_seqs = (fseq1_t*)realloc(m_seqs, m_seqs_size * sizeof(fseq1_t));
  //   }
    
  //   if (!m_seqs)
  //     return false;
    
  //   // make sure seq and qual are even valid (if qual provided)
  //   if (!qual.empty() && !seq.empty())
  //     if (seq.length() != qual.length())
  // 	return false;
  //   //    if (strlen(qual) && seq && qual) 
  //   //if (strlen(seq) != strlen(qual))
  //   //return false;
  //   if (seq.empty())
  //     return false;
    
  //     //      if (!strlen(seq))
  //     //return false;

  //   fseq1_t *s;
    
  //   s = &m_seqs[n_seqs];
    
  //   s->seq   = strndup(seq.data(), seq.size()); //seq.c_str());
  //   s->qual = 0;
  //   if (!qual.empty()) {
  //     //if (strlen(qual)) {
  //     s->qual  = strndup(qual.data(), qual.size()); //qual.c_str());
  //   }
    
  //   s->l_seq = seq.length(); //strlen(seq);
  //   n_seqs++;

  //   m_names.push_back(std::string(name));
  //   //m_names.push_back(name); //strdup(name));
    
  //   assert(m_names.size() == n_seqs);
    
  //   return true;
  // }

  bool BFC::GetSequence(std::string& s, std::string& q) {
    //if (m_idx >= n_seqs)
    assert(m_seqs.size() == m_names.size());
    if (m_idx >= m_seqs.size())
      return false;
    s = std::string(m_seqs.at(m_idx).seq);
    q = m_names.at(m_idx);
    std::transform(s.begin(), s.end(), s.begin(), ::toupper);
    ++m_idx;
    return true;
  }

  // correct a single sequence
  /*  bool BFC::CorrectSequence(std::string& str, const std::string& q) {

    assert(n_seqs == 0);
    assert(m_names.size() == 0);

    m_seqs = (fseq1_t*)malloc(1 * sizeof(fseq1_t));
    n_seqs = 1;
    
    //uint64_t size = 0;
    //for (std::vector<char*>::const_iterator r = v.begin(); r != v.end(); ++r) {
    //    for (auto& r : v) {
    fseq1_t *s;
    s = &m_seqs[0];
    s->seq   = strdup(str.c_str());
    s->qual  = q.empty() || q.length() != str.length() ? NULL : strdup(q.c_str()); 
    s->l_seq = str.length();

    // do the error correction of this one sequence
    ErrorCorrect();

    // add a dummy name
    //m_names.push_back(strdup("1"));
    m_names.push_back("1");

    // send to uppercase, and return 
    std::string cstr = std::string(m_seqs[0].seq);
    std::transform(cstr.begin(), cstr.end(), cstr.begin(), ::toupper);
    str = cstr;

    clear();

    return true;

    }*/

  void free_char(char*& c) {
    if (c) {
      free (c);
      c = NULL;
    }
  }

  void BFC::ClearReads() {
    assert(m_names.size() == m_seqs.size());
    for (size_t i = 0; i < m_seqs.size(); ++i) {
      free_char(m_seqs[i].seq);
      free_char(m_seqs[i].qual);
    }
    m_seqs.clear();
    m_names.clear();
    m_idx = 0;
  }
  

  void BFC::Train() {

    // clear the old if there
    if (ch) {
      bfc_ch_destroy(ch);
      ch = nullptr;
    }
    
    // Initialize Fermi-lite options
    fml_opt_init(&fml_opt);

    // Auto-learn kmer size if not set
    if (kmer <= 0) {
      fml_opt_adjust(&fml_opt, m_seqs.size(), m_seqs.data());
      kmer = fml_opt.ec_k;
    }

    // initialize BFC options
    //////
    
    // Compute total sequence length
    uint64_t tot_len = 0;
    for (auto const& s : m_seqs) tot_len += s.l_seq;
    // l_pre = min(tot_len - 8, 20)
    bfc_opt.l_pre = (tot_len > 8)
      ? static_cast<int>(std::min<uint64_t>(tot_len - 8, 20))
      : 0;
    
    //OLD
    //for (size_t i = 0; i < n_seqs; ++i) 
    //  tot_len += m_seqs[i].l_seq; // compute total length
    //bfc_opt.l_pre = tot_len - 8 < 20? tot_len - 8 : 20;
    
    //  setup the counting of kmers
    // Reset error-correction state
    std::memset(&es, 0, sizeof(es));    
    //OLD //memset(&es, 0, sizeof(ec_step_t));
    //kmer is learned before this
    
    bfc_opt.k = kmer;
    
    //es.opt = &bfc_opt, es.n_seqs = n_seqs, es.seqs = m_seqs, es.flt_uniq = flt_uniq;
    
    // hold count info. also called bfc_ch_s. Composed of
    //    int k
    //    int l_pre
    //    cnthash_t **h
    //        h is of size 1<<l_pre (2^l_pre). It is array of hash tables
    //        h[i] is initialized with kh_init(cnt) which makes a cnthash_t
    // bfc_ch_t *ch; // set in BFC.h
    
    // Perform k-mer counting
    ch = fml_count(
		   m_seqs.size(),
		   m_seqs.data(),
      bfc_opt.k,
      bfc_opt.q,
      bfc_opt.l_pre,
      bfc_opt.n_threads
    );
    //ch = fml_count(n_seqs, m_seqs, bfc_opt.k, bfc_opt.q, bfc_opt.l_pre, bfc_opt.n_threads);

#ifdef DEBUG_BFC
    // size of random hash value
    khint_t k;
    int* ksize = (int*)calloc(1<<ch->l_pre, sizeof(int));
    for (int i = 0; i < (1<<ch->l_pre); ++i) {
      for (k = kh_begin(ch->h[i]); k != kh_end(ch->h[i]); ++k)
        ++ksize[i];
      fprintf(stderr, "K: %d S: %d\n", i, ksize[i]);
    }
#endif
  }

  void BFC::ErrorCorrect() {
    
    assert(kmer > 0);

    // 1) Prepare the ec_step_t
    es.ch       = ch;
    es.opt      = &bfc_opt;
    es.n_seqs   = static_cast<int>(m_seqs.size());
    es.seqs     = m_seqs.data();
    es.flt_uniq = flt_uniq;

    //OLD
    //es.ch = ch;
    //es.opt = &bfc_opt;
    //es.n_seqs = n_seqs;
    //es.seqs = m_seqs;
    //es.flt_uniq = flt_uniq;

    // 2) Build the k-mer histogram
    // histogram of kmer occurences
    uint64_t hist[256];
    // diff histogram of kmers??
    uint64_t hist_high[64];
    std::fill(std::begin(hist), std::end(hist), 0);
    std::fill(std::begin(hist_high), std::end(hist_high), 0);
    
    int mode = bfc_ch_hist(es.ch, hist, hist_high);
    
    // make the histogram?
    // es.ch is unchanged (const)
    //OLD//int mode = bfc_ch_hist(es.ch, hist, hist_high);
    
    // 3) Accumulate sum_k and tot_k
    uint64_t sum_k = 0; ///< total valid kmer count
    uint64_t tot_k = 0; ///< total number of kmers
    for (int i = fml_opt.min_cnt; i < 256; ++i) {
      sum_k += hist[i];
      tot_k  += i * hist[i];
    }

    //OLD
    //for (int i = fml_opt.min_cnt; i < 256; ++i) 
    //  sum_k += hist[i], tot_k += i * hist[i];    

#ifdef DEBUG_BFC
    std::cerr << " sum_k " << sum_k << " tot_k " << tot_k << std::endl;
    fprintf(stderr, "MODE: %d\n", mode);
    for (int i = fml_opt.min_cnt; i < 256; ++i) {
      fprintf(stderr, "hist[%d]: %d\n",i,hist[i]);
    }
    for (int i = fml_opt.min_cnt; i < 64; ++i) {
      fprintf(stderr, "hist_high[%d]: %d\n",i,hist_high[i]);
    }
#endif

    // 4) Estimate coverage & set min_cov
    kcov = sum_k ? static_cast<float>(tot_k) / sum_k : 0.0f;
    int raw_min = static_cast<int>(BFC_EC_MIN_COV_COEF * kcov + 0.499f);
    bfc_opt.min_cov = std::clamp(raw_min, fml_opt.min_cnt, fml_opt.max_cnt);
    
    // 5) Run the actual correction pass
    kmer_correct(&es, mode, ch);

    //OLD
    //kcov = (float)tot_k / sum_k;
    //bfc_opt.min_cov = (int)(BFC_EC_MIN_COV_COEF * kcov + .499);
    //bfc_opt.min_cov = bfc_opt.min_cov < fml_opt.max_cnt? bfc_opt.min_cov : fml_opt.max_cnt;
    //bfc_opt.min_cov = bfc_opt.min_cov > fml_opt.min_cnt? bfc_opt.min_cov : fml_opt.min_cnt;

#ifdef DEBUG_BFC
    fprintf(stderr, "kcov: %f mincov: %d  mode %d \n", kcov, bfc_opt.min_cov, mode);  
#endif

    //OLD
    // do the actual error correction
    //kmer_correct(&es, mode, ch);

    return;


  }

}
