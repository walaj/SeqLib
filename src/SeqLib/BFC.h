#ifndef SEQLIB_BFC_H__
#define SEQLIB_BFC_H__

extern "C" {
  #include "fermi-lite/bfc.h"
  #include "fermi-lite/fml.h"
}

#include "SeqLib/BamRecord.h"

namespace SeqLib {

/** Class to perform error-correction using BFC algorithm
 *
 * BFC is designed and implemented by Heng Li (https://github.com/lh3/bfc). 
 * From Heng: It is a variant of the classical spectrum alignment algorithm introduced
 * by Pevzner et al (2001). It uses an exhaustive search to find a k-mer path 
 * through a read that minimizeds a heuristic objective function jointly considering
 * penalities on correction, quality and k-mer support.
 */
  class BFC {

  public:
    /** Construct a new BFC engine */
    BFC() {
      bfc_opt_init(&m_opt);
    }

    ~BFC() {
      //if (m_opt)
      //free(m_opt);
    }

    /** NOT IMPLEMENTED YET */
    std::vector<std::pair<std::string,std::string>> ErrorCorrect(BamRecordVector& brv);

  private:

    bfc_opt_t m_opt;

    // histogram of kmer occurences
    uint64_t hist[256];

    // diff histogram of kmers??
    uint64_t hist_high[64];

    uint64_t tot_len = 0;

    uint64_t sum_k = 0; // total valid kmer count (kmers above min_count) ? 

    // total number of kmers?
    uint64_t tot_k = 0;

    float kcov = 0;

  };

   }

#endif
