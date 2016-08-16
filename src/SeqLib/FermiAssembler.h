#ifndef SEQLIB_FERMI_H__
#define SEQLIB_FERMI_H__

#include <string>
#include <cstdlib>
#include <iostream>

#include "SeqLib/BamRecord.h"
#include "fermi-lite/fml.h"

namespace SeqLib {
  
  /** Stores an indexed reference genome
   *
   * RefGenome is currently used as an interface to obtain
   * sequences from the reference given an interval.
   */
  class FermiAssembler {

  public:

    FermiAssembler ();
    
    ~FermiAssembler();

    void AddReads(const BamRecordVector& brv);

    void ClearReads();

    void CorrectReads();

  private:

    // reads to assemble
    fseq1_t *m_seqs = 0;

    // number of reads
    size_t n_seqs = 0;
    
    // options
    fml_opt_t opt;
    

    
  };
  

}

#endif
