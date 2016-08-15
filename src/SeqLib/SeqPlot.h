#ifndef SEQLIB_CONTIG_PLOT_H__
#define SEQLIB_CONTIG_PLOT_H__

#include "SeqLib/BamRecord.h"

namespace SeqLib {
  
  /** Object for creating ASCII alignment plots
   */

  class ReadPlot {

  public:
    
    /** Create an empty plot */
    ReadPlot() {}

  private: 

    // reads that align to the contig
    BamRecordVector m_reads;

  };

}



#endif
