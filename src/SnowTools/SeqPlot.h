#ifndef SNOWTOOLS_CONTIG_PLOT_H__
#define SNOWTOOLS_CONTIG_PLOT_H__

#include "SnowTools/BamRead.h"

namespace SnowTools {
  
  /** Object for creating ASCII alignment plots
   */

  class ReadPlot {

  public:
    
    /** Create an empty plot */
    ReadPlot() {}

  private: 

    // reads that align to the contig
    std::vector<BamRead> m_reads;

  };

}



#endif
