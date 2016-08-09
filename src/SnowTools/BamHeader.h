#ifndef SNOWTOOLS_BAM_HEADER_H__
#define SNOWTOOLS_BAM_HEADER_H__

#include "htslib/hts.h"
#include "htslib/sam.h"
#include "htslib/bgzf.h"
#include "htslib/kstring.h"

#include <string>

namespace SnowTools {
  
  /** Store a header to a BAM file 
   *
   * Stores a BAM header, which also acts as a dictionary of 
   * reference sequences, with names and lengths.
   */
  class BamHeader {

  public:

    /** Initializes a new empty BamHeader with no data
     * 
     * @note No memory is allocated here
     */
    BamHeader() { h = nullptr; };
    
    /** Initialize a BamHeader from a string containing
     * a BAM header in human-readable form (e.g. @PG ... )
     * @param Text of a BAM header, with newlines separating lines
     */
    BamHeader(const std::string& hdr);

  private:

    bam_hdr_t * h;
    

  };
  
}


#endif
