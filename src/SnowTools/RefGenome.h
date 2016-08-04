#ifndef SNOWTOOLS_REF_GENOME_H__
#define SNOWTOOLS_REF_GENOME_H__

#include <string>
#include <cstdlib>
#include <iostream>

#include "htslib/faidx.h"


namespace SnowTools {
  
  /** Stores an indexed reference genome
   *
   * RefGenome is currently used as an interface to obtain
   * sequences from the reference given an interval.
   */
  class RefGenome {

  public:

    /** Load a reference genome
     * @param file Path to an indexed reference genome
     */
    RefGenome(const std::string& file);

    /** Create an empty RefGenome object */
    RefGenome() { index = nullptr; }
    
    /** Destroy malloc'ed faidx_t index inside object
     */
    ~RefGenome() { if (index) fai_destroy(index); }
    
    /** Query a region to get the sequence
     * @param chr_name name of the chr to query
     * @param p1 position 1
     * @param p2 position 2
     */
    std::string queryRegion(const std::string& chr_name, int32_t p1, int32_t p2);

    /** Load an indexed reference sequence */
    void retrieveIndex(const std::string& file);
    
    /** Check if reference has been loaded 
     */
    bool empty() const { return (index == nullptr); }
    
  private:

    faidx_t * index;

  };
  

}

#endif
