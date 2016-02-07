#ifndef SNOWTOOLS_SVVCF_H__
#define SNOWTOOLS_SVVCF_H__

#include <string>
#include <unordered_map>

#include "SnowTools/GenomicRegionCollection.h"
#include "SnowTools/RefGenome.h"
#include "SnowTools/BamRead.h"

namespace SnowTools {

  /** Represent a VCF containing structural variants.
   *
   * [In development]
   */
  class SVVCF {

  public:

    SVVCF() {}

    SVVCF(const std::string& vcf_file, bam_hdr_t* h = nullptr);
    
    size_t size() const;

    std::string getSurroundingContigFromReference(size_t i, RefGenome& ref, const bam_hdr_t *h, BamReadVector& brv, int width = 100) const;

  private:

    GRC m_side1, m_side2;
    std::unordered_map<std::string, int> m_hash;
    std::vector<std::string> m_insertions, m_homology;
    std::vector<std::string> m_name;
  };


};


#endif
