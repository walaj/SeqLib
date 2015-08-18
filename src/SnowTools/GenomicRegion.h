#ifndef SNOWTOOLS_GENOMIC_REGION_H__
#define SNOWTOOLS_GENOMIC_REGION_H__

#include <vector>
#include <iostream>
#include <cstdint>
#include <utility>
#include <list>
#include <cstring>
#include <memory>

#include "SnowTools/SnowUtils.h"
#include "SnowTools/SnowToolsCommon.h"

#include "htslib/sam.h"

/** 
 */
namespace SnowTools {

class GenomicRegion {

  template<typename T> friend class GenomicRegionCollection;
  
 public:

  /** Construct an "empty" GenomicRegion at (chr -1), pos 0, width = 1
   */
  GenomicRegion() : chr(-1), pos1(0), pos2(0) {};

  /** Construct a GenomicRegion at a specific start and end location 
   * @param t_chr Chromosome id  (chr1 = 0, etc)
   * @param t_pos1 Start position
   * @param t_pos2 End position
   * @param strand true for positive, false for negative
  */
  GenomicRegion(int32_t t_chr, uint32_t t_pos1, uint32_t t_pos2, char t_strand = true);

  /** Construct a GenomicRegion from a string
   */
  //GenomicRegion(std::string t_chr, std::string t_pos1, std::string t_pos2);
  
  /** Construct a GenomicRegion from a set of strings */
  GenomicRegion(const std::string& tchr, const std::string& tpos1, const std::string& tpos2, bam_hdr_t* h = NULL);

  /** Construct a GenomicRegion from a samtools style region string.
   *
   * This calls the samtools-like parser, which accepts in form "chr7:10,000-11,100".
   * Note that this requires that a pointer to the BAM header be provided as well 
   * to convert the text representation of the chr to the id number.
   */
  GenomicRegion(const std::string& reg, bam_hdr_t* h);

  /** Return a string representation of just the first base-pair 
   */
  std::string pointString() const;

  static int32_t chrToNumber(std::string ref);

  static std::string chrToString(int32_t ref);

  /** Randomize the position of this GenomicRegion on the genome
   * 
   * Creates a GenomicRegion with pos1 = pos2. Simulates a random value
   * with val <= genome_size_XY and then converts to GenomicRegion
   */
  void random(uint32_t seed = 0);

  /** Does this GenomicRegion represent a valid region? */
  bool valid() const { return chr >= 0; }

  /** Check if the GenomicRegion is empty (aka chr -1)
   */
  bool isEmpty() const;

  /** Find the distance between two GenomicRegion objects
   */
  int distance(const GenomicRegion &gr) const;

  // define how these are to be sorted
  bool operator < (const GenomicRegion& b) const;
  bool operator==(const GenomicRegion& b) const;
  bool operator<=(const GenomicRegion &b) const;
  
  // determine if something overlaps with centromere 
  int centromereOverlap() const;

  // check if there is an overlap
  int getOverlap(const GenomicRegion gr) const;

  friend std::ostream& operator<<(std::ostream& out, const GenomicRegion& gr);

  std::string toString() const;

  std::string ChrName(const bam_hdr_t* h) const;

  void pad(int32_t pad);
  
  int width() const;

  int32_t chr = 0;
  int32_t pos1 = 0;
  int32_t pos2 = 0;
  char strand = '*';

  //  bam_hdr_t * m_hdr = nullptr;

 private:

  //char strand = '*';
  //std::string id;
  //int mapq = 0;
  
};

typedef std::vector<GenomicRegion> GenomicRegionVector;
typedef GenomicRegion GR;

}


#endif
