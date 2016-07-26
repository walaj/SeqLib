#ifndef SNOWTOOLS_GENOMIC_REGION_H__
#define SNOWTOOLS_GENOMIC_REGION_H__

#include <vector>
#include <iostream>
#include <cstdint>
#include <utility>
#include <list>
#include <cstring>
#include <memory>

#include "SnowTools/SnowToolsCommon.h"
#include "SnowTools/SnowUtils.h"

#include "htslib/sam.h"

/** 
 */
namespace SnowTools {

  /** @brief Container for an interval on the genome 
   */
class GenomicRegion {

  template<typename T> friend class GenomicRegionCollection;
  
 public:

  /** Construct an "empty" GenomicRegion at (chr -1), pos 0, width = 1
   */
  GenomicRegion() : chr(-1), pos1(0), pos2(0) {};

  /** Construct a GenomicRegion at a specific start and end location 
   * @param t_chr Chromosome id  (chr1 = 0, etc)
   * @param t_pos1 Start position
   * @param t_pos2 End position. Must be >= start position.
   * @param strand. +, -, or * (default is *)
   * @exception throws an invalid_argument exception if pos2 < pos1
   * @exception throws an invalid_argument exception if char not one of +, - , *
  */
  GenomicRegion(int32_t t_chr, int32_t t_pos1, int32_t t_pos2, char t_strand = '*');

  /** Construct a GenomicRegion from a set of strings */
  GenomicRegion(const std::string& tchr, const std::string& tpos1, const std::string& tpos2, bam_hdr_t* h = NULL);

  /** Construct a GenomicRegion from a samtools style region string.
   *
   * This calls the samtools-like parser, which accepts in form "chr7:10,000-11,100".
   * Note that this requires that a pointer to the BAM header be provided as well 
   * to convert the text representation of the chr to the id number.
   * @param reg Samtools-style string (e.g. "1:1,000,000-2,000,000")
   * @param h Pointer to BAM header that will be used to convert chr string to ref id
   * @exception throws an invalid_argument exception if cannot parse correctly
   */
  GenomicRegion(const std::string& reg, bam_hdr_t* h);

  /** Return a string representation of just the first base-pair 
   * e.g. 1:10,000
   */
  std::string pointString() const;

  /** Convert a chromosome number to a string using default ordering (1-Y)
   * Assumes a 1-based ordering (1, ...), not zero-based.
   * e.g. chrToString(10) return "11"
   * @param ref Reference ID to convert
   * @exception throws an invalid_argument exception if ref < 0
   */
  static std::string chrToString(int32_t ref);

  /** Randomize the position of this GenomicRegion on the genome
   * 
   * Creates a GenomicRegion with pos1 = pos2. Simulates a random value
   * with val <= genome_size_XY and then converts to GenomicRegion
   * @note Seed is set before-hand at any time with srand
   */
  void random();

  /** Returns true if chr id >= 0, false otherwise
   */
  bool valid() const { return chr >= 0; }

  /** Check if the GenomicRegion is empty (aka chr -1 and pos1=pos2=0)   */
  bool isEmpty() const;

  /** Find the absolute distance between start of two GenomicRegion objects 
   * 
   * If chr1 != chr2, then -1 is returned
   * @param gr GenomicRegion object to compare with
   */
  int32_t distanceBetweenStarts(const GenomicRegion &gr) const;

  /** Find the absolute distance between ends of two GenomicRegion objects 
   * 
   * If chr1 != chr2, then -1 is returned
   * @param gr GenomicRegion object to compare with
   */
  int32_t distanceBetweenEnds(const GenomicRegion &gr) const;

  /** Output as a string, with chr ID bumped up by one to make ID 0
   * print as "1". Add commas to pos
   */
  std::string toPrettyString() const;

  // define how these are to be sorted
  bool operator < (const GenomicRegion& b) const;
  bool operator==(const GenomicRegion& b) const;
  bool operator<=(const GenomicRegion &b) const;
  
  /** Check if the GenomicRegion has a complete or partial overlap
   * If the argument contains the calling object, returns 3
   * If the argument is contained in the calling object, returns 2
   * If the argument overlaps partially the calling object, returns 1
   * If the argument and calling object do not overlap, returns 0
   * @param gr GenomicRegion to compare against
   */
  int getOverlap(const GenomicRegion& gr) const;

  friend std::ostream& operator<<(std::ostream& out, const GenomicRegion& gr);

  std::string toString() const;

  /** Extract the chromosome name as a string 
   * @param h BAM header with h->target_name field
   * @exception throws an invalid_argument exception if ref id >= h->n_targets
   */
  std::string ChrName(const bam_hdr_t* h = nullptr) const;

  /** Pad the object to make larger or smaller
   * @param pad Amount to pad by.
   * @exception throws an out_of_bounds if for pad < -width/2
   */
  void pad(int32_t pad);
  
  int width() const;

  int32_t chr = 0;
  //int32_t chr:30, strand:2;
  int32_t pos1 = 0;
  int32_t pos2 = 0;
  char strand = '*';

 private:

};

typedef std::vector<GenomicRegion> GenomicRegionVector;

}


#endif
