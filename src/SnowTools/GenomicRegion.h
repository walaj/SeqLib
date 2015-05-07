#ifndef SNOWTOOLS_GENOMIC_REGION_H__
#define SNOWTOOLS_GENOMIC_REGION_H__

#include <vector>
#include <iostream>
#include <cstdint>
#include <utility>
#include <list>

#include "SnowTools/SnowUtils.h"
#include "SnowTools/SnowToolsCommon.h"

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
  GenomicRegion(int32_t t_chr, uint32_t t_pos1, uint32_t t_pos2, bool t_strand = true);

  /** Construct a GenomicRegion from a string
   */
  GenomicRegion(std::string t_chr, std::string t_pos1, std::string t_pos2);

  static int32_t chrToNumber(std::string ref);
  static std::string chrToString(int32_t ref);

  static uint32_t posToBigPos(int32_t refid, uint32_t pos);

  /** Randomize the position of this GenomicRegion on the genome
   * 
   * Creates a GenomicRegion with pos1 = pos2. Simulates a random value
   * with val <= genome_size_XY and then converts to GenomicRegion
   */
  void random();

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
  void pad(uint32_t pad);
  
  int width() const;

  int32_t chr = 0;
  uint32_t pos1 = 0;
  uint32_t pos2 = 0;
  bool strand = true; // true = pos, false = neg

 private:


  //char strand = '*';
  //std::string id;
  //int mapq = 0;
  
};

typedef std::vector<GenomicRegion> GenomicRegionVector;
typedef GenomicRegion GR;

}


#endif
