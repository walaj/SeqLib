#ifndef SNOWTOOLS_BREAKPOINT_H__
#define SNOWTOOLS_BREAKPOINT_H__

#include <cstdlib>
#include <string>
#include <unordered_map>

#include "SnowTools/GenomicRegion.h"
#include "SnowTools/GenomicRegionCollection.h"
#include "SnowTools/DiscordantCluster.h"
#include "SnowTools/BamRead.h"

//#define GET_COVERAGE 4

#include "SnowTools/STCoverage.h"

namespace SnowTools {

struct BreakPoint;
typedef std::vector<BreakPoint> BPVec;
typedef std::unordered_map<std::string, size_t> PON;

/**
 *
 */
void runRefilterBreakpoints(int argc, char** argv);

/**
 *
 */
void parseBreakOptions(int argc, char** argv);
 
struct BreakPoint {

  static std::string header() { 
    return "chr1\tpos1\tstrand1\tchr2\tpos2\tstrand2\tspan\tmapq1\tmapq2\tnsplit\ttsplit\tndisc\ttdisc\tncigar\ttcigar\thomology\tinsertion\tcontig\tnumalign\tconfidence\tevidence\tpon_samples\trepeat_seq\tnormal_cov\ttumor_cov\tnormal_allelic_fraction\ttumor_allelic_fraction\tblacklist\treads"; 
  }
  
  // reads spanning this breakpoint
  BamReadVector reads;
  
  // discordant reads supporting this aseembly bp
  DiscordantCluster dc;
  
  // breakpoints on the reference
  GenomicRegion gr1;
  GenomicRegion gr2;
  
  // string of read names concatenated together
  std::string read_names;

  // total coverage at that position
  size_t tcov = 0;
  size_t ncov = 0;

  // total coverage supporting the variant at that position
  size_t tcov_support = 0;
  size_t ncov_support = 0;

  int mapq1 = 0;
  int mapq2 = 0;
  
  int cpos1 = 0;  
  int cpos2 = 0;
  
  int nm1 = 0;
  int nm2 = 0;

  std::string seq;

  std::string cname;

  std::string insertion = "";
  std::string homology = "";

  std::string repeat_seq = "";

  std::string id1;
  std::string id2;
  int matchlen1 = 0;
  int matchlen2 = 0;

  size_t tcigar = 0;
  size_t ncigar = 0;

  size_t min_end_align_length = 0; // minimum length of alignment on end. real short are not to be trusted
  
  //char strand1;
  //char strand2;

  bool isSomatic = false;
  bool isGermline = false;

  size_t pon = 0;

  //unsigned mapq1; 
  //unsigned mapq2; 

  unsigned tsplit1 = 0;
  unsigned tsplit2 = 0;

  unsigned nsplit1 = 0;
  unsigned nsplit2 = 0;

  size_t nsplit = 0;
  size_t tsplit = 0;

  unsigned tall = 0;
  unsigned nall = 0; 

  //int nm1 = 0;
  //int nm2 = 0;

  unsigned num_dups = 0;
   
  //Window window;
  GenomicRegion window;

  //int span;

  unsigned num_align = 0;

  bool part_of_local = false;

  bool local1 = false;
  bool local2 = false;

  std::string evidence = "";
  std::string confidence = "";

  bool isindel = false;

  bool blacklist = false;

  /** Construct a breakpoint from a cluster of discordant reads
   */
  BreakPoint(const DiscordantCluster& tdc);

  
  BreakPoint() {
    gr1.pos1 = 0;
    gr2.pos1 = 0;
    gr1.chr = 0;
    gr2.chr = 0;
    gr1.pos2 = 0;
    gr2.pos2 = 0;
    
  }
  
  BreakPoint(std::string &line);

  /*! Return a string with information useful for printing at the 
   * command line as Snowman runs 
   * @return string with minimal information about the BreakPoint.
   */
  std::string toPrintString() const;
  
  static void readPON(std::string &file, std::unique_ptr<PON> &pmap);

  void __combine_with_discordant_cluster(DiscordantClusterMap& dmap);
  
  /*! @function determine if the breakpoint has split read support
   * @param reference to a vector of read smart pointers that have been aligned to a contig
   * @discussion Note: will cause an error if the AL tag not filled in for the reads. 
   * The AL tag is filled in by AlignedContig::alignReadsToContigs.
   */
  void splitCoverage(BamReadVector &bav);

  /*! Determines if the BreakPoint overlays a blacklisted region. If 
   * and overlap is found, sets the blacklist bool to true.
   *
   * Note that currently this only is set for the pos1 of indels.
   * If the BreakPoint object is not an indel, no action is taken. 
   * @param grm An interval tree map created from a BED file containing blacklist regions
   */
  void checkBlacklist(GRC &grv);

  /*! Compute the allelic fraction (tumor and normal) for this BreakPoint.
   *
   * The allelic fraction is computed by taking the base-pair level coverage
   * as the denominator, and the max of number of split reads and number of 
   * cigar supporting reads as the numerator. Note that because, theoretically
   * but rarely, the number of split reads could be > 0 while the bp-level coverage
   * at a variant could be exactly zero. This is because unmapped reads could be called split
   * reads but are not counted in the coverage calculation. In such a case, the allelic fraction is
   * set to -1. By the same argument, the allelic fraction could rarely be > 1.
   * @param t_cov Base-pair level Coverage object, with coverage for all reads from Tumor bam(s).
   * @param n_cov Base-pair level Coverage object, with coverage for all reads from Normal bam(s).
   */
  void addAllelicFraction(STCoverage * t_cov, STCoverage * n_cov);
  
  /*! @function get the span of the breakpoints (in bp). -1 for interchrom
   * @return int distance between breakpoints
   */
  int getSpan() const { 
    if (isindel && insertion == "" )// deletion
      return (abs((int)gr1.pos1 - (int)gr2.pos1) - 1);
    if (isindel)
      return (insertion.length()); // insertion
    if (gr1.chr == gr2.chr)
      return abs((int)gr1.pos1-(int)gr2.pos1);
    else
      return -1;
  }

  /*! @function check the breakpoint against a panel of normals
   * @param Panel of normals hash
   * @return number of normal samples with this variant
   */
  int checkPon(std::unique_ptr<PON> &p);
  
  /*! @function get a unique string representation of this breakpoint.
   * Format for indel is chr_breakpos_type (eg. 0_134134_I)
   * @return string with breakpoint info
   */
  std::string getHashString() const;

  bool hasMinimal() const;

  std::string toString() const; 
 
  bool sameBreak(BreakPoint &bp) const;

  void order();

  bool isEmpty() const { return (gr1.pos1 == 0 && gr2.pos1 == 0); }

  // return whether a bp is good to move on
  bool isGoodSomatic(int mapq, size_t tsplit_cutoff, size_t nsplit_cutoff) const;

  std::string toFileString(bool noreads = false);
  
  bool hasDiscordant() const;

  // return whether a bp is good to move on
  bool isGoodGermline(int mapq, size_t allsplit) const;

  bool operator==(const BreakPoint& bp) const;

  // define how to sort these 
  bool operator < (const BreakPoint& bp) const { 
    return (gr1 < bp.gr1) || (gr1 == bp.gr1 && gr2 < bp.gr2) || 
      (gr1 == bp.gr1 && gr2 == bp.gr2 && nsplit > bp.nsplit) || // more read support should go first, do to property of std::unique
      (gr1 == bp.gr1 && gr2 == bp.gr2 && tsplit > bp.tsplit) || 
      (gr1 == bp.gr1 && gr2 == bp.gr2 && dc.ncount > bp.dc.ncount) || 
      (gr1 == bp.gr1 && gr2 == bp.gr2 && dc.tcount > bp.dc.tcount);
    //(bp.gr1.ref == refID1 && bp.pos1 > pos1) || // low pos is first
      //  (bp.refID1 == refID1 && bp.pos1 == pos1 && bp.pos2 == pos2 && nsplit1 > bp.nsplit1) || // if same, check nsplit
      // (bp.refID1 == refID1 && bp.pos1 == pos1 && bp.pos2 == pos2 && nsplit1 == bp.nsplit1 && tsplit1 > bp.tsplit1); // if also same, check tsplit
  }
  friend std::ostream& operator<<(std::ostream& out, const BreakPoint& bp) { out << bp.toString(); return out; }
  
  // print to file
  //void printToFile(ofstream &of, const BamAlignmentVector &bamreads);
  
};

}

#endif
