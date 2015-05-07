#ifndef SNOWTOOLS_BAM_WALKER_H__
#define SNOWTOOLS_BAM_WALKER_H__

#include <cassert>

#include "SnowTools/MiniRules.h"
#include "SnowTools/GenomicRegion.h"
#include "SnowTools/GenomicRegionCollection.h"

#include "SnowTools/HTSTools.h"
#include "SnowTools/SnowUtils.h"

// Phred score transformations
inline int char2phred(char b) {
  uint8_t v = b;
  assert(v >= 33);
  return v - 33;
}

namespace SnowTools {

/////////////// 
// Hold read counts
//////////////
struct ReadCount {

  int keep = 0;
  int total = 0;
  
  int percent () const {
    int perc  = SnowTools::percentCalc<int>(keep, total); 
    return perc;
  }

  string totalString() const {
    return SnowTools::AddCommas<int>(total);
  }

  string keepString() const {
    return SnowTools::AddCommas<int>(keep);
  }

};

/** Walk along a BAM or along BAM regions and stream in/out reads
 */
class BamWalker {

 public:

  /** Construct a new BamWalker for streaming data in and streaming
   * out to a new BAM file.
   *
   * This constructor will open the BAM file and read the header. 
   * It will also open the out BAM file.
   */
  
  BamWalker(const std::string& in, const std::string& out);

  /** Construct a new BamWalker for reading a BAM
   */
  BamWalker(const std::string& in);

  /** Construct an empty BamWalker */
  BamWalker() {}

  /** Destroy a BamWalker and close all connections to the BAM 
   * 
   * Calling the destructor will take care of all of the C-style dealloc
   * calls required within HTSlib to close a BAM or SAM file. 
   */
  ~BamWalker() {

    if (fp)
      bgzf_close(fp);
    if (br)
      bam_hdr_destroy(br);
    if (hts_itr)
      hts_itr_destroy(hts_itr);
    if (idx)
      hts_idx_destroy(idx);
    if (fop)
      sam_close(fop);

  }

  struct timespec start;

  void printRuleCounts(unordered_map<string, size_t> &rm) const;
  
  /** Set a part of the BAM to walk.
   *
   * This will set the BAM pointer to the given region.
   * @param gp Location to point the BAM to
   */
  void setBamWalkerRegion(const GenomicRegion& gp);

  /** Set up multiple regions. Overwrites current regions. 
   * 
   * This will set the BAM pointer to the first element of the
   * input list.
   * @param grv Set of location to point BAM to
   */
  void setBamWalkerRegions(const GenomicRegionVector& grv);

  // create the index file for the output bam
  //void MakeIndex();

  /** Print a run-time message to stdout.
   *
   * Prints a message about all of the reads that have been visited, and informaiton
   * about the current read
   */
  void printRuntimeMessage(const ReadCount &rc_main, const Read &r) const;

  /** Print out some basic info about this walker, 
   * including Minz0iRules
   */
  friend std::ostream& operator<<(std::ostream& out, const BamWalker& b);

  /** Open a BAM file for streaming in
   */
  bool OpenReadBam(const std::string& bam);

  /** Open a BAM file for streaming out
   */
  bool OpenWriteBam(const std::string& bam);

  /** Pass a MiniRules script to the BAM.
   * 
   * This will call the constructor of MiniRulesCollection, and 
   * parse the provides rules and add it as a rule set to this BamWalker.
   * @param rules A string of rules, or a file pointing to a rules script
   */
  void SetMiniRulesCollection(const std::string& rules);
  
  /** Explicitly provide a MiniRulesCollection to this BamWalker
   */
  void SetMiniRulesCollection(const MiniRulesCollection& mr) { m_mr = mr; }

  /** Retrieve the next read from the BAM.
   *
   * If a MiniRulesCollection is defined for this BAM
   * will grab the next valid read.
   * r Read to fill with data
   * ref String identifying which rule this read passed. Empty if not passed
   * @return true if the next read is available
   */
  bool GetNextRead(Read &r, std::string& ref);

  /** Write an alignment to the output BAM file 
   * @param r The BamRead to save
   */
  void WriteAlignment(Read &r);

  /** Return the MiniRulesCollection object used by this BamWalker
   */
  const MiniRulesCollection& GetMiniRulesCollection() const { return m_mr; }

  std::string m_in;
  std::string m_out;
  MiniRulesCollection m_mr;

 protected:

  // point index to this region of bam
  bool __set_region(const GenomicRegion& gp);

  //open bam, true if success
  bool __open_BAM_for_reading();

  // open m_out, true if success
  bool __open_BAM_for_writing();
  
  // define the regions to walk
  size_t m_region_idx = 0;
  GenomicRegionVector m_region;

  int m_verbose = 1;

  size_t m_reads_seen = 0;
  
  size_t m_reads_seen_valid = 0;

  // hts
  BGZF * fp = 0;
  hts_idx_t * idx = 0; // hts_idx_load(bamfile.c_str(), HTS_FMT_BAI);
  hts_itr_t * hts_itr = 0; // sam_itr_queryi(idx, 3, 60000, 80000);
  bam_hdr_t * br = 0;

  samFile* fop = 0;
  //fp = sam_open(fn, mode);

};


}
#endif 


