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

// from samtools
inline char *samfaipath(const char *fn_ref)
{
  char *fn_list = 0;
  if (fn_ref == 0) return 0;
  fn_list = (char*)calloc(strlen(fn_ref) + 5, 1);
  strcat(strcpy(fn_list, fn_ref), ".fai");
  if (access(fn_list, R_OK) == -1) { // fn_list is unreadable
    std::cerr << "ERROR: Cannot read the index file for CRAM read/write" << std::endl;
  }
  return fn_list;
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

  std::string totalString() const {
    return SnowTools::AddCommas<int>(total);
  }

  std::string keepString() const {
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
  //BamWalker(const std::string& in, const std::string& out);

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
    __close_read_bam();
    __close_write_bam();
  }

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

  /** Create the index file for the output bam in BAI format.
   *
   * This will make a call to HTSlib bam_index_build for the output file.
   */
  void MakeIndex();

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

  /** Set a flag to say if we should print reads to stdout
   */
  void setStdout() { fop = sam_open("-", "w"); }

  /** Set a flag to say if we should print reads to CRAM format
   */
  void setCram(const std::string& out, const std::string& ref) { 
    m_out = out;
    fop = sam_open(m_out.c_str(), "wc"); 
    if (!fop) {
      std::cerr << "!!!\n!!!\n!!!\nCannot open CRAM file for writing. Will try BAM next. File: " <<  m_out << std::endl;
      return;
    }

    // need to open reference for CRAM writing 
    char* fn_list = samfaipath(ref.c_str());
    if (fn_list) {
      if (hts_set_fai_filename(fop, fn_list) != 0) {
	fprintf(stderr, "Failed to use reference \"%s\".\n", fn_list);
      }
    } else {
      std::cerr << "Failed to get the reference for CRAM compression" << std::endl;
    }
    m_print_header = true; 
  }

  void setPrintHeader() { 
    m_print_header = true;
  }

  /** Set the output bam to remove all alignment tags */
  void setStripAllTags() { m_strip_all_tags = true; }

  /** Set a list of tags to strip */
  void setStripTags(const std::string& list);
  
  /** Set to have verbose actions */
  void setVerbose() { m_verbose = true; }

 protected:

  // for stdout mode, print header?
  bool m_print_header = false;

  // point index to this region of bam
  bool __set_region(const GenomicRegion& gp);

  struct timespec start;

  //open bam, true if success
  bool __open_BAM_for_reading();

  // close this bam file and delete all pointers
  void __close_read_bam();
  void __close_write_bam();

  // open m_out, true if success
  bool __open_BAM_for_writing();
  
  // define the regions to walk
  size_t m_region_idx = 0;
  GenomicRegionVector m_region;

  bool m_strip_all_tags = false;

  // hts
  BGZF * fp = 0;
  hts_idx_t * idx = 0; // hts_idx_load(bamfile.c_str(), HTS_FMT_BAI);
  hts_itr_t * hts_itr = 0; // sam_itr_queryi(idx, 3, 60000, 80000);
  bam_hdr_t * br = 0;

  htsFile* fop = 0;

  // which tags to strip
  std::vector<std::string> m_tag_list;

  // 
  bool m_verbose = false;
};


}
#endif 


