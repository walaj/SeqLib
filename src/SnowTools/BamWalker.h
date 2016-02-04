#ifndef SNOWTOOLS_BAM_WALKER_H__
#define SNOWTOOLS_BAM_WALKER_H__

#include <cassert>
#include <memory>

#include "SnowTools/MiniRules.h"

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

  uint64_t keep = 0;
  uint64_t total = 0;
  
  /** Return the percent of total reads kept
   */
  int percent () const {
    int perc  = SnowTools::percentCalc<uint64_t>(keep, total); 
    return perc;
  }

  /** Return the total reads visited as a comma-formatted string
   */
  std::string totalString() const {
    return SnowTools::AddCommas<uint64_t>(total);
  }

  /** Return the kept reads as a comma-formatted string
   */
  std::string keepString() const {
    return SnowTools::AddCommas<uint64_t>(keep);
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
  void setBamWalkerRegion(const GenomicRegion& gp, std::shared_ptr<hts_idx_t> passed_idx = std::shared_ptr<hts_idx_t>(nullptr));

  /** Set up multiple regions. Overwrites current regions. 
   * 
   * This will set the BAM pointer to the first element of the
   * input list.
   * @param grv Set of location to point BAM to
   */
  void setBamWalkerRegions(const GenomicRegionVector& grv, std::shared_ptr<hts_idx_t> passed_idx = std::shared_ptr<hts_idx_t>(nullptr));

  /** Create the index file for the output bam in BAI format.
   *
   * This will make a call to HTSlib bam_index_build for the output file.
   */
  void makeIndex();

  /** Print a run-time message to stdout.
   *
   * Prints a message about all of the reads that have been visited, and informaiton
   * about the current read
   */
  void printRuntimeMessage(const ReadCount &rc_main, const BamRead &r) const;

  /** Print out some basic info about this walker, 
   * including Minz0iRules
   */
  friend std::ostream& operator<<(std::ostream& out, const BamWalker& b);

  /** Open a BAM file for streaming in
   */
  bool OpenReadBam(const std::string& bam);

  /** Open an NGS file direct from stdin */
  //bool OpenReadBam(FILE * stdin);

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
   * rule bool identifying if this read passed the rules
   * @return true if the next read is available
   */
  bool GetNextRead(BamRead &r, bool& rule);

  /** Set regions to not be read from
   */
  void addBlacklist(GRC& bl);

  /** Write an alignment to the output BAM file 
   * @param r The BamRead to save
   */
  void writeAlignment(BamRead &r);

  /** Return the MiniRulesCollection object used by this BamWalker
   */
  const MiniRulesCollection& GetMiniRulesCollection() const { return m_mr; }

  /** Send the counts for passed rules to a file */
  void MiniRulesToFile(const std::string& file) const { m_mr.countsToFile(file); }

  /** Set the BamWalker to count reads for all rules */
  void setCountAllRules() { m_mr.m_fall_through = true; }


  std::string m_in;
  std::string m_out;
  MiniRulesCollection m_mr;

  /** Set a flag to say if we should print reads to stdout
   */
  void setStdout();

  /** Set a flag to say if we should print reads to CRAM format
   */
  void setCram(const std::string& out, const std::string& ref);

  void setPrintHeader() { 
    m_print_header = true;
  }

  /** Set the output bam to remove all alignment tags */
  void setStripAllTags() { m_strip_all_tags = true; }

  /** Set a list of tags to strip */
  void setStripTags(const std::string& list);
  
  /** Set to have verbose actions */
  void setVerbose() { m_verbose = true; }

  /** Return the MiniRulesCollection as a string */
  std::string displayMiniRulesCollection() const;

  /** Return a pointer to the BAM header */
  bam_hdr_t * header() const { return br.get(); };

  /** Explicitly provide the output BAM a header */
  void SetWriteHeader(bam_hdr_t* hdr);

  /** Set the limit for total number of reads seen */
  void setReadLimit(int lim) { m_limit = lim; m_num_reads_seen = 0; }

  /** Set the limit for total number of reads kept */
  void setReadKeepLimit(int lim) { m_keep_limit = lim; m_num_reads_kept = 0; }
  
  /** Reset all the counters and regions, but keep the loaded index */
  void resetAll();

 protected:

  bool m_region_fail = false;

  // for stdout mode, print header?
  bool m_print_header = false;

  // limit to number of reads to be read in
  int m_limit = -1;
  int m_keep_limit = -1;
  
  //void __check_regions_blacklist();

  // point index to this region of bam
  bool __set_region(const GenomicRegion& gp, std::shared_ptr<hts_idx_t> passed_idx = std::shared_ptr<hts_idx_t>(nullptr));

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
  std::shared_ptr<BGZF> fp;
  //BGZF * fp = 0;
  //hts_idx_t * idx = 0; // hts_idx_load(bamfile.c_str(), HTS_FMT_BAI);
  std::shared_ptr<hts_idx_t> idx;
  std::shared_ptr<hts_itr_t> hts_itr;
  //hts_itr_t * hts_itr = 0; // sam_itr_queryi(idx, 3, 60000, 80000);
  std::shared_ptr<bam_hdr_t> br;
  std::shared_ptr<bam_hdr_t> hdr_write;
  //bam_hdr_t * br = 0;

  std::shared_ptr<htsFile> fop;
  //htsFile* fop = nullptr;

  // which tags to strip
  std::vector<std::string> m_tag_list;

  // blacklist
  GRC blacklist;

  // 
  bool m_verbose = false;

  // read counter
  int m_num_reads_seen = 0;
  int m_num_reads_kept = 0;
};


}
#endif 


