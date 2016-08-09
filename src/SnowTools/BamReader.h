#ifndef SNOWTOOLS_BAM_READER_H__
#define SNOWTOOLS_BAM_READER_H__

#include <cassert>
#include <memory>

#include "SnowTools/ReadFilter.h"
#include "SnowTools/BamWalker.h"

namespace SnowTools {

/** Walk along a BAM or along BAM regions and stream in/out reads
 */
class BamReader {

 public:

  /** Construct a new BamReader for reading a BAM/SAM/CRAM */
  BamReader(const std::string& in);

  /** Construct an empty BamReader */
  BamReader();

  /** Destroy a BamReader and close all connections to the BAM 
   * 
   * Calling the destructor will take care of all of the C-style dealloc
   * calls required within HTSlib to close a BAM or SAM file. 
   */
  ~BamReader() { }

  /** Set a part of the BAM to walk.
   *
   * This will set the BAM pointer to the given region.
   * @param gp Location to point the BAM to
   */
  void setBamReaderRegion(const GenomicRegion& gp, std::shared_ptr<hts_idx_t> passed_idx = std::shared_ptr<hts_idx_t>(nullptr));

  /** Set up multiple regions. Overwrites current regions. 
   * 
   * This will set the BAM pointer to the first element of the
   * input list.
   * @param grv Set of location to point BAM to
   */
  void setBamReaderRegions(const GenomicRegionVector& grv, std::shared_ptr<hts_idx_t> passed_idx = std::shared_ptr<hts_idx_t>(nullptr));

  /** Create a string representation of 
   * all of the regions to walk
   */
  std::string printRegions() const;

  /** Print a run-time message to stdout.
   *
   * Prints a message about all of the reads that have been visited, and informaiton
   * about the current read
   */
  void printRuntimeMessage(const ReadCount &rc_main, const BamRead &r) const;

  /** Print out some basic info about this walker, 
   * including Minz0iRules
   */
  friend std::ostream& operator<<(std::ostream& out, const BamReader& b);

  /** Open a BAM file for streaming in
   */
  bool OpenReadBam(const std::string& bam);

  /** Pass a ReadFilter script to the BAM.
   * 
   * This will call the constructor of ReadFilterCollection, and 
   * parse the provides rules and add it as a rule set to this BamReader.
   * @param rules A string of rules, or a file pointing to a rules script
   */
  void SetReadFilterCollection(const std::string& rules);
  
  /** Explicitly provide a ReadFilterCollection to this BamReader
   */
  void SetReadFilterCollection(const ReadFilterCollection& mr) { m_mr = mr; }

  /** Retrieve the next read from the BAM.
   *
   * If a ReadFilterCollection is defined for this BAM
   * will grab the next valid read.
   * r Read to fill with data
   * rule bool identifying if this read passed the rules
   * @return true if the next read is available
   */
  bool GetNextRead(BamRead &r, bool& rule);

  /** Return the ReadFilterCollection object used by this BamReader
   */
  const ReadFilterCollection& GetReadFilterCollection() const { return m_mr; }

  /** Send the counts for passed rules to a file 
   * @param Path to file to write the read statistics to
   */
  void ReadFilterToFile(const std::string& file) const { m_mr.countsToFile(file); }

  /** Set the BamReader to count reads for all rules */
  void setCountAllRules() { m_mr.m_fall_through = true; }

  std::string m_in;

  ReadFilterCollection m_mr;

  /** Set to have verbose actions */
  void setVerbose() { m_verbose = true; }

  /** Return the ReadFilterCollection as a string */
  std::string displayReadFilterCollection() const;

  /** Return a pointer to the BAM header */
  bam_hdr_t * header() const { return br.get(); };

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
  
  // point index to this region of bam
  bool __set_region(const GenomicRegion& gp, std::shared_ptr<hts_idx_t> passed_idx = std::shared_ptr<hts_idx_t>(nullptr));

  // open bam, true if success
  bool __open_BAM_for_reading();

  // define the regions to walk
  size_t m_region_idx = 0;
  GenomicRegionVector m_region;

  struct timespec start;

  // hts
  std::shared_ptr<BGZF> fp;
  std::shared_ptr<hts_idx_t> idx;
  std::shared_ptr<hts_itr_t> hts_itr;
  std::shared_ptr<bam_hdr_t> br;

  bool m_verbose = false;

  // read counter
  int m_num_reads_seen = 0;
  int m_num_reads_kept = 0;
};


}
#endif 


