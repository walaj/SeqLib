#ifndef SEQLIB_BAM_READER_H__
#define SEQLIB_BAM_READER_H__

#include <cassert>
#include <memory>

#include "SeqLib/ReadFilter.h"
#include "SeqLib/BamWalker.h"

// notes
// htsFile->format.format
//    // sam = 3
//    // bam = 4
//    // cram = 6

namespace SeqLib {

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
   * @return true if the region is found in the index
   */
  bool setBamReaderRegion(const GenomicRegion& gp, std::shared_ptr<hts_idx_t> passed_idx = std::shared_ptr<hts_idx_t>(nullptr));

  /** Explicitly set a reference genome to be used to decode CRAM file.
   * If no reference is specified, will automatically load from
   * file pointed to in CRAM header using the @SQ tags. 
   * @note This function is useful if the reference path pointed
   * to by the UR field of @SQ is not on your system, and you would
   * like to explicitly provide one.
   * @param ref Path to an index reference genome
   * @return Returns true if reference loaded.
   */
  bool SetCramReference(const std::string& ref);
  
  /** Set up multiple regions. Overwrites current regions. 
   * 
   * This will set the BAM pointer to the first element of the
   * input list.
   * @param grv Set of location to point BAM to
   * @return true if the regions are found in the index
   */
  bool setBamReaderRegions(const GenomicRegionVector& grv, std::shared_ptr<hts_idx_t> passed_idx = std::shared_ptr<hts_idx_t>(nullptr));

  /** Create a string representation of 
   * all of the regions to walk
   */
  std::string printRegions() const;

  /** Print a run-time message to stdout.
   *
   * Prints a message about all of the reads that have been visited, and informaiton
   * about the current read
   */
  void printRuntimeMessage(const ReadCount &rc_main, const BamRecord &r) const;

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
   * parse the provides rules and add it as a rule set to this BamRecorder.
   * @param rules A JSON string of filters, or a file pointing to a filter script
   */
  void SetReadFilterCollection(const std::string& rules);
  
  /** Explicitly provide a ReadFilterCollection to this BamReader
   */
  void SetReadFilterCollection(const ReadFilterCollection& mr); 

  /** Retrieve the next read from the BAM.
   *
   * If a ReadFilterCollection is defined for this BAM
   * will grab the next valid read.
   * r Read to fill with data
   * rule bool identifying if this read passed the rules
   * @return true if the next read is available
   */
  bool GetNextRead(BamRecord &r, bool& rule);

  /** Return the ReadFilterCollection object used by this BamReader
   */
  const ReadFilterCollection& GetReadFilterCollection() const { return m_mr; }

  /** Set the BamReader to count reads for all rules */
  void setCountAllRules() { m_mr.CheckAllFilters(); }

  /** Set to have verbose actions */
  void setVerbose() { m_verbose = true; }

  /** Return the ReadFilterCollection as a string */
  std::string displayReadFilterCollection() const;
  
  const BamHeader& Header() const { return m_hdr; }
  
  /** Reset all the counters and regions, but keep the loaded index */
  void resetAll();

 protected:

  std::string m_in; ///< file name

  ReadFilterCollection m_mr; ///< filter collection

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
  std::shared_ptr<htsFile> fp_htsfile;
  std::shared_ptr<hts_idx_t> idx;
  std::shared_ptr<hts_itr_t> hts_itr;

  BamHeader m_hdr;

  bool m_verbose = false;

  std::string m_cram_reference;

};


}
#endif 


