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

/** Walk along a BAM/SAM/CRAM/STDIN and stream in reads
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

  /** Set a part of the BAM/CRAM (ie able to be indexed) to walk.
   *
   * This will set the BAM pointer to the given region.
   * @param gp Location to point the BAM to
   * @return true if the region is found in the index
   */
  bool SetRegion(const GenomicRegion& gp);

  /** Explicitly set a reference genome to be used to decode CRAM file.
   * If no reference is specified, will automatically load from
   * file pointed to in CRAM header using the @SQ tags. 
   * @note This function is useful if the reference path pointed
   * to by the UR field of @SQ is not on your system, and you would
   * like to explicitly provide one.
   * @param ref Path to an index reference genome
   * @return Returns true if reference loaded.
   * @exception Throws an invalid_argument if reference cannot be loaded
   */
  bool SetCramReference(const std::string& ref);
  
  /** Set up multiple regions. Overwrites current regions. 
   * 
   * This will set the BAM pointer to the first element of the
   * input list.
   * @param grv Set of location to point BAM to
   * @return true if the regions are found in the index
   */
  bool SetMultipleRegions(const GRC& grv);

  /** Return if the reader has opened the file */
  bool IsOpen() const { return fp_htsfile != 0; }

  /** Create a string representation of 
   * all of the regions to walk
   */
  std::string PrintRegions() const;

  /** Print out some basic info about this reader */
  friend std::ostream& operator<<(std::ostream& out, const BamReader& b);

  /** Open a BAM file for streaming in
   */
  bool Open(const std::string& bam);

  /** Retrieve the next read from the BAM.
   *
   * If a ReadFilterCollection is defined for this BAM
   * will grab the next valid read.
   * r Read to fill with data
   * rule bool identifying if this read passed the rules
   * @return true if the next read is available
   */
  bool GetNextRecord(BamRecord &r);

  /** Return the ReadFilterCollection as a string */
  std::string displayReadFilterCollection() const;
  
  const BamHeader& Header() const { return m_hdr; }
  
  /** Reset all the counters and regions, but keep the loaded index */
  void Reset();

 protected:

  std::string m_in; ///< file name

  // point index to this region of bam
  bool __set_region(const GenomicRegion& gp);

  // open bam, true if success
  bool __open_BAM_for_reading();

  // define the regions to walk
  size_t m_region_idx = 0;

  // regions to access
  GRC m_region;

  // hts
  std::shared_ptr<htsFile> fp_htsfile;
  std::shared_ptr<hts_idx_t> idx;
  std::shared_ptr<hts_itr_t> hts_itr;

  BamHeader m_hdr;

  // hold the reference for CRAM reading
  std::string m_cram_reference;

};


}
#endif 


