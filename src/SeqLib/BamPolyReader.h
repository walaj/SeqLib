#ifndef SEQLIB_BAM_POLYREADER_H__
#define SEQLIB_BAM_POLYREADER_H__

#include <cassert>
#include <memory>

#include "SeqLib/ReadFilter.h"
#include "SeqLib/BamWalker.h"

namespace SeqLib {
 
  // store file accessors for single BAM
  struct _Bam {
    
  _Bam(const std::string& m) : m_in(m) {}
    ~_Bam() {}
    
    std::shared_ptr<htsFile> fp;
    std::shared_ptr<hts_idx_t> idx;
    std::shared_ptr<hts_itr_t> hts_itr;
    std::string m_in;
    BamHeader m_hdr;

    BamRecord next_read;

    bool empty = true;
    
    bool open_BAM_for_reading();

    std::string id;

    // hold the reference for CRAM reading
    std::string m_cram_reference;

  };
  
/** Walk along a collection of BAM/SAM/CRAM stream in reads
 */
class BamPolyReader {

 public:

  /** Construct an empty BamPolyReader */
  BamPolyReader();

  /** Destroy a BamPolyReader and close all connections to the BAMs 
   * 
   * Calling the destructor will take care of all of the C-style dealloc
   * calls required within HTSlib to close a BAM or SAM file. 
   */
  ~BamPolyReader() { }

  /** Explicitly set a reference genome to be used to decode CRAM file.
   * If no reference is specified, will automatically load from
   * file pointed to in CRAM header using the @SQ tags. 
   * @note This function is useful if the reference path pointed
   * to by the UR field of @SQ is not on your system, and you would
   * like to explicitly provide one.
   * @param ref Path to an index reference genome
   */
  void SetCramReference(const std::string& ref);

  /** Set a part of the BAM to walk.
   *
   * This will set the BAM pointer to the given region.
   * @param gp Location to point the BAM to
   * @return true if the region is found in the index
   */
  bool SetRegion(const GenomicRegion& gp);

  /** Set up multiple regions. Overwrites current regions. 
   * 
   * This will set the BAM pointer to the first element of the
   * input list.
   * @param grc Set of location to point BAM to
   * @return true if the regions are found in the index
   */
  bool SetMultipleRegions(const GRC& grc);

  /** Create a string representation of 
   * all of the regions to walk
   */
  std::string PrintRegions() const;

  /** Print out some basic info about this reader */
  friend std::ostream& operator<<(std::ostream& out, const BamPolyReader& b);

  /** Open a BAM/SAM/CRAM/STDIN file for streaming in 
   * @param bam Path to a SAM/CRAM/BAM file, or "-" for stdin
   * @return True if open was successful
   */
  bool Open(const std::string& bam);

  /** Open a set of BAM/SAM/CRAM/STDIN files for streaming in 
   * @param bams Path to a vector fo SAM/CRAM/BAM files, or "-" for stdin
   * @return True if open was successful
   */
  bool Open(const std::vector<std::string>& bams);

  /** Retrieve the next read from the available input streams.
   * @note Will chose the read with the lowest left-alignment position
   * from the available streams.
   * @param r Read to fill with data
   * @return true if the next read is available
   */
  bool GetNextRecord(BamRecord &r);

  /** Reset all the counters and regions, but keep the loaded index */
  void Reset();

  /** Return a header to the first file */
  BamHeader Header() const { if (m_bams.size()) return m_bams[0].m_hdr; return BamHeader(); }

 protected:

  // point index to this region of bam
  bool __set_region(const GenomicRegion& gp);

  // which region are we on
  size_t m_region_idx = 0;

  // regions to walk
  GRC m_region;

  // store the file pointers etc to BAM files
  std::vector<_Bam> m_bams;

  // hold the reference for CRAM reading
  std::string m_cram_reference;

};


}
#endif 


