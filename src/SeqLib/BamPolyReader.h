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

  };
  
/** Walk along a BAM or along BAM regions and stream in/out reads
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

  /** Set a part of the BAM to walk.
   *
   * This will set the BAM pointer to the given region.
   * @param gp Location to point the BAM to
   * @return true if the region is found in the index
   */
  bool setBamReaderRegion(const GenomicRegion& gp);

  /** Set up multiple regions. Overwrites current regions. 
   * 
   * This will set the BAM pointer to the first element of the
   * input list.
   * @param grv Set of location to point BAM to
   * @return true if the regions are found in the index
   */
  bool setBamReaderRegions(const GenomicRegionVector& grv);

  /** Create a string representation of 
   * all of the regions to walk
   */
  std::string printRegions() const;

  /** Print out some basic info about this walker, 
   * including Minz0iRules
   */
  friend std::ostream& operator<<(std::ostream& out, const BamPolyReader& b);

  /** Open a BAM file for streaming in
   */
  bool OpenReadBam(const std::string& bam);

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

  /** Reset all the counters and regions, but keep the loaded index */
  void resetAll();

  /** Return a header to the first file */
  BamHeader Header() const { if (m_bams.size()) return m_bams[0].m_hdr; return BamHeader(); }

 protected:

  ReadFilterCollection m_mr; ///< filter collection

  // point index to this region of bam
  bool __set_region(const GenomicRegion& gp);

  // open bam, true if success
  bool __open_BAM_for_reading();

  // define the regions to walk
  size_t m_region_idx = 0;

  GenomicRegionVector m_region;

  std::vector<_Bam> m_bams;

};


}
#endif 


