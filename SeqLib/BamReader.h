#ifndef SEQLIB_BAM_POLYREADER_H__
#define SEQLIB_BAM_POLYREADER_H__

#include <cassert>
#include <memory>
#include <unordered_map>

#include "SeqLib/ReadFilter.h"
#include "SeqLib/BamWalker.h"

// forward declare this from hts.c
extern "C" {
int hts_useek(htsFile *file, long uoffset, int where);
}

class BamReader;

namespace SeqLib {
 
  // store file accessors for single BAM
  class _Bam {

    friend class BamReader;

  public:
    _Bam(const std::string& m) : m_in(m) {}
    _Bam() {}
    ~_Bam() {}

    std::string GetFileName() const { return m_in; }

  private:

    // do the read loading
    bool __load_read(BamRecord& r);

    void reset() {
      empty = true;
      mark_for_closure = false;
      m_region_idx = 0;
      
      // close and reopen, only way I know how to do right now
      //fp = nullptr; // calls destructor on htsFile
      //hts_itr = nullptr;
      
      // re-read it in
      //fp = std::shared_ptr<htsFile>(hts_open(m_in.c_str(), "r"), htsFile_delete()); 
      
      // other attempts at a true rewind
      //hts_useek(fp.get(), 0, 0);
      //ks_rewind((kstream_t*)fp->fp.voidp);
      //hts_itr = nullptr; // reset to beginning

    }

    // close this bam
    bool close() {
      if (!fp)
	return false;
      fp = nullptr; // calls destructor actually
      idx = nullptr;
      hts_itr = nullptr;

      empty = true;
      mark_for_closure = false;
      m_region_idx = 0;

      return true;
    }

    GRC* m_region; // local copy of region

    std::shared_ptr<htsFile> fp;     // BAM file pointer
    std::shared_ptr<hts_idx_t> idx;  // bam index
    std::shared_ptr<hts_itr_t> hts_itr; // iterator to index location
    std::string m_in;                   // file name
    BamHeader m_hdr;                    // the BAM header

    // the next read "slotted" for this BAM
    BamRecord next_read;

    // the next read "slot" is empty
    bool empty = true;
    
    // if set to true, then won't even attempt to lookup read
    bool mark_for_closure = false;
    
    // open the file pointer
    bool open_BAM_for_reading();

    // hold the reference for CRAM reading
    std::string m_cram_reference;

    // which region are we on
    size_t m_region_idx = 0;

    // point index to this region of bam
    bool __set_region(const GenomicRegion& gp);

    

  };
  
/** Stream in reads from multiple BAM/SAM/CRAM or stdin */
class BamReader {

 public:

  /** Construct an empty BamReader */
  BamReader();

  /** Destroy a BamReader and close all connections to the BAMs 
   * 
   * Calling the destructor will take care of all of the C-style dealloc
   * calls required within HTSlib to close a BAM or SAM file. 
   */
  ~BamReader() { }

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
   * @note This clears all other regions and resets the index
   * pointer to this location
   * @param gp Location to point the BAM to
   * @return true if the region is found in the index
   */
  bool SetRegion(const GenomicRegion& gp);

  /** Set up multiple regions. Overwrites current regions. 
   * 
   * This will set the BAM pointer to the first element of the
   * input list.
   * @note This clears all other regions and resets the index
   * pointer to the first element of grc
   * @param grc Set of location to point BAM to
   * @return true if the regions are found in the index
   */
  bool SetMultipleRegions(const GRC& grc);

  /** Return if the reader has opened the first file */
  bool IsOpen() const { if (m_bams.size()) return m_bams.begin()->second.fp != 0; return false; }

  /** Return if the reader has opened the file
   * @param f Name of file to check
   */
  bool IsOpen(const std::string& f) const { 
    std::unordered_map<std::string, _Bam>::const_iterator ff = m_bams.find(f);
    if (ff == m_bams.end())
      return false;
    return ff->second.fp != 0; 
  }

  /** Close all of the BAMs */
  bool Close();

  /** Close a particular BAM/CRAM/SAM
   * @param f Particular file to close
   * @return True if BAM is found and is closable (eg no already closed)
   */
  bool Close(const std::string& f);

  /** Reset the given BAM/SAM/CRAM to the begining, but keep the loaded indicies and file-pointers 
   * @param f Name of file to reset
   * @return Returns false if this BAM is not found in object
   * @note Unlike Reset(), this version will NOT reset the regions, since other BAMs may still be
   * using them.
   */
  bool Reset(const std::string& f);

  /** Return a vector of all of the BAM/SAM/CRAMs in this reader */
  std::vector<std::string> ListFiles() const {
    std::vector<std::string> out;
    for (auto& i : m_bams)
      out.push_back(i.first);
    return out;
  }

  /** Create a string representation of 
   * all of the regions to walk
   */
  std::string PrintRegions() const;

  /** Print out some basic info about this reader */
  friend std::ostream& operator<<(std::ostream& out, const BamReader& b);

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

  /** Reset all the regions, but keep the loaded indicies and file-pointers */
  void Reset();

  /** Return a header to the first file */
  BamHeader Header() const { if (m_bams.size()) return m_bams.begin()->second.m_hdr; return BamHeader(); }

 protected:

  // regions to walk
  GRC m_region;

  // store the file pointers etc to BAM files
  std::unordered_map<std::string, _Bam> m_bams;

  // hold the reference for CRAM reading
  std::string m_cram_reference;

};


}
#endif 


