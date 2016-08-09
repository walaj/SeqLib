#ifndef SNOWTOOLS_BAM_WRITER_H__
#define SNOWTOOLS_BAM_WRITER_H__

#include <cassert>
#include <memory>

namespace SnowTools {

/** Walk along a BAM or along BAM regions and stream in/out reads
 */
class BamWriter  {

 public:

  /** Construct a new BamWriter for writing a BAM/SAM/CRAM
   * @param in Name of BAM/SAM/CRAM file to write
   */
  BamWriter(const std::string& in);

  /** Construct an empty BamWriter */
  BamWriter() {}

  /** Destroy a BamWriter and close all connections to the BAM 
   * 
   * Calling the destructor will take care of all of the C-style dealloc
   * calls required within HTSlib to close a BAM or SAM file. 
   */
  ~BamWriter() {}

  /** Create the index file for the output bam in BAI format.
   *
   * This will make a call to HTSlib bam_index_build for the output file.
   */
  void makeIndex();

  /** Print out some basic info about this walker, 
   * including Minz0iRules
   */
  friend std::ostream& operator<<(std::ostream& out, const BamWriter& b);

  /** Open a BAM file for streaming out
   */
  bool OpenWriteBam(const std::string& bam);

  /** Write an alignment to the output BAM file 
   * @param r The BamRead to save
   */
  void writeAlignment(BamRead &r);

  std::string m_in;
  std::string m_out;

  /** Set a flag to say if we should print reads to stdout */
  void setStdout();

  /** Set a flag to say if we should print reads to CRAM format
   * @param out Output CRAM file to write to
   * @param ref File with the reference genome used for compression
   */
  void setCram(const std::string& ref);

  /** If set to true, will print header in output 
   * @param val Set whether to print the hader
   */
  void setPrintHeader(bool val) { 
    m_print_header = val;
  }
  
  /** Return a pointer to the BAM header */
  bam_hdr_t * header() const { return br.get(); };

  /** Explicitly provide the output BAM a header. 
   * 
   * This will create a copy of the bam_hdr_t for this BamWriter.
   */
  void SetWriteHeader(bam_hdr_t* hdr);

 protected:

  // for stdout mode, print header?
  bool m_print_header = false;

  // open m_out, true if success
  bool __open_BAM_for_writing();
  
  // hts
  std::shared_ptr<BGZF> fp;
  std::shared_ptr<hts_idx_t> idx;
  std::shared_ptr<hts_itr_t> hts_itr;
  std::shared_ptr<bam_hdr_t> br;
  std::shared_ptr<bam_hdr_t> hdr_write;

  std::shared_ptr<htsFile> fop;

  // which tags to strip
  std::vector<std::string> m_tag_list;

};


}
#endif 


