#ifndef SEQLIB_BAM_WRITER_H__
#define SEQLIB_BAM_WRITER_H__

#include <cassert>
#include <memory>
#include "SeqLib/BamRecord.h"

namespace SeqLib {

  const int BAM = 4;
  const int SAM = 3;
  const int CRAM = 6;

/** Walk along a BAM or along BAM regions and stream in/out reads
 */
class BamWriter  {

 public:

  /** Construct an empty BamWriter to write BAM */
  BamWriter() {}

  /** Construct an empty BamWriter and specify output format 
   * @param o One of SeqLib::BAM, SeqLib::CRAM, SeqLib::SAM, SeqLib::STDOUT
   */
  BamWriter(int o);

  /** Destroy a BamWriter and close all connections to the BAM 
   * 
   * Calling the destructor will take care of all of the C-style dealloc
   * calls required within HTSlib to close a BAM or SAM file. 
   */
  ~BamWriter() {}

  /** Write the BAM header */
  void WriteHeader() const;

  /** Provide a header to this writer 
   * @param h Header for this writer. Copies contents
   */
  void SetHeader(const SeqLib::BamHeader& h);

  /** Close a file explitily. This is required before indexing with makeIndex.
   * @note If not called, BAM will close properly on object destruction
   * @exception Throws a runtime_error if BAM already closed or was never opened
   */
  void Close();

  /** Create the index file for the output bam in BAI format.
   *
   * This will make a call to HTSlib bam_index_build for the output file. 
   * @exception Throws a runtime_error if sam_index_build2 exits with < 0 status
   */
  void BuildIndex() const;

  /** Print out some basic info about this writer */
  friend std::ostream& operator<<(std::ostream& out, const BamWriter& b);

  /** Open a BAM file for streaming out.
   * @param f Path to the output BAM file
   * @exception Throws a runtime_error if cannot write BAM
   * @note Calling this function will immediately write the BAM with its header
   */
  void Open(const std::string& f);

  /** Write an alignment to the output BAM file 
   * @param r The BamRecord to save
   * @exception Throws a runtime_error if cannot write alignment
   */
  void WriteRecord(BamRecord &r);

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

  /** Return the BAM header */
  BamHeader Header() const { return hdr; };

 protected:

  // path to output file
  std::string m_out; 

  // open m_out, true if success
  void __open_BAM_for_writing();
  
  // output format
  std::string output_format = "wb";
  
  // hts
  std::shared_ptr<htsFile> fop;

  // header
  SeqLib::BamHeader hdr;
  
};


}
#endif 


