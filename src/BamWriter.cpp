#include "SeqLib/BamWalker.h"
#include "SeqLib/BamWriter.h"

//#define DEBUG_WALKER 1

namespace SeqLib {

  void BamWriter::SetHeader(const SeqLib::BamHeader& h) {
    hdr = h;
  }

  void BamWriter::WriteHeader() const {
    
    if (hdr.isEmpty())
      throw std::runtime_error("BamWriter::WriteHeader - No header supplied. Provide with SetWriteHeader");

    // hts open the writer
    //if (!fop) { // default is bam. if already set by flag, don't reopen
    //  assert(m_out.length());
    //  fop = std::shared_ptr<htsFile>(sam_open(m_out.c_str(), "wb"), sam_write_delete());
    // }
    
    if (!fop) 
      throw std::runtime_error("BamWriter::WriteHeader - Output not open for writing. Open with Open()");
    
    if (sam_hdr_write(fop.get(), hdr.get()) < 0) 
      throw std::runtime_error("BamWriter - Cannot write header. sam_hdr_write exited with < 0");
    
  }
  
  void BamWriter::SetWriteHeader(bam_hdr_t* hdr) {
    hdr_write = std::shared_ptr<bam_hdr_t>(bam_hdr_dup(hdr), bam_hdr_delete()); 
  }

  void BamWriter::CloseBam() {

    if (!fop)
      throw std::runtime_error("BamWriter::CloseBam() - Trying to close BAM that is already closed or never opened");

    fop = nullptr; // this clears shared_ptr, calls sam_close
  }

void BamWriter::makeIndex() const {
  
  // throw an error if BAM is not already closed
  if (fop)
    throw std::runtime_error("BamWriter::makeIndex - Trying to index open BAM. Close first with CloseBam()");

  if (m_out.empty())
    throw std::runtime_error("Trying to make index, but no BAM specified");    
  
  std::cerr << " BUILDING " << m_out.c_str() << std::endl;
  std::cerr << " INDEX OUT " << sam_index_build(m_out.c_str(), 0) << std::endl;
  // call to htslib to build bai index
  if (sam_index_build(m_out.c_str(), 0) < 0) // 0 is "min_shift", which is 0 for bai index
    throw std::runtime_error("BamWriter::makeIndex - Failed to create index");
  std::cerr << " FAILED TRINY TO BUILD INDEX FOR " << m_out << std::endl;  

}

  void BamWriter::Open(const std::string& f) {

    // don't reopen
    if (fop)
      return;

    m_out = f;

    // hts open the writer
    fop = std::shared_ptr<htsFile>(sam_open(m_out.c_str(), output_format.c_str()), sam_write_delete());

    if (!fop)
      throw std::runtime_error("BamWriter::Open - Cannot open output file: " + f);
  }

  BamWriter::BamWriter(int o) {

    switch(o) {
    case BAM : output_format = "wb";
    case CRAM : output_format = "wc";
    case SAM : output_format = "w";
    }

  }
  

void BamWriter::writeAlignment(BamRecord &r)
{
  if (!fop) {
    throw std::runtime_error("BamWriter::writeAlignment - Cannot write BamRecord. Did you forget to open the Bam for writing (OpenWriteBam)?");
  } else {
    if (sam_write1(fop.get(), hdr.get(), r.raw()) < 0)
      throw std::runtime_error("BamWriter::writeAlignment - Cannot write BamRecord. sam_write1 exited with < 0");      
  }
}

std::ostream& operator<<(std::ostream& out, const BamWriter& b)
{
  
  out << "Write "; 
  if (!b.fop)
    out << "BAM" << std::endl;
  else if (b.fop->format.format == 4)
    out << "BAM" << std::endl;
  else if (b.fop->format.format == 6)
    out << "CRAM" << std::endl;
  else if (b.fop->format.format == text_format)
    out << "SAM" << std::endl;

  out << ":" << b.m_out;
  return out;
}

void BamWriter::setCram(const std::string& ref) {

  fop = std::shared_ptr<htsFile>(sam_open(m_out.c_str(), "wc"), sam_write_delete()); 
  if (!fop) {
    std::cerr << "!!!\n!!!\n!!!\nCannot open CRAM file for writing. Will try BAM next. File: " <<  m_out << std::endl;
    return;
  }
  
  // need to open reference for CRAM writing 
  char* fn_list = samfaipath(ref.c_str());
  if (fn_list) {
    if (hts_set_fai_filename(fop.get(), fn_list) != 0) {
      fprintf(stderr, "Failed to use reference \"%s\".\n", fn_list);
    }
  } else {
    std::cerr << "Failed to get the reference for CRAM compression" << std::endl;
  }
  m_print_header = true; 
}

void BamWriter::setStdout() {
  fop = std::shared_ptr<htsFile>(sam_open("-", "w"), sam_write_delete()); 
}

}
