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

    if (!fop) 
      throw std::runtime_error("BamWriter::WriteHeader - Output not open for writing. Open with Open()");
    
    if (sam_hdr_write(fop.get(), hdr.get()) < 0) 
      throw std::runtime_error("BamWriter - Cannot write header. sam_hdr_write exited with < 0");
    
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
  
  // call to htslib to build bai index
  if (sam_index_build(m_out.c_str(), 0) < 0) // 0 is "min_shift", which is 0 for bai index
    throw std::runtime_error("BamWriter::makeIndex - Failed to create index");

}

  void BamWriter::Open(const std::string& f) {

    // don't reopen
    if (fop)
      return;

    m_out = f;

    // hts open the writer
    fop = std::shared_ptr<htsFile>(hts_open(m_out.c_str(), output_format.c_str()), htsFile_delete());

    if (!fop)
      throw std::runtime_error("BamWriter::Open - Cannot open output file: " + f);
  }

  BamWriter::BamWriter(int o) {

    switch(o) {
    case BAM :  output_format = "wb"; break;
    case CRAM : output_format = "wc"; break;
    case SAM :  output_format = "w"; break;
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

bool BamWriter::SetCramReference(const std::string& ref) {

  if (!fop)
    return false;

  // need to open reference for CRAM writing 
  char* fn_list = samfaipath(ref.c_str());
  if (fn_list) {
    if (hts_set_fai_filename(fop.get(), fn_list) != 0) {
      fprintf(stderr, "Failed to use reference \"%s\".\n", fn_list);
    }
  } else {
    std::cerr << "Failed to get the reference for CRAM compression" << std::endl;
  }
}

}
