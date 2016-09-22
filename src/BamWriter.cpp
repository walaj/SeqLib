#include "SeqLib/BamWalker.h"
#include "SeqLib/BamWriter.h"

#include <stdexcept>

//#define DEBUG_WALKER 1

namespace SeqLib {

  void BamWriter::SetHeader(const SeqLib::BamHeader& h) {
    hdr = h;
  }

  bool BamWriter::WriteHeader() const {
    
    if (hdr.isEmpty()) {
      //throw std::runtime_error("BamWriter::WriteHeader - No header supplied. Provide with SetWriteHeader");
      std::cerr << "BamWriter::WriteHeader - No header supplied. Provide with SetWriteHeader" << std::endl;
      return false;
    }

    if (!fop) {
      std::cerr << "BamWriter::WriteHeader - Output not open for writing. Open with Open()" << std::endl;
      //throw std::runtime_error("BamWriter::WriteHeader - Output not open for writing. Open with Open()");
      return false;
    }
    
    if (sam_hdr_write(fop.get(), hdr.get()) < 0) {
      //throw std::runtime_error("Cannot write header. sam_hdr_write exited with < 0");
      std::cerr << "Cannot write header. sam_hdr_write exited with < 0" << std::endl;
      return false;
    }

    return true;
    
  }
  
  bool BamWriter::Close() {

    if (!fop)
      return false;

    fop.reset(); //tr1
    //fop = NULL; // this clears shared_ptr, calls sam_close (c++11)

    return true;
  }

bool BamWriter::BuildIndex() const {
  
  // throw an error if BAM is not already closed
  if (fop) {
    //throw std::runtime_error("Trying to index open BAM. Close first with Close()");
    std::cerr << "Trying to index open BAM. Close first with Close()" << std::endl;
    return false;
  }

  if (m_out.empty()) {
    std::cerr << "Trying to make index, but no BAM specified" << std::endl;
    //throw std::runtime_error("Trying to make index, but no BAM specified");    
    return false;
  }
  
  // call to htslib to build bai index
  if (sam_index_build(m_out.c_str(), 0) < 0) { // 0 is "min_shift", which is 0 for bai index
    std::cerr << "Failed to create index";
    return false;
  }

  return true;

}

  bool BamWriter::Open(const std::string& f) {

    // don't reopen
    if (fop)
      return false;

    m_out = f;

    // hts open the writer
    fop = SeqPointer<htsFile>(hts_open(m_out.c_str(), output_format.c_str()), htsFile_delete());

    if (!fop) {
      return false;
      //throw std::runtime_error("BamWriter::Open - Cannot open output file: " + f);
    }

    return true;
  }

  BamWriter::BamWriter(int o) {

    switch(o) {
    case BAM :  output_format = "wb"; break;
    case CRAM : output_format = "wc"; break;
    case SAM :  output_format = "w"; break;
    default : throw std::invalid_argument("Invalid writer type");
    }

  }
  

bool BamWriter::WriteRecord(const BamRecord &r)
{
  if (!fop) {
    //throw std::runtime_error("BamWriter::writeAlignment - Cannot write BamRecord. Did you forget to open the Bam for writing (OpenWriteBam)?");
    return false;
  } else {
    if (sam_write1(fop.get(), hdr.get(), r.raw()) < 0)
      return false;
      //throw std::runtime_error("BamWriter::writeAlignment - Cannot write BamRecord. sam_write1 exited with < 0");      
  }

  return true;
}

std::ostream& operator<<(std::ostream& out, const BamWriter& b)
{
  if (b.fop)
    out << "Write format: " << b.fop->format.format;
  out << " Write file " << b.m_out; 
  return out;
}

  //does not return false if file not found
bool BamWriter::SetCramReference(const std::string& ref) {

  if (!fop)
    return false;

  // need to open reference for CRAM writing 
  char* fn_list = samfaipath(ref.c_str()); // eg ref = my.fa  returns my.fa.fai
  if (fn_list) {
    int status = hts_set_fai_filename(fop.get(), fn_list);
    if (status != 0) {
      fprintf(stderr, "Failed to use reference \"%s\".\n", fn_list);
      return false;
    }
  } else {
    std::cerr << "Failed to get the reference for CRAM compression" << std::endl;
    return false;
  }

  return true;
}

}
