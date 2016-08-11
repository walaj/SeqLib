#include "SnowTools/BamWalker.h"
#include "SnowTools/BamWriter.h"

//#define DEBUG_WALKER 1

namespace SnowTools {

void BamWriter::SetWriteHeader(bam_hdr_t* hdr) {
  hdr_write = std::shared_ptr<bam_hdr_t>(bam_hdr_dup(hdr), bam_hdr_delete()); 
};

  void BamWriter::CloseBam() {

    if (!fop)
      throw std::runtime_error("BamWriter::CloseBam() - Trying to close BAM that is already closed or never opened");

    sam_close(fop.get());

    // clear the shared_ptr to output BAM
    fop = nullptr;
  }

void BamWriter::makeIndex() const {
  
  // throw an error if BAM is not already closed
  if (fop)
    throw std::runtime_error("BamWriter::makeIndex - Trying to index open BAM. Close first with CloseBam()");
  
  // call to htslib to build bai index
  if (bam_index_build(m_out.c_str(), 0) < 0) // 0 is "min_shift", which is 0 for bai index
    throw std::runtime_error("BamWriter::makeIndex - Failed to create index");
  

}
void BamWriter::OpenWriteBam(const std::string& bam) 
{
  m_out = bam;

  __open_BAM_for_writing();

  return;
}


BamWriter::BamWriter(const std::string& f) : m_out(f)
{
  __open_BAM_for_writing();
}


void BamWriter::__open_BAM_for_writing() 
{

  // hts open the writer
  if (!fop) { // default is bam. if already set by flag, don't reopen
    assert(m_out.length());
    fop = std::shared_ptr<htsFile>(sam_open(m_out.c_str(), "wb"), sam_write_delete());
    m_print_header = true;
  }

  if (!fop) 
    throw std::runtime_error("BamWriter - Cannot open BAM for writing");

  // if no write header, set as read
  if (!hdr_write)
    hdr_write = br;

  // hts write the header
  if (m_print_header) {
    if (sam_hdr_write(fop.get(), hdr_write.get()) < 0) {
      throw std::runtime_error("BamWriter - Cannot write header. sam_hdr_write exited with < 0");
    }
  }


}

void BamWriter::writeAlignment(BamRecord &r)
{
  if (!fop) {
    throw std::runtime_error("BamWriter::writeAlignment - Cannot write BamRecord. Did you forget to open the Bam for writing (OpenWriteBam)?");
  } else {
    if (sam_write1(fop.get(), hdr_write.get(), r.raw()) < 0)
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
