#include "SnowTools/BamWalker.h"
#include "SnowTools/BamWriter.h"

//#define DEBUG_WALKER 1

namespace SnowTools {

void BamWriter::SetWriteHeader(bam_hdr_t* hdr) {
  hdr_write = std::shared_ptr<bam_hdr_t>(bam_hdr_dup(hdr), bam_hdr_delete()); 
};

void BamWriter::makeIndex() {
  
  if (!fop)
    std::cerr << "WARNING: Trying to close write BAM when none specified" << std::endl;

  //__close_write_bam();
  
  // call to htslib to build bai index
  bam_index_build(m_out.c_str(), 0); // 0 is "min_shift", which is 0 for bai index

}
bool BamWriter::OpenWriteBam(const std::string& bam) 
{
  m_out = bam;

  return __open_BAM_for_writing();
}


BamWriter::BamWriter(const std::string& in) : m_in(in)
{
  // open for reading
  if (!__open_BAM_for_writing())
    throw 20;
}

bool BamWriter::__open_BAM_for_writing() 
{

  // hts open the writer
  if (!fop) { // default is bam. if already set by flag, don't reopen
    assert(m_out.length());
    fop = std::shared_ptr<htsFile>(sam_open(m_out.c_str(), "wb"), sam_write_delete());
    m_print_header = true;
  }

  if (!fop) {
    std::cerr << "Error: Cannot open BAM for writing " << m_out << std::endl;
    return false;
  }

  // if no write header, set as read
  if (!hdr_write)
    hdr_write = br;

  // hts write the header
  if (m_print_header) {
    sam_hdr_write(fop.get(), hdr_write.get());      
  }

  return true;

}

void BamWriter::writeAlignment(BamRead &r)
{

  //for (auto& i : m_tag_list)
  //  r.RemoveTag(i.c_str());
    
  if (!fop)
    std::cerr << "BamWriter ERROR in writeAlignment. Did you forget to open the Bam for writing (OpenWriteBam)? Skipping write"  << std::endl;
  else
    sam_write1(fop.get(), hdr_write.get(), r.raw());
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
