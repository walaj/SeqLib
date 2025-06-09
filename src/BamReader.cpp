#include "SeqLib/BamReader.h"
#include <iostream>
#include <htslib/sam.h>
#include <htslib/hts.h>
#include <cstdlib>

namespace SeqLib {

bool BamReader::Open(const std::string& path) {

  // store the string path
  path_ = path;
  
  // Open the BAM/SAM/CRAM file ("-" for stdin)
  fp_.reset(hts_open(path.c_str(), "r"));
  if (!fp_) {
    std::cerr << "BamReader::Open - failed to open '" << path << "'" << std::endl;
    return false;
  }

  // Read and wrap the header
  bam_hdr_t* r = sam_hdr_read(fp_.get());
  if (!r) {
    std::cerr << "BamReader::Open - failed to read header from '" << path << "'" << std::endl;
    fp_.reset();
    return false;
  }
  hdr_ = BamHeader(r); // copies it over to BamHeader
  bam_hdr_destroy(r);  // deletes it since already just copied it over  

  // Load index for random access
  idx_.reset(sam_index_load(fp_.get(), path.c_str()));
  // iterator remains null until SetRegion/SetRegions

  // Clear any previous region state
  regions_.clear();
  region_idx_ = 0;
  itr_.reset();
  
  return true;
}

void BamReader::Close() {
  // drop any active iterator and index
  itr_.reset();
  idx_.reset();
  // drop the header and file handle
  fp_.reset();
  hdr_ = BamHeader();
  // clear region list and reset position
  regions_.clear();
  region_idx_ = 0;
}

void BamReader::Reset() {

  // easier to close and reset
  this->Close();
  this->Open(path_);
  region_idx_ = 0;
}  

bool BamReader::SetRegion(const GenomicRegion& region) {
  if (!fp_ || !idx_)
    return false;

  // clear any existing regions, then set this one
  regions_.clear();
  regions_.add(region);
  region_idx_ = 0;

  // arm the iterator on the requested region
  itr_.reset(sam_itr_queryi(
    idx_.get(),
    region.chr,
    region.pos1,
    region.pos2
  ));

  return static_cast<bool>(itr_);
}

bool BamReader::SetRegions(const GRC& regions) {
  if (!fp_ || !idx_ || !regions.size())
    return false;

  // copy all desired regions
  regions_ = regions;
  region_idx_ = 0;

  // arm iterator on the first region
  const auto& rg = regions_[0];
  itr_.reset(sam_itr_queryi(
    idx_.get(),
    rg.chr,
    rg.pos1,
    rg.pos2
  ));

  return static_cast<bool>(itr_);
}

std::optional<BamRecord> BamReader::Next() {
  if (!fp_) 
    return std::nullopt;

  // Allocate a fresh bam1_t*
  bam1_t* raw = bam_init1();
  int ret = -1;

  if (itr_) {
    // Read from the current region iterator
    ret = sam_itr_next(fp_.get(), itr_.get(), raw);
    if (ret < 0) {
      // Move on to the next region, if any
      ++region_idx_;
      while (region_idx_ < regions_.size()) {
        const auto& rg = regions_[region_idx_];
        itr_.reset(sam_itr_queryi(
          idx_.get(),
          rg.chr,
          rg.pos1,
          rg.pos2
        ));
        if (!itr_) { 
          ++region_idx_; 
          continue; 
        }
        // Try reading from the newly armed iterator
        ret = sam_itr_next(fp_.get(), itr_.get(), raw);
        if (ret > 0) 
          break;
        ++region_idx_;
      }
    }
  }
  else {
    // No region restriction: read sequentially
    ret = sam_read1(fp_.get(), hdr_.get_(), raw);
  }

  // If nothing was read (EOF or error), clean up and signal end
  if (ret < 0) {
    bam_destroy1(raw);
    return std::nullopt;
  }

  // Wrap the raw pointer in a BamRecord (takes ownership) and return it
  return BamRecord(raw);
}

  const BamHeader& BamReader::Header() const {
    if (!fp_) 
      throw std::runtime_error("BamReader::Header() called before Open()");
    return hdr_;
  }
  
  void BamReader::SetCramReference(const std::string& cram) {
    
    namespace fs = std::filesystem;
    if (!fs::exists(cram) || !fs::is_regular_file(cram) || !(std::ifstream(cram).good())) {
      throw std::runtime_error("File '" + cram + "' does not exist or is not readable.");
    }
    
    cram_reference_ = cram;
  }

  std::string BamReader::PrintRegions() const {
    std::ostringstream ss;
    for (const auto& region : regions_) {
      ss << region << '\n';
    }
    return ss.str();
  }

  std::ostream& operator<<(std::ostream& out, const BamReader& b) {
    out << ": " << b.path_ << '\n';
    
    const auto& regions = b.regions_;
    if (!regions.IsEmpty() && regions.size() < 20) {
      out << " ------- BamReader Regions ----------\n";
      for (const auto& region : regions) {
	out << region << '\n';
        }
    } else if (regions.size() >= 20) {
      int total_width = 0;
      for (const auto& region : regions) {
	total_width += region.Width();
      }
      out << " ------- BamReader Regions ----------\n";
      out << " -- " << regions.size() << " regions covering "
	  << AddCommas(total_width) << " bp of sequence\n";
    } else {
      out << " - BamReader - Walking whole genome -\n";
    }
    
    out << " ------------------------------------";
    return out;
  }
  
}
