#include "SnowTools/BamWalker.h"
#include "SnowTools/BamReader.h"

//#define DEBUG_WALKER 1

namespace SnowTools {

// set the bam region
bool BamReader::__set_region(const GenomicRegion& gp, std::shared_ptr<hts_idx_t> passed_idx) {
  
  assert(m_in.length());
  
#ifdef DEBUG_WALKER
  std::cerr << "trying to set the region for "<< m_in << " and on region " << gp << std::endl;
#endif
  
  // store the pre-loaded idx rather than re-load
  if (!idx && passed_idx) 
    idx = passed_idx;
  
  //HTS set region
  if (!idx) 
    idx = std::shared_ptr<hts_idx_t>(hts_idx_load(m_in.c_str(), HTS_FMT_BAI), idx_delete());
  
  if (gp.chr >= br->n_targets) {
    m_region_fail= true;
    std::cerr << "Failed to set region on " << gp << ". Chr ID is bigger than n_targets=" << br->n_targets << std::endl;
    return false;
  }
  
  if (!idx) {
    std::cerr << "Failed to load index file for file " << m_in << std::endl;
    std::cerr << "...suggest rebuilding index with samtools index" << std::endl;
    m_region_fail = true;
    exit(EXIT_FAILURE);
    return false;
  }
  hts_itr = std::shared_ptr<hts_itr_t>(sam_itr_queryi(idx.get(), gp.chr, gp.pos1, gp.pos2), hts_itr_delete());

  if (!hts_itr) {
    std::cerr << "Error: Failed to set region: " << gp << std::endl; 
    m_region_fail = true;
    return false;
  }

  m_region_fail = false;
  return true;
}

void BamReader::resetAll() {

  m_region_fail = false;
  m_region_idx = 0;
  m_region = GenomicRegionVector();
  m_num_reads_seen = 0;
  m_num_reads_kept = 0;

}

void BamReader::setBamReaderRegion(const GenomicRegion& g, std::shared_ptr<hts_idx_t> passed_idx)
{
  m_region.clear();
  m_region.push_back(g);
  m_region_idx = 0; // rewind it
  //__check_regions_blacklist(); // sets m_region
  if (m_region.size())
    __set_region(m_region[0], passed_idx);
  else
    m_region.push_back(GenomicRegion(-1,-1,-1));
}

void BamReader::setBamReaderRegions(const GenomicRegionVector& grv, std::shared_ptr<hts_idx_t> passed_idx) 
{
  if (grv.size() == 0)
    {
      std::cerr << "Warning: Trying to set an empty bam region"  << std::endl;
      return;
    }
  m_region = grv;
  m_region_idx = 0; // rewind it
  //__check_regions_blacklist(); // sets m_region
  if (m_region.size())
    __set_region(m_region[0], passed_idx);
  else
    m_region.push_back(GenomicRegion(-1,-1,-1));
}

bool BamReader::OpenReadBam(const std::string& bam) 
{
  m_in = bam;

  return __open_BAM_for_reading();
}

BamReader::BamReader(const std::string& in) : m_in(in)
{
  // open for reading
  if (!__open_BAM_for_reading())
    throw 20;
}

BamReader::BamReader() {}

bool BamReader::__open_BAM_for_reading()
{

  assert(m_in.length());
  
  // HTS open the reader
  fp = std::shared_ptr<BGZF>(bgzf_open(m_in.c_str(), "r"), bgzf_delete()); 
  
  if (!fp) {
    std::cerr << "Error using HTS reader on opening NGS file " << m_in << std::endl;
    exit(EXIT_FAILURE);
  } 

  br = std::shared_ptr<bam_hdr_t>(bam_hdr_read(fp.get()), bam_hdr_delete());
  
  if (!br) {
    std::cerr << "Error using HTS reader on opening NGS file " << m_in << std::endl;
    return false;
  }
  
  return true;

}

void BamReader::SetReadFilterCollection(const std::string& rules)
{

  // construct the minirules
  m_mr = ReadFilterCollection(rules, br.get());

  // check that it worked
  if (!m_mr.size()) {
    //std::cerr << "No ReadFilter were successfully parsed" << std::endl;
    //throw 20;
  }
}

void BamReader::printRuntimeMessage(const ReadCount &rc_main, const BamRead &r) const {

  char buffer[100];
  std::string posstring = AddCommas<int>(r.Position());
  sprintf (buffer, "Reading read %11s at position %2s:%-11s. Kept %11s (%2d%%) [running count across whole BAM]",  
	   rc_main.totalString().c_str(), r.ChrName().c_str(), posstring.c_str(),  
	   rc_main.keepString().c_str(), rc_main.percent());
  printf ("%s | ",buffer);

#ifndef __APPLE__
  displayRuntime(start);
#endif
  std::cerr << std::endl;
  
}

bool BamReader::GetNextRead(BamRead& r, bool& rule)
{
  
  if (m_region_fail) {
    std::cerr << "BamReader::GetNextRead - Since region failed to be set, refuse to read in reads" << std::endl;
    return false;
  }

  void* dum = 0;
  bam1_t* b = bam_init1(); 

  int32_t valid;
  if (hts_itr == 0) {
    valid = bam_read1(fp.get(), b);    
    if (valid < 0) { 

#ifdef DEBUG_WALKER
      std::cerr << "ended reading on null hts_itr" << std::endl;
#endif
      bam_destroy1(b); 
      return false;
    } 
  } else {

    valid = hts_itr_next(fp.get(), hts_itr.get(), b, dum);
  }

  if (valid <= 0) { // read not found
    do {

#ifdef DEBUG_WALKER
      std::cerr << "Failed read, trying next region. Moving counter to " << m_region_idx << " of " << m_region.size() << " FP: "  << fp << " hts_itr " << std::endl;
#endif

      // try next region, return if no others to try
      ++m_region_idx; // increment to next region
      if (m_region_idx >= m_region.size()) {
	bam_destroy1(b);
#ifdef DEBUG_WALKER
	std::cerr << "reached end of regions " << std::endl;
#endif
	return false;
      }
      // next region exists, try it
      __set_region(m_region[m_region_idx]);
      valid = hts_itr_next(fp.get(), hts_itr.get(), b, dum);
    } while (valid <= 0); // keep trying regions until works
  }
  
  r.assign(b); // = std::shared_ptr<bam1_t> (b, free_delete());

  // check if it passed the rules
  rule = m_mr.isValid(r);

  ++m_num_reads_seen;
  if (rule)
    ++m_num_reads_kept;

  // hit the limit, no more raeding
  if (m_limit >= 0 && m_num_reads_seen > m_limit)
    return false;
  if (m_keep_limit >= 0 && m_num_reads_kept > m_keep_limit)
    return false;


  return true;
}

std::string BamReader::printRegions() const {

  std::stringstream ss;
  for (auto& r : m_region)
    ss << r << std::endl;
  return(ss.str());

}

std::ostream& operator<<(std::ostream& out, const BamReader& b)
{
  out << "Read"; 
  out << ":" << b.m_in << std::endl; 
  out << b.m_mr << std::endl;
  if (b.m_region.size() && b.m_region.size() < 20) {
    out << " ------- BamReader Regions ----------" << std::endl;;
    for (auto& i : b.m_region)
      out << i << std::endl;
  } 
  else if (b.m_region.size() >= 20) {
    int wid = 0;
    for (auto& i : b.m_region)
      wid += i.width();
    out << " ------- BamReader Regions ----------" << std::endl;;
    out << " -- " << b.m_region.size() << " regions covering " << AddCommas(wid) << " bp of sequence"  << std::endl;
  }
  else 
    out << " - BamReader - Walking whole genome -" << std::endl;

  out <<   " ------------------------------------";
  return out;
}

std::string BamReader::displayReadFilterCollection() const 
{
  std::stringstream ss;
  ss << m_mr;
  return ss.str();
}


}
