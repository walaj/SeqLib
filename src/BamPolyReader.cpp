#include "SeqLib/BamPolyReader.h"

//#define DEBUG_WALKER 1

namespace SeqLib {

  void BamPolyReader::SetReadFilterCollection(const ReadFilterCollection& mr) {
    m_mr = mr;
  }

// set the bam region
bool BamPolyReader::__set_region(const GenomicRegion& gp) {
  
#ifdef DEBUG_WALKER
  std::cerr << "trying to set the region for "<< m_in << " and on region " << gp << std::endl;
#endif
  
  //HTS set region for all of the indicies
  for (auto& b : m_bams) {
    if (!b.idx) {
      b.idx = std::shared_ptr<hts_idx_t>(hts_idx_load(b.m_in.c_str(), HTS_FMT_BAI), idx_delete());
    }
    if (!b.idx) {
      std::cerr << "Failed to load index file for file " << b.m_in << std::endl;
      std::cerr << "...suggest rebuilding index with samtools index" << std::endl;
    }
    
    if (gp.chr >= b.m_hdr.NumSequences()) {
      std::cerr << "Failed to set region on " << gp << ". Chr ID is bigger than n_targets=" << b.m_hdr.NumSequences() << std::endl;
      return false;
    }
 
    b.hts_itr = std::shared_ptr<hts_itr_t>(sam_itr_queryi(b.idx.get(), gp.chr, gp.pos1, gp.pos2), hts_itr_delete());
    if (!b.hts_itr)
      std::cerr << "Error: Failed to set region: " << gp << std::endl; 
  }

  return true;
}

void BamPolyReader::resetAll() {

  m_region_idx = 0;
  m_region = GenomicRegionVector();

}

bool BamPolyReader::setBamReaderRegion(const GenomicRegion& g)
{
  m_region.clear();
  m_region.push_back(g);
  m_region_idx = 0; // rewind it

  if (m_region.size())
    return __set_region(m_region[0]);

  m_region.push_back(GenomicRegion(-1,-1,-1));
  return false;
  
}

  bool BamPolyReader::setBamReaderRegions(const GenomicRegionVector& grv) 
{
  if (grv.size() == 0)
    {
      std::cerr << "Warning: Trying to set an empty bam region"  << std::endl;
      return false;
    }
  m_region = grv;
  m_region_idx = 0; // rewind it
  //__check_regions_blacklist(); // sets m_region
  if (m_region.size())
    return __set_region(m_region[0]);

  m_region.push_back(GenomicRegion(-1,-1,-1));
  return false;
}

  bool BamPolyReader::OpenReadBam(const std::string& bam) {
    m_bams.push_back(_Bam(bam));
    return m_bams.back().open_BAM_for_reading();
  }
  
BamPolyReader::BamPolyReader() {}

bool _Bam::open_BAM_for_reading()
{

  // HTS open the reader
  fp = m_in == "-" ? std::shared_ptr<BGZF>(bgzf_fdopen(fileno(stdin), "r")) : std::shared_ptr<BGZF>(bgzf_open(m_in.c_str(), "r"), bgzf_delete()); 
  
  if (!fp) 
    return false; 

  //br = std::shared_ptr<bam_hdr_t>(bam_hdr_read(fp.get()), bam_hdr_delete()g);
  bam_hdr_t * hdr = bam_hdr_read(fp.get());
  m_hdr = BamHeader(hdr); // calls BamHeader(bam_hdr_t), makes a copy
  
  if (!m_hdr.get()) 
    return false;
  
  if (hdr)
    bam_hdr_destroy(hdr);

  return true;

}

bool BamPolyReader::GetNextRead(BamRecord& r, bool& rule)
{
  
  return true;
}

std::string BamPolyReader::printRegions() const {

  std::stringstream ss;
  for (auto& r : m_region)
    ss << r << std::endl;
  return(ss.str());

}

std::ostream& operator<<(std::ostream& out, const BamPolyReader& b)
{
  for (auto& bam : b.m_bams) 
    out << ":" << bam.m_in << std::endl; 

  if (b.m_region.size() && b.m_region.size() < 20) {
    out << " ------- BamPolyReader Regions ----------" << std::endl;;
    for (auto& i : b.m_region)
      out << i << std::endl;
  } 
  else if (b.m_region.size() >= 20) {
    int wid = 0;
    for (auto& i : b.m_region)
      wid += i.width();
    out << " ------- BamPolyReader Regions ----------" << std::endl;;
    out << " -- " << b.m_region.size() << " regions covering " << AddCommas(wid) << " bp of sequence"  << std::endl;
  }
  else 
    out << " - BamPolyReader - Walking whole genome -" << std::endl;

  out <<   " ------------------------------------";
  return out;
}

}
