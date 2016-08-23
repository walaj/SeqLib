#include "SeqLib/BamWalker.h"
#include "SeqLib/BamReader.h"

extern "C" {
#include "htslib/cram/cram.h"
#include "htslib/cram/cram_io.h"
}

//#define DEBUG_WALKER 1

namespace SeqLib {

  void BamReader::SetReadFilterCollection(const ReadFilterCollection& mr) {
    m_mr = mr;
  }

// set the bam region
bool BamReader::__set_region(const GenomicRegion& gp, std::shared_ptr<hts_idx_t> passed_idx) {
  
#ifdef DEBUG_WALKER
  std::cerr << "trying to set the region for "<< m_in << " and on region " << gp << std::endl;
#endif
  
  // store the pre-loaded idx rather than re-load
  if (!idx && passed_idx) 
    idx = passed_idx;
  
  //HTS set region
  if (!idx) {
    idx = std::shared_ptr<hts_idx_t>(hts_idx_load(m_in.c_str(), HTS_FMT_BAI), idx_delete());
  }
  
  if (gp.chr >= m_hdr.NumSequences()) {
    std::cerr << "Failed to set region on " << gp << ". Chr ID is bigger than n_targets=" << m_hdr.NumSequences() << std::endl;
    return false;
  }
  
  if (!idx) {
    std::cerr << "Failed to load index file for file " << m_in << std::endl;
    std::cerr << "...suggest rebuilding index with samtools index" << std::endl;
    return false;
  }
  hts_itr = std::shared_ptr<hts_itr_t>(sam_itr_queryi(idx.get(), gp.chr, gp.pos1, gp.pos2), hts_itr_delete());

  if (!hts_itr) {
    std::cerr << "Error: Failed to set region: " << gp << std::endl; 
    return false;
  }

  return true;
}

void BamReader::resetAll() {

  m_region_idx = 0;
  m_region = GenomicRegionVector();

}

bool BamReader::setBamReaderRegion(const GenomicRegion& g, std::shared_ptr<hts_idx_t> passed_idx)
{
  m_region.clear();
  m_region.push_back(g);
  m_region_idx = 0; // rewind it

  if (m_region.size())
    return __set_region(m_region[0], passed_idx);

  m_region.push_back(GenomicRegion(-1,-1,-1));
  return false;
  
}

bool BamReader::setBamReaderRegions(const GenomicRegionVector& grv, std::shared_ptr<hts_idx_t> passed_idx) 
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
    return __set_region(m_region[0], passed_idx);

  m_region.push_back(GenomicRegion(-1,-1,-1));
  return false;
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
    throw std::runtime_error("BamReader: Cannot read file: " + m_in);
}

BamReader::BamReader() {}

bool BamReader::__open_BAM_for_reading()
{

  // HTS open the reader
  //fp = (m_in == "-" ) ? std::shared_ptr<BGZF>(bgzf_fdopen(fileno(stdin), "r")) : std::shared_ptr<BGZF>(bgzf_open(m_in.c_str(), "r"), bgzf_delete()); 
  //const char format = 'b';
  fp_htsfile = std::shared_ptr<htsFile>(hts_open(m_in.c_str(), "r")); 

  // open cram reference
  if (!m_cram_reference.empty()) {
    char * m_cram_reference_cstr = strdup(m_cram_reference.c_str());
    int ret = cram_load_reference(fp_htsfile->fp.cram, m_cram_reference_cstr);
    free(m_cram_reference_cstr);
    if (ret < 0) 
      throw std::invalid_argument("Could not read reference genome " + m_cram_reference + " for CRAM opt");
  }
  
  if (!fp_htsfile) 
    return false; 

  bam_hdr_t * hdr = sam_hdr_read(fp_htsfile.get());
  m_hdr = BamHeader(hdr); // calls BamHeader(bam_hdr_t), makes a copy

  if (!m_hdr.get()) 
    return false;
  
  if (hdr)
    bam_hdr_destroy(hdr);

  return true;

}

void BamReader::SetReadFilterCollection(const std::string& rules)
{

  // construct the minirules
  m_mr = ReadFilterCollection(rules, m_hdr.get());

  // check that it worked
  if (!m_mr.size()) {
    //std::cerr << "No ReadFilter were successfully parsed" << std::endl;
    //throw 20;
  }
}

void BamReader::printRuntimeMessage(const ReadCount &rc_main, const BamRecord &r) const {

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

  bool BamReader::SetCramReference(const std::string& ref) {
    
    m_cram_reference = ref;

  }

bool BamReader::GetNextRead(BamRecord& r, bool& rule)
{
  
  void* dum = 0;
  bam1_t* b = bam_init1(); 

  int32_t valid;
  if (hts_itr == 0) {
    valid = sam_read1(fp_htsfile.get(), m_hdr.get_(), b);    
    if (valid < 0) { 

#ifdef DEBUG_WALKER
      std::cerr << "ended reading on null hts_itr" << std::endl;
#endif
      bam_destroy1(b); 
      return false;
    } 
  } else {
    
    //changed to sam from hts_itr_next
    valid = sam_itr_next(fp_htsfile.get(), hts_itr.get(), b);
  }

  if (valid < 0) { // read not found
    do {

#ifdef DEBUG_WALKER
      std::cerr << "Failed read, trying next region. Moving counter to " << m_region_idx << " of " << m_region.size() << " FP: "  << fp_htsfile << " hts_itr " << std::endl;
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
      valid = sam_itr_next(fp_htsfile.get(), hts_itr.get(), b);
    } while (valid <= 0); // keep trying regions until works
  }
  
  r.assign(b); // = std::shared_ptr<bam1_t> (b, free_delete());

  // check if it passed the rules
  rule = m_mr.isValid(r);

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
