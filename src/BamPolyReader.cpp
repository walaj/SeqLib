#include "SeqLib/BamPolyReader.h"


//#define DEBUG_WALKER 1

namespace SeqLib {

// set the bam region
bool BamPolyReader::__set_region(const GenomicRegion& gp) {
  
#ifdef DEBUG_WALKER
  std::cerr << "trying to set the region for "<< m_in << " and on region " << gp << std::endl;
#endif
  
  //HTS set region for all of the indicies
  for (auto& b : m_bams) {
    if ( (b.fp->format.format == 4 || b.fp->format.format == 6) && !b.idx)  // BAM (4) or CRAM (6)
      b.idx = std::shared_ptr<hts_idx_t>(sam_index_load(b.fp.get(), b.m_in.c_str()), idx_delete());
    
    if (!b.idx) {
      std::cerr << "Failed to load index file for file " << b.m_in << std::endl;
      std::cerr << "...suggest rebuilding index with samtools index" << std::endl;
    }
    
    if (gp.chr >= b.m_hdr.NumSequences()) {
      std::cerr << "Failed to set region on " << gp << ". Chr ID is bigger than n_targets=" << b.m_hdr.NumSequences() << std::endl;
      return false;
    }

    // should work for BAM or CRAM
    b.hts_itr = std::shared_ptr<hts_itr_t>(sam_itr_queryi(b.idx.get(), gp.chr, gp.pos1, gp.pos2), hts_itr_delete());

    if (!b.hts_itr)
      std::cerr << "Error: Failed to set region: " << gp << std::endl; 
  }

  return true;
}

void BamPolyReader::Reset() {

  m_region_idx = 0;
  m_region = GRC();

}

bool BamPolyReader::SetRegion(const GenomicRegion& g)
{
  m_region.clear();
  m_region.add(g);
  m_region_idx = 0; // rewind it

  if (m_region.size())
    return __set_region(m_region[0]);

  m_region.add(GenomicRegion(-1,-1,-1));
  return false;
  
}

  bool BamPolyReader::SetMultipleRegions(const GRC& grc) 
{
  if (grc.size() == 0) {
      std::cerr << "Warning: Trying to set an empty bam region"  << std::endl;
      return false;
    }

  m_region = grc;
  m_region_idx = 0; // rewind it
  
  if (m_region.size())
    return __set_region(m_region[0]);
  
  m_region.add(GenomicRegion(-1,-1,-1));
  return false;
}

  bool BamPolyReader::Open(const std::string& bam) {
    
    // id will be bam file name, unless
    // its here mutliple times, then name + random num
    std::string id = bam;
    for (auto& b : m_bams)
      if (b.id == bam)
	id = id + ":" + std::to_string(rand() % 10000);
    m_bams.push_back(_Bam(bam));
    m_bams.back().id = id;
    return m_bams.back().open_BAM_for_reading();
  }

  bool BamPolyReader::Open(const std::vector<std::string>& bams) {
    
    bool pass = true;
    for (auto& i : bams)
      pass = pass && Open(i);
    return pass;
  }
  
BamPolyReader::BamPolyReader() {}

  bool _Bam::open_BAM_for_reading()
{

  // HTS open the reader
  fp = std::shared_ptr<htsFile>(hts_open(m_in.c_str(), "r"), htsFile_delete()); 

  // open cram reference
  if (!m_cram_reference.empty()) {
    char * m_cram_reference_cstr = strdup(m_cram_reference.c_str());
    int ret = cram_load_reference(fp->fp.cram, m_cram_reference_cstr);
    free(m_cram_reference_cstr);
    if (ret < 0) 
      throw std::invalid_argument("Could not read reference genome " + m_cram_reference + " for CRAM opt");
  }

  if (!fp) 
    return false; 

  bam_hdr_t * hdr = sam_hdr_read(fp.get());
  m_hdr = BamHeader(hdr); // calls BamHeader(bam_hdr_t), makes a copy
  
  if (!m_hdr.get()) 
    return false;
  
  if (hdr)
    bam_hdr_destroy(hdr);

  return true;

}

  void BamPolyReader::SetCramReference(const std::string& ref) {
    m_cram_reference = ref;
    for (auto& b : m_bams)
      b.m_cram_reference = ref;
    
  }

bool BamPolyReader::GetNextRecord(BamRecord& r)
{

  bool found = false;

  // loop the files and load the next read
  // for the one that was emptied last
  for (auto& bam : m_bams) {

    if (!bam.empty)
      continue; // read not loaded, make it
    
    bam1_t* b = bam_init1(); 
    int32_t valid;
    
    if (bam.hts_itr == 0) {
      valid = sam_read1(bam.fp.get(), bam.m_hdr.get_(), b);    
      if (valid < 0) { 
	
#ifdef DEBUG_WALKER
	std::cerr << "ended reading on null hts_itr" << std::endl;
#endif
	goto endloop;
      } 
    } else {
      
      //changed to sam from hts_itr_next
      valid = sam_itr_next(bam.fp.get(), bam.hts_itr.get(), b);
    }
    
    if (valid < 0) { // read not found
      do {
	
#ifdef DEBUG_WALKER
	std::cerr << "Failed read, trying next region. Moving counter to " << m_region_idx << " of " << m_region.size() << " FP: "  << fp_htsfile << " hts_itr " << std::endl;
#endif
	
	// try next region, return if no others to try
	++m_region_idx; // increment to next region
	if (m_region_idx >= m_region.size()) 
	  goto endloop;
	
	// next region exists, try it
	__set_region(m_region[m_region_idx]);
	valid = sam_itr_next(bam.fp.get(), bam.hts_itr.get(), b);
      } while (valid <= 0); // keep trying regions until works
      
      
    }
    
    bam.empty = false;
    bam.next_read.assign(b); // = std::shared_ptr<bam1_t> (b, free_delete());
    continue;

    // couldn't find a valid read anywhere, move to next BAM
  endloop:
    bam_destroy1(b);
    bam.empty = true;
    continue;
  }

  // choose the one to return
  // sort based on chr and left-most alignment pos. Same as samtools
  int min_chr = INT_MAX;
  int min_pos = INT_MAX;
  size_t hit = 0, count = 0;
  for (auto& bam : m_bams) {
    ++count;

    // dont check if already marked for removal
    if (bam.empty) 
      continue;
    
    found = true;
    if (bam.next_read.ChrID() < min_chr || 
	(bam.next_read.Position() <  min_pos && bam.next_read.ChrID() == min_chr)) {
      min_pos = bam.next_read.Position();
      min_chr = bam.next_read.ChrID();
      r = bam.next_read;
      hit = count - 1;
    }
  }
  
  // mark the one we just found as empty
  if (found)
    m_bams[hit].empty = true;

  return found;
}
  
std::string BamPolyReader::PrintRegions() const {

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
