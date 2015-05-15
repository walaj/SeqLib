#include "SnowTools/BamWalker.h"

//#define DEBUG_WALKER 1

using namespace SnowTools;

// set the bam region
bool BamWalker::__set_region(const GenomicRegion& gp) {

  assert(m_in.length());

#ifdef DEBUG_WALKER
    std::cerr << "trying to set the region for "<< m_in << " and on region " << gp << std::endl;
#endif

  //HTS set region
  if (!idx) {
    idx = hts_idx_load(m_in.c_str(), HTS_FMT_BAI);
#ifdef DEBUG_WALKER
    std::cerr << "loading the index for "<< m_in << " and has value " << idx << std::endl;
#endif
  }
  hts_itr = sam_itr_queryi(idx, gp.chr, gp.pos1, gp.pos2);
  if (!hts_itr) {
    std::cerr << "Error: Failed to set region: " << gp << std::endl; 
    return false;
  }

  return true;
}

void BamWalker::setBamWalkerRegion(const GenomicRegion& g)
{
  m_region.clear();
  m_region.push_back(g);
  m_region_idx = 0; // rewind it
  __set_region(g);
}

void BamWalker::setBamWalkerRegions(const GenomicRegionVector& grv) 
{
  assert(grv.size());
  m_region = grv;
  m_region_idx = 0; // rewind it
  //m_region.rewind();
  __set_region(grv[0]);
}

void BamWalker::__close_read_bam() {

  if (fp)
    bgzf_close(fp);
  if (br)
    bam_hdr_destroy(br);
  if (hts_itr)
    hts_itr_destroy(hts_itr);
  if (idx)
    hts_idx_destroy(idx);

  fp = NULL;
  br = NULL;
  hts_itr = NULL;
  idx = NULL;

}

void BamWalker::__close_write_bam() 
{
  if (fop)
    sam_close(fop);
  fop = NULL;
  
}
// closes the BamWriter and makes an index file
void BamWalker::MakeIndex() {
  
  __close_write_bam();
  
  // call to htslib to buiild bai index
  bam_index_build(m_out.c_str(), 0); // 0 is "min_shift", which is 0 for bai index

}

bool BamWalker::OpenReadBam(const std::string& bam) 
{
  m_in = bam;

  return __open_BAM_for_reading();
}

bool BamWalker::OpenWriteBam(const std::string& bam) 
{
  m_out = bam;

  return __open_BAM_for_writing();
}


BamWalker::BamWalker(const std::string& in) : m_in(in)
{
  // open for reading
  if (!__open_BAM_for_reading())
    throw 20;

}

bool BamWalker::__open_BAM_for_reading()
{

  assert(m_in.length());
  
  // HTS open the reader
  const char rflag = 'r';
  fp = bgzf_open(m_in.c_str(), &rflag); 
  
  if (!fp) {
    std::cerr << "Error using HTS reader on opening " << m_in << std::endl;
    exit(EXIT_FAILURE);
  }
  br = bam_hdr_read(fp);
  
  if (!br) {
    std::cerr << "Error using HTS reader on opening " << m_in << std::endl;
    return false;
    //exit(EXIT_FAILURE);
  }
  
  return true;

}

bool BamWalker::__open_BAM_for_writing() 
{

  // hts open the writer
  if (!fop) { // default is bam. if already set by flag, don't reopen
    assert(m_out.length());
    fop = sam_open(m_out.c_str(), "wb");
    m_print_header = true;
  }
  if (!fop) {
    std::cerr << "Error: Cannot open BAM for writing " << m_out << std::endl;
    return false;
    //exit(EXIT_FAILURE);
  }

  // hts write the header
  if (m_print_header)
    sam_hdr_write(fop, br);

  return true;

}

// this is also opens the files for reading/writing and checks
// that they are readable/writable
/*BamWalker::BamWalker(const std::string& in, const std::string& out) : m_in(in), m_out(out)
{

  // open for reading
  if (!__open_BAM_for_reading())
    throw 20;

  // open for writing
  if (out.length())
    if (!__open_BAM_for_writing())
      throw 20;

      }*/

void BamWalker::SetMiniRulesCollection(const std::string& rules)
{

  // construct the minirules
  m_mr = MiniRulesCollection(rules, br);

  // check that it worked
  if (!m_mr.size()) {
    std::cerr << "No MiniRules were successfully parsed" << std::endl;
    //throw 20;
  }
}

void BamWalker::printRuntimeMessage(const ReadCount &rc_main, const Read &r) const {

  char buffer[100];
  std::string posstring = SnowTools::AddCommas<int>(r_pos(r));
  sprintf (buffer, "Reading read %11s at position %2s:%-11s. Kept %11s (%2d%%) [running count across whole BAM]",  
	   rc_main.totalString().c_str(), GenomicRegion::chrToString(r_id(r)).c_str(), posstring.c_str(),  
	   rc_main.keepString().c_str(), rc_main.percent());
  printf ("%s | ",buffer);
  SnowTools::displayRuntime(start);
  std::cerr << std::endl;
  
}

bool BamWalker::GetNextRead(Read& r, bool& rule)
{
  
  void* dum = 0;
  bam1_t* b = bam_init1(); 
  if (hts_itr == 0) { 
    if (bam_read1(fp, b) < 0) { 

#ifdef DEBUG_WALKER
      std::cerr << "ended reading on null hts_itr" << std::endl;
#endif
      bam_destroy1(b); 
      return false;
    } 
  } 

  int32_t valid;
  if (hts_itr == 0)
    valid = bam_read1(fp, b);
  else
    valid = hts_itr_next(fp, hts_itr, b, dum);

  if (valid <= 0) { // read not found
    do {

#ifdef DEBUG_WALKER
      std::cerr << "Failed read, trying next region. FP: "  << fp << " hts_itr " << 
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
      valid = hts_itr_next(fp, hts_itr, b, dum);
    } while (valid <= 0); // keep trying regions until works
  }
  
  r = std::shared_ptr<bam1_t> (b, free_delete());

  // check if it passed the rules
  rule = m_mr.isValid(r);

  return true;
}

void BamWalker::WriteAlignment(Read &r)
{
  if (m_strip_all_tags)
    removeAllTags(r);
  for (auto& i : m_tag_list)
    r_remove_tag(r,i.c_str());
    
  sam_write1(fop, br, r.get());
}

std::ostream& SnowTools::operator<<(std::ostream& out, const BamWalker& b)
{
  out << " -- In Bam:  " << b.m_in << std::endl;
  out << " -- Out Bam: " << b.m_out << std::endl;
  out << " -- Output format: ";
  if (!b.fop)
    out << "BAM" << std::endl;
  else if (b.fop->format.format == 4)
    out << "BAM" << std::endl;
  else if (b.fop->format.format == 6)
    out << "CRAM" << std::endl;
  else if (b.fop->format.format == text_format)
    out << "SAM" << std::endl;

  out << b.m_mr << std::endl;
  if (b.m_region.size()) {
    out << " ------- BamWalker Regions ----------" << std::endl;;
    for (auto& i : b.m_region)
      out << i << std::endl;
  }
  else 
    out << " - BamWalker - Walking whole genome -" << std::endl;

  out <<   " ------------------------------------";
  return out;
}

void BamWalker::setStripTags(const std::string& list)
{
  std::istringstream iss(list);
  std::string val;
  while(std::getline(iss, val, ',')) 
    {
      m_tag_list.push_back(val);
    }
}

std::string BamWalker::displayMiniRulesCollection() const 
{
  std::stringstream ss;
  ss << m_mr;
  return ss.str();
}
