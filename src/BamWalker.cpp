#include "SnowTools/BamWalker.h"

#define DEBUG_WALKER 1

using namespace SnowTools;

// set the bam region
bool BamWalker::__set_region(const GenomicRegion& gp) {

  assert(m_in.length());

#ifdef DEBUG_WALKER
    std::cout << "trying to set the region for "<< m_in << " and on region " << gp << std::endl;
#endif

  //HTS set region
  if (!idx) {
    idx = hts_idx_load(m_in.c_str(), HTS_FMT_BAI);
#ifdef DEBUG_WALKER
    std::cout << "loading the index for "<< m_in << " and has value " << idx << std::endl;
#endif
  }
  hts_itr = sam_itr_queryi(idx, gp.chr, gp.pos1, gp.pos2);
  if (!hts_itr) {
    std::cerr << "Error: Failed to set region: " << gp << endl; 
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

// closes the BamWriter and makes an index file
/*void BamWalker::MakeIndex() {

#ifdef HAVE_BAMTOOLS
  m_writer->Close();
  
  // open the file 
  BamReader reader;
  if (!reader.Open(m_out)) {
    cerr << "Error: Could not open the output BAM to create index " << m_out << endl;
    exit(EXIT_FAILURE);
  }

  // create the index
  if (!reader.CreateIndex()) {
    cerr << "Error: Could not create the output BAM index for " << m_out << endl;
    exit(EXIT_FAILURE);
  }

  reader.Close();
#endif
}*/

bool BamWalker::OpenReadBam(const std::string& bam) 
{
  m_in = bam;

  return __open_BAM_for_reading();
}

bool BamWalker::OpenWriteBam(const std::string& bam) 
{
  m_in = bam;

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
    cerr << "Error using HTS reader on opening " << m_in << endl;
    exit(EXIT_FAILURE);
  }
  br = bam_hdr_read(fp);
  
  if (!br) {
    cerr << "Error using HTS reader on opening " << m_in << endl;
    return false;
    //exit(EXIT_FAILURE);
  }
  
  return true;

}

bool BamWalker::__open_BAM_for_writing() 
{

  assert(m_out.length());

  // hts open the writer
  fop = sam_open(m_out.c_str(), "wb");
  if (!fop) {
    cerr << "Error: Cannot open BAM for writing " << m_out << endl;
    return false;
    //exit(EXIT_FAILURE);
  }

  // hts write the header
  sam_hdr_write(fop, br);

  return true;

}

// this is also opens the files for reading/writing and checks
// that they are readable/writable
BamWalker::BamWalker(const std::string& in, const std::string& out) : m_in(in), m_out(out)
{

  // open for reading
  if (!__open_BAM_for_reading())
    throw 20;

  // open for writing
  if (!__open_BAM_for_writing())
    throw 20;

}

void BamWalker::SetMiniRulesCollection(const std::string& rules)
{
  // construct the minirules
  m_mr = MiniRulesCollection(rules);

  // check that it worked
  if (!m_mr.size()) {
    std::cerr << "No MiniRules were successfully parsed" << std::endl;
    throw 20;
  }
}

void BamWalker::printRuntimeMessage(const ReadCount &rc_main, const Read &r) const {

  char buffer[100];
  string posstring = SnowTools::AddCommas<int>(r_pos(r));
  sprintf (buffer, "Reading read %11s at position %2s:%-11s. Kept %11s (%2d%%) [running count across whole BAM]",  
	   rc_main.totalString().c_str(), GenomicRegion::chrToString(r_id(r)).c_str(), posstring.c_str(),  
	   rc_main.keepString().c_str(), rc_main.percent());
  printf ("%s | ",buffer);
  SnowTools::displayRuntime(start);
  cout << endl;
  
}

void BamWalker::printRuleCounts(unordered_map<string, size_t> &rm) const {

  size_t total = 0;
  for (auto& i : rm)
    total += i.second;
  for (auto& i : rm) {
    cout << "  " << i.first << ":" << i.second << "(" << SnowTools::percentCalc<size_t>(i.second, total) << "%)" << endl;
  }
  
}

bool BamWalker::GetNextRead(Read& r, std::string& rule)
{
  
  void* dum = 0;
  bam1_t* b = bam_init1(); 
  if (hts_itr == 0) { 
    if (bam_read1(fp, b) < 0) { 

#ifdef DEBUG_WALKER
      std::cout << "ended reading on null hts_itr" << std::endl;
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
      std::cout << "Failed read, trying next region. FP: "  << fp << " hts_itr " << 
#endif

      // try next region, return if no others to try
      ++m_region_idx; // increment to next region
      if (m_region_idx >= m_region.size()) {
	bam_destroy1(b);
#ifdef DEBUG_WALKER
	std::cout << "reached end of regions " << std::endl;
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
  sam_write1(fop, br, r.get());
}

std::ostream& SnowTools::operator<<(std::ostream& out, const BamWalker& b)
{
  out << " -- In Bam:  " << b.m_in << std::endl;
  out << " -- Out Bam: " << b.m_out << std::endl;
  out << b.m_mr << std::endl;
  if (b.m_region.size()) {
    out << " ------- Regions ---------" << std::endl;;
    for (auto& i : b.m_region)
      out << i << std::endl;
  }
  else 
    out << " -- Walking whole genome" << std::endl;

  return out;
}
