#include "SnowTools/BamWalker.h"

//#define DEBUG_WALKER 1

using namespace SnowTools;

struct idx_delete {
  void operator()(void* x) { hts_idx_destroy((hts_idx_t*)x); }
};

struct hts_itr_delete {
  void operator()(void* x) { hts_itr_destroy((hts_itr_t*)x); }
};

struct bgzf_delete {
  void operator()(void* x) { bgzf_close((BGZF*)x); }
};

struct bam_hdr_delete {
  void operator()(void* x) { bam_hdr_destroy((bam_hdr_t*)x); }
};

struct sam_write_delete {
  void operator()(void* x) { sam_close((htsFile*)x); }
};


void BamWalker::SetWriteHeader(bam_hdr_t* hdr) { 
  hdr_write = std::shared_ptr<bam_hdr_t>(hdr, bam_hdr_delete()); 
};

// set the bam region
bool BamWalker::__set_region(const GenomicRegion& gp, std::shared_ptr<hts_idx_t> passed_idx) {
  
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

void BamWalker::resetAll() {

  m_region_fail = false;
  m_region_idx = 0;
  m_region = GenomicRegionVector();
  m_num_reads_seen = 0;
  m_num_reads_kept = 0;

}

void BamWalker::setBamWalkerRegion(const GenomicRegion& g, std::shared_ptr<hts_idx_t> passed_idx)
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

void BamWalker::setBamWalkerRegions(const GenomicRegionVector& grv, std::shared_ptr<hts_idx_t> passed_idx) 
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

void BamWalker::__close_read_bam() {

  /*  if (fp)
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
  */
}

void BamWalker::__close_write_bam() 
{

  //if (fop)
  //  sam_close(fop.get());
  //fop = nullptr;
  
}
// closes the BamWriter and makes an index file
void BamWalker::makeIndex() {
  
  if (!fop)
    std::cerr << "WARNING: Trying to close write BAM when none specified" << std::endl;

  __close_write_bam();
  
  // call to htslib to buiild bai index
  bam_index_build(m_out.c_str(), 0); // 0 is "min_shift", which is 0 for bai index

}

/*bool BamWalker::OpenReadBam(FILE * stdin) 
{
  
  return __open_BAM_for_reading();
  }*/


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

bool BamWalker::__open_BAM_for_writing() 
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
  m_mr = MiniRulesCollection(rules, br.get());

  // check that it worked
  if (!m_mr.size()) {
    //std::cerr << "No MiniRules were successfully parsed" << std::endl;
    //throw 20;
  }
}

void BamWalker::printRuntimeMessage(const ReadCount &rc_main, const BamRead &r) const {

  char buffer[100];
  std::string posstring = SnowTools::AddCommas<int>(r.Position());
  sprintf (buffer, "Reading read %11s at position %2s:%-11s. Kept %11s (%2d%%) [running count across whole BAM]",  
	   rc_main.totalString().c_str(), r.ChrName().c_str(), posstring.c_str(),  
	   rc_main.keepString().c_str(), rc_main.percent());
  printf ("%s | ",buffer);
  SnowTools::displayRuntime(start);
  std::cerr << std::endl;
  
}

bool BamWalker::GetNextRead(BamRead& r, bool& rule)
{
  
  if (m_region_fail) {
    std::cerr << "BamWalker::GetNextRead - Since region failed to be set, refuse to read in reads" << std::endl;
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

void BamWalker::writeAlignment(BamRead &r)
{

  if (m_strip_all_tags)
    r.RemoveAllTags();

  for (auto& i : m_tag_list)
    r.RemoveTag(i.c_str());
    
  if (!fop)
    std::cerr << "BamWalker ERROR in writeAlignment. Did you forget to open the Bam for writing (OpenWriteBam)? Skipping write"  << std::endl;
  else
    sam_write1(fop.get(), hdr_write.get(), r.raw());
}

std::string BamWalker::printRegions() const {

  std::stringstream ss;
  for (auto& r : m_region)
    ss << r << std::endl;
  return(ss.str());

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
  if (b.m_region.size() && b.m_region.size() < 20) {
    out << " ------- BamWalker Regions ----------" << std::endl;;
    for (auto& i : b.m_region)
      out << i << std::endl;
  } 
  else if (b.m_region.size() >= 20) {
    int wid = 0;
    for (auto& i : b.m_region)
      wid += i.width();
    out << " ------- BamWalker Regions ----------" << std::endl;;
    out << " -- " << b.m_region.size() << " regions covering " << AddCommas(wid) << " bp of sequence"  << std::endl;
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

/*void BamWalker::__check_regions_blacklist() 
{

  // does not work until get complement implemented

  // check if it overlaps with blacklist
  if (blacklist.size()) 
    {
      GRC regs;
      for (auto& i : m_region)
	regs.add(i);
         //GRC out = regs.complement(blacklist);
      m_region = out.asGenomicRegionVector();
    }

    }*/

 /*void BamWalker::addBlacklist(GRC& bl) 
{
  blacklist = bl;
  }*/

void BamWalker::setCram(const std::string& out, const std::string& ref) {
  m_out = out;
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

void BamWalker::setStdout() {
  fop = std::shared_ptr<htsFile>(sam_open("-", "w"), sam_write_delete()); 
}
