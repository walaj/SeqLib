#include "SnowTools/BamRead.h"

#include <cassert>

namespace SnowTools {

  struct free_delete {
    void operator()(void* x) { bam_destroy1((bam1_t*)x); }
  };
  
  void BamRead::init() {
    bam1_t* f = bam_init1();
    b = std::shared_ptr<bam1_t>(f, free_delete());
  }

  void BamRead::assign(bam1_t* a) { 
    b = std::shared_ptr<bam1_t>(a, free_delete()); 
  }

  GenomicRegion BamRead::asGenomicRegion() const {
    return GenomicRegion(b->core.tid, b->core.pos, PositionEnd());
  }

  std::string BamRead::Sequence() const {
    uint8_t * p = bam_get_seq(b);
    std::string out(b->core.l_qseq, 'N');
    for (int32_t i = 0; i < b->core.l_qseq; ++i) 
      out[i] = BASES[bam_seqi(p,i)];
    return out;
    
  }

  void BamRead::SmartAddTag(const std::string& tag, const std::string& val)
  {
    // get the old tag
    assert(tag.length());
    assert(val.length());
    std::string tmp = GetZTag(tag);
    if (!tmp.length()) 
      {
	AddZTag(tag, val);
	return;
      }
    
    // append the tag
    tmp += "x"  + val;
    
    // remove the old tag
    RemoveTag(tag.c_str());
    
    // add the new one
    assert(tmp.length());
    AddZTag(tag, tmp);
  }

  void BamRead::clearSeqQualAndTags() {

    int new_size = b->core.l_qname + ((b)->core.n_cigar<<2);// + 1; ///* 0xff seq */ + 1 /* 0xff qual */;
    b->data = (uint8_t*)realloc(b->data, new_size);
    b->l_data = new_size;
    b->core.l_qseq = 0;
  }

  void BamRead::SetSequence(const std::string& seq) {

    int new_size = b->l_data - ((b->core.l_qseq+1)>>1) - b->core.l_qseq + ((seq.length()+1)>>1) + seq.length();    
    int old_aux_spot = (b->core.n_cigar<<2) + b->core.l_qname + ((b->core.l_qseq + 1)>>1) + b->core.l_qseq;
    int old_aux_len = bam_get_l_aux(b); //(b->core.n_cigar<<2) + b->core.l_qname + ((b->core.l_qseq + 1)>>1) + b->core.l_qseq;

    // copy out all the old data
    uint8_t* oldd = (uint8_t*)malloc(b->l_data);
    memcpy(oldd, b->data, b->l_data);
    
    // clear out the old data and alloc the new amount
    free(b->data);
    b->data = (uint8_t*)calloc(new_size, sizeof(uint8_t)); 
    
    // add back the qname and cigar
    memcpy(b->data, oldd, b->core.l_qname + (b->core.n_cigar<<2));

    // update the sizes
    // >>1 shift is because only 4 bits needed per ATCGN base
    b->l_data = new_size; //b->l_data - ((b->core.l_qseq + 1)>>1) - b->core.l_qseq + ((seq.length()+1)>>1) + seq.length();
    b->core.l_qseq = seq.length();
    
    // allocate the sequence
    uint8_t* m_bases = b->data + b->core.l_qname + (b->core.n_cigar<<2);
    int slen = seq.length();

    for (int i = 0; i < slen; ++i) {
	
      // bad idea but works for now
      uint8_t base = 15;
      if (seq.at(i) == 'A')
	base = 1;
      else if (seq.at(i) == 'C')
	base = 2;
      else if (seq.at(i) == 'G')
	base = 4;
      else if (seq.at(i) == 'T')
	base = 8;
      
      m_bases[i >> 1] &= ~(0xF << ((~i & 1) << 2));   ///< zero out previous 4-bit base encoding
      m_bases[i >> 1] |= base << ((~i & 1) << 2);  ///< insert new 4-bit base encoding
    }

    // add in a NULL qual
    uint8_t* s = bam_get_qual(b);
    s[0] = 0xff;

    // add the aux data
    uint8_t* t = bam_get_aux(b);
    memcpy(t, oldd + old_aux_spot, old_aux_len);

    // reset the max size
    b->m_data = b->l_data;
    
  }
  
  void BamRead::SetQname(const std::string& n)
  {
    // copy out the non-qname data
    size_t nonq_len = b->l_data - b->core.l_qname;
    uint8_t* nonq = (uint8_t*)malloc(nonq_len);
    memcpy(nonq, b->data + b->core.l_qname, nonq_len);

    // clear the old data and alloc the new amount 
    free(b->data);
    b->data = (uint8_t*)calloc(nonq_len + n.length() + 1, 1);
    
    // add in the new qnamev
    memcpy(b->data, (uint8_t*)n.c_str(), n.length() + 1); // +1 for \0

    // update the sizes
    b->l_data = b->l_data - b->core.l_qname + n.length() + 1;
    b->core.l_qname = n.length() + 1;    
    
    // copy over the old data
    memcpy(b->data + b->core.l_qname, nonq, nonq_len);
    free(nonq);

    // reset the max size
    b->m_data = b->l_data;
  }

  /*void BamRead::SetSequence(std::string s)
  {

    // change the size to accomodate new sequence. Clear the quality string
    //std::cout << "osize " << b->l_data << " calcsize " << (b->core.l_qseq + b->core.l_qname + (b->core.n_cigar<<2)) << std::endl;
    b->data = (uint8_t*)realloc(b->data, b->core.l_qname + s.length() + (b->core.n_cigar<<2));
    
    // copy in the new sequence
    memcpy(b->data + b->core.l_qname + (b->core.n_cigar<<2), (uint8_t*)s.c_str(), s.length());

    }*/

  /*  std::string BamRead::toSam(bam_hdr_t* h) const 
  {
    kstring_t *str;
    sam_format1(h, b, str);
    return std::string(str->s);
    }*/

  std::string BamRead::QualitySequence() const {
    std::string seq = GetZTag("GV");
    if (!seq.length()) 
      seq = Sequence();
    return seq;
  }

  std::ostream& operator<<(std::ostream& out, const BamRead &r)
  {
    //kstring_t *str;
    //sam_format1(hdr, b, str);
    //out << str->s;

    out << bam_get_qname(r.b) << "\t" << r.b->core.flag
	<< "\t" << (r.b->core.tid+1) << "\t" << r.b->core.pos 
	<< "\t" << r.b->core.qual << "\t" << r.CigarString() 
	<< "\t" << (r.b->core.mtid+1) << "\t" << r.b->core.mpos << "\t" 
        << r.b->core.isize 
	<< "\t" << r.Sequence()/* << "\t" << r.Qualities()*/;
    return out;
      
    
  }

  int32_t BamRead::CountSecondaryAlignments() const 
  {
    int xp_count = 0;
    
    // xa tag
    std::string xar_s = GetZTag("XA");
    //r_get_Z_tag(r, "XA", xar_s);
    if (xar_s.length()) {
      xp_count += std::count(xar_s.begin(), xar_s.end(), ';');
    }
    
    // xp tag
    std::string xpr_s = GetZTag("XP");
    //r_get_Z_tag(r, "XP", xpr_s);
    if (xpr_s.length()) {
      xp_count += std::count(xpr_s.begin(), xpr_s.end(), ';');
    }

    return xp_count;
    
  }

  int32_t BamRead::CountNBases() const {
    uint8_t* p = bam_get_seq(b); 
    int32_t n = 0;
    for (int ww = 0; ww < b->core.l_qseq; ww++)
      if (bam_seqi(p,ww) == 15) 
	++n; 
    return n;
  }

  std::string BamRead::QualityTrimmedSequence(int32_t qualTrim, int32_t& startpoint) const {

    int endpoint = -1; //seq.length();
    startpoint = 0;
    int i = 0; 
    
    uint8_t * qual = bam_get_qual(b.get());
    
    // if there is no quality score, return whole thing
    if (qual[0] == 0xff) {
      startpoint = 0;
      return Sequence();
    }

    
    // get the start point (loop forward)
    while(i < b->core.l_qseq) {
      int ps = qual[i];
      if (ps >= qualTrim) {
          startpoint = i;
          break;
	}
	++i;
    }

    // get the end point (loop backwards)
    i = b->core.l_qseq - 1; //seq.length() - 1;
    while(i >= 0) {

      int ps = qual[i];
      
      if (ps >= qualTrim) { //ps >= qualTrim) {
	endpoint = i + 1; // endpoint is one past edge
	break;
      }
      --i;
    }

    // check that they aren't all bad
    if (startpoint == 0 && endpoint == -1) 
      return "";
    
    std::string output = std::string(endpoint-startpoint, 'N');
    try { 
      uint8_t * p = bam_get_seq(b);
      for (int32_t i = startpoint; i < (endpoint - startpoint); ++i) 
	output[i] = BASES[bam_seqi(p,i)];
    } catch (...) {
      std::cerr << "Trying to subset string in BamRead::QualityTrimRead out of bounds. String: " << Sequence() << " start " << startpoint << " length " << (endpoint - startpoint) << std::endl;
    }

    return output;

  }

  void BamRead::AddZTag(std::string tag, std::string val) {
    if (tag.empty() || val.empty())
      return;
    bam_aux_append(b.get(), tag.data(), 'Z', val.length()+1, (uint8_t*)val.c_str());
  }

  std::string BamRead::GetZTag(const std::string& tag) const {
    uint8_t* p = bam_aux_get(b.get(),tag.c_str());
    if (!p)
      return std::string();
    char* pp = bam_aux2Z(p);
    if (!pp) 
      return std::string();
    return std::string(pp);
  }

  
  // get a string tag that might be separted by "x"
  std::vector<std::string> BamRead::GetSmartStringTag(const std::string& tag) const {
    
    std::vector<std::string> out;
    std::string tmp = GetZTag(tag);

    if (tmp.empty())
      return std::vector<std::string>();
    
    if (tmp.find("x") != std::string::npos) {
      std::istringstream iss(tmp);
      std::string line;
      while (std::getline(iss, line, 'x')) {
	out.push_back(line);
      }
    } else {
      out.push_back(tmp);
    }
    
    assert(out.size());
    return out;
    
  }
  
  
  std::vector<int> BamRead::GetSmartIntTag(const std::string& tag) const {
    
    std::vector<int> out;
    std::string tmp;
    
    tmp = GetZTag(tag);
    //r_get_Z_tag(a, tag.c_str(), tmp);
    //assert(tmp.length());
    if (tmp.empty())
      return std::vector<int>();
    
    if (tmp.find("x") != std::string::npos) {
      std::istringstream iss(tmp);
      std::string line;
      while (std::getline(iss, line, 'x')) {
	try { out.push_back(stoi(line)); } catch (...) { std::cerr << "Failed to read parsed int tag " << tag << " for value " << tmp << " with line " << line << std::endl; std::exit(EXIT_FAILURE); }
      }
    } else {
      try { out.push_back(stoi(tmp)); } catch (...) { std::cerr << "Failed to read int tag " << tag << " for value " << tmp << std::endl; std::exit(EXIT_FAILURE); }
    }
    
    assert(out.size());
    return out;
    
  }
  
  
}
