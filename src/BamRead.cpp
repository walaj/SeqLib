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

  void BamRead::SmartAddTag(const std::string& tag, const std::string& val)
  {
    // get the old tag
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
    AddZTag(tag, tmp);
    
  }

  
  void BamRead::SetQname(std::string n)
  {
    // copy out the non-qname data
    size_t nonq_len = b->l_data - b->core.l_qname;
    uint8_t* nonq = (uint8_t*)malloc(nonq_len);
    memcpy(nonq, b->data + b->core.l_qname, nonq_len);

    // clear the old data and alloc the new amount 
    free(b->data);
    b->data = (uint8_t*)calloc(nonq_len + n.length() + 1, 1);
    
    // add in the new qname
    memcpy(b->data, (uint8_t*)n.c_str(), n.length() + 1); // +1 for \0

    // update the sizes
    b->l_data = b->l_data - b->core.l_qname + n.length() + 1;
    b->core.l_qname = n.length() + 1;    
    
    // copy over the old data
    memcpy(b->data + b->core.l_qname, nonq, nonq_len);
    free(nonq);
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

  std::ostream& operator<<(std::ostream& out, const BamRead &r)
  {
    //kstring_t *str;
    //sam_format1(hdr, b, str);
    //out << str->s;

    out << bam_get_qname(r.b) << "\t" << r.b->core.flag
	<< "\t" << r.b->core.tid << "\t" << r.b->core.pos 
	<< "\t" << r.b->core.qual << "\t" << r.CigarString() 
	<< "\t" << "=" << "\t" << 0 
	<< "\t" << r.Sequence() << "\t" << "*";
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
    
    std::string output = "";
    try { 
      output = Sequence().substr(startpoint, endpoint - startpoint);
    } catch (...) {
      std::cerr << "Trying to subset string in BamRead::QualityTrimRead out of bounds. String: " << Sequence() << " start " << startpoint << " length " << (endpoint - startpoint) << std::endl;
    }

    return output;

  }
  
  // get a string tag that might be separted by "x"
  std::vector<std::string> BamRead::GetSmartStringTag(const std::string& tag) const {
    
    std::vector<std::string> out;
    std::string tmp = GetZTag(tag);
    assert(tmp.length());
    
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
    assert(tmp.length());
    
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
