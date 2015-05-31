#include "SnowTools/HTSTools.h"
#include <cassert>
#include <sstream>
#include <iostream>

namespace SnowTools {

// Trim the sequence by removing low quality bases from either end
int32_t qualityTrimRead(int qualTrim, int32_t &startpoint, Read &r) {

    int endpoint = -1; //seq.length();
    startpoint = 0;
    int i = 0; 

    uint8_t * qual = bam_get_qual(r.get());

    // get the start point (loop forward)
    while(i < r->core.l_qseq) {
      int ps = qual[i];
      if (ps >= qualTrim) {
          startpoint = i;
          break;
	}
	i++;
    }

    // get the end point (loop backwards)
    i = r->core.l_qseq - 1; //seq.length() - 1;
    while(i >= 0) {

      int ps = qual[i];
      
      if (ps >= qualTrim) { //ps >= qualTrim) {
	endpoint = i + 1; // endpoint is one past edge
	break;
      }
      i--;
    }
    
    // check that they aren't all bad
    if (startpoint == 0 && endpoint == -1) 
      return 0;


    return (endpoint - startpoint);

}

void removeAllTags(Read& a) 
{
  size_t keep = (a->core.n_cigar<<2) + a->core.l_qname + ((a->core.l_qseq + 1)>>1) + a->core.l_qseq;
  a->data = (uint8_t*)realloc(a->data, keep); // free the end, which has aux data
  a->l_data = keep;
}

// get a string tag that might be separted by "x"
std::vector<std::string> GetStringTag(const Read& a, const std::string tag) {
  
  std::vector<std::string> out;
  std::string tmp;
  
  r_get_Z_tag(a, tag.c_str(), tmp);
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

// add a tag that might already be there, separete by 'x'
void SmartAddTag(Read &a, const std::string tag, const std::string val) {
  
  std::string tmp;
  r_get_Z_tag(a, tag.c_str(), tmp);

  if (tmp.length()) {
    tmp += "x"  + val;
    r_remove_tag(a, tag.c_str());
    r_add_Z_tag(a, tag.c_str(), tmp);
  } else { // normal with no x
    r_add_Z_tag(a, tag.c_str(), val);
  }

}

// get an integer tag that might be separted by "x"
std::vector<int> GetIntTag(const Read& a, const std::string tag) {
  
  std::vector<int> out;
  std::string tmp;
  
  r_get_Z_tag(a, tag.c_str(), tmp);
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

// return a sequence from the reference
std::string getRefSequence(const GenomicRegion &gr, faidx_t * fi) {

  int len;
  std::string chrstring = GenomicRegion::chrToString(gr.chr);
  char * seq = faidx_fetch_seq(fi, const_cast<char*>(chrstring.c_str()), gr.pos1-1, gr.pos2-1, &len);
  
  if (seq) {
    return std::string(seq);
  } else {
    //cout << "Failed to get reference sequence at " << gr << endl;
    return "LOAD_FAIL";
  }

}


}

