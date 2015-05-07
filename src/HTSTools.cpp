#include "SnowTools/HTSTools.h"
#include <cassert>
#include <sstream>
#include <iostream>

// Trim the sequence by removing low quality bases from either end
int32_t SnowTools::qualityTrimRead(int qualTrim, int32_t &startpoint, Read &r) {

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



// get a string tag that might be separted by "x"
std::vector<std::string> SnowTools::GetStringTag(const Read& a, const std::string tag) {
  
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
void SnowTools::SmartAddTag(Read &a, const std::string tag, const std::string val) {
  
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
std::vector<int> SnowTools::GetIntTag(const Read& a, const std::string tag) {
  
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

void SnowTools::rcomplement(std::string &a) {

  std::reverse(&a[0], &a[a.size()]);
  std::string::iterator it = a.begin();
  for (; it != a.end(); it++)
    if (*it == 'A')
      *it = 'T';
    else if (*it == 'T')
      *it = 'A';
    else if (*it == 'C')
      *it = 'G';
    else
      *it = 'C';
}
