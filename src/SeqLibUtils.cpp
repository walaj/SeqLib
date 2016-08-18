#include "SeqLib/SeqLibUtils.h"
#include <random>


std::string SeqLib::scrubString(const std::string& toscrub, const std::string& toremove) 
{
  if (toscrub.empty() || toremove.empty())
    return toscrub;

  std::string::size_type i = toscrub.find(toremove);
  if (i == std::string::npos)
    return toscrub;
  
  std::string ts = toscrub;
  while (i != std::string::npos) {
    ts.erase(i, toremove.length());
    i = ts.find(toremove);
  }
  return ts;
}

int SeqLib::weightedRandom(const std::vector<double>& cs) {

  // get a weighted random number
  size_t al = 0;
  double rand_val = rand() % 1000;
  while (al < cs.size()) {
    if (rand_val <= cs[al] * 1000) 
      return al;
    ++al;
  }
  return al;
}

