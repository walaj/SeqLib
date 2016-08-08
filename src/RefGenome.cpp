#include "SnowTools/RefGenome.h"

#include "htslib/khash.h"
#include <stdexcept>
#include "SnowTools/SnowUtils.h"

namespace SnowTools {

  void RefGenome::retrieveIndex(const std::string& file) {

    // clear the old one
    if (index)  
      fai_destroy(index);
    
    index = nullptr;
    
    // check that its readable
    if (!read_access_test(file)) {
      throw std::invalid_argument("RefGenome: file not found - " + file);
    }
    
    // load it in
    index = fai_load(file.c_str());

  }
  
  RefGenome::RefGenome(const std::string& file ) {

    retrieveIndex(file);
    
  }

  std::string RefGenome::queryRegion(const std::string& chr_name, int32_t p1, int32_t p2) const {
    
    // check that we ahve a loaded index
    if (!index) 
      throw std::invalid_argument("RefGenome::queryRegion index not loaded");

    // check input is OK
    if (p1 > p2)
      throw std::invalid_argument("RefGenome::queryRegion p1 must be <= p2");
    if (p1 < 0)
      throw std::invalid_argument("RefGenome::queryRegion p1 must be >= 0");
    if (p2 < 0)
      throw std::invalid_argument("RefGenome::queryRegion p2 must be >= 0");

    // set the iterator to the chr position
        /*    

    khiter_t iter;
    iter = kh_get_s(index->hash, const_cast<char*>(chr_name.c_str()));
	  if (iter == kh_end(fai->hash))
      throw std::out_of_range("RefGenome::queryRegion - Could not find chr " + chr_name + " in faidx hash table");
    faidx1_t val;
    val = kh_value(index->hash, iter);

    if (val.len <= p1)
      throw std::out_of_range("RefGenome::queryRegion - p1 val of " + std::to_string(p1) + " is greater than chr len of " + std::to_string(val.len) +
			      " for chr " << chr_name); 
    */
    int len; 
    // why is this not thread safe???
    char * f = faidx_fetch_seq(index, const_cast<char*>(chr_name.c_str()), p1, p2, &len);

    if (!f)
      throw std::invalid_argument("RefGenome::queryRegion - Could not find valid sequence");
    
    std::string out(f);
    
    free(f); 

    if (out.empty())
      throw std::invalid_argument("RefGenome::queryRegion - Returning empty query on " + chr_name + ":" + std::to_string(p1) + "-" + std::to_string(p2));
    
    return (out);

  }

}
