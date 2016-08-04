#include "SnowTools/RefGenome.h"

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

  std::string RefGenome::queryRegion(const std::string& chr_name, int32_t p1, int32_t p2) {
    
    // check that we ahve a loaded index
    if (!index) {
      throw std::invalid_argument("RefGenome::queryRegion index not loaded");
    }

    int len;
    char * f = faidx_fetch_seq(index, const_cast<char*>(chr_name.c_str()), p1, p2, &len);

    if (!f)
      throw std::invalid_argument("RefGenome::queryRegion - Could not find valid sequence");

    std::string out(f);

    free(f);

    return (out);

  }

}
