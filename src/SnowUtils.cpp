#include "SnowTools/SnowUtils.h"
#include <random>


void SnowTools::genRandomVals(uint32_t &i1, uint32_t &i2, const uint32_t &max) {
  
  std::random_device rd; // obtain a random number from hardware
  std::mt19937 eng(rd()); // seed the generator
  std::uniform_int_distribution<uint32_t> distr(0, max); // define the range
  i1 = distr(eng);
  i2 = distr(eng);

}

void SnowTools::genRandomValue(uint32_t &i, const uint32_t &max) {
  
  std::random_device rd; // obtain a random number from hardware
  std::mt19937 eng(rd()); // seed the generator
  std::uniform_int_distribution<uint32_t> distr(0, max); // define the range
  i = distr(eng);
}
