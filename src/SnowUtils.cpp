#include "SnowTools/SnowUtils.h"
#include <random>

std::vector<double> SnowTools::getWeightedSum(const std::vector<double>& c) {

  double sum = 0;
  for (auto& i : c)
    sum += i;
  
  if (sum == 0)
    return {};

  std::vector<double> wsum;
  double tsum = 0;
  for (auto& i : c) {
    tsum += i;
    wsum.push_back(tsum / sum);
  }

  return wsum;
  

}

int SnowTools::weightedRandom(const std::vector<double>& cs) {

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

void SnowTools::genRandomVals(uint32_t &i1, uint32_t &i2, const uint32_t &max, uint32_t seed) {

  std::random_device rd; // obtain a random number from hardware
  if (seed == 0) {
    seed = rd();
  }

  std::mt19937 eng(seed); // seed the generator
  std::uniform_int_distribution<uint32_t> distr(0, max); // define the range
  i1 = distr(eng);
  i2 = distr(eng);

}

void SnowTools::genRandomValue(uint32_t &i, const uint32_t &max, uint32_t seed) {
  
  std::random_device rd; // obtain a random number from hardware
  if (seed == 0) {
    seed = rd();
  }
  
  std::mt19937 eng(seed); // seed the generator
  std::uniform_int_distribution<uint32_t> distr(0, max); // define the range
  i = distr(eng);
}
 
