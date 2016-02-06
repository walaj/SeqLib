#include "SnowTools/SnowUtils.h"
#include <random>

/*std::vector<double> SnowTools::getWeightedSum(const std::vector<double>& c) {

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
  

  }*/

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


std::string SnowTools::getFileName(const std::string& s) {

    char sep = '/';
#ifdef _WIN32
    sep = '\\';
#endif
    size_t i = s.rfind(sep, s.length());
    if (i != std::string::npos) {
      return(s.substr(i+1, s.length() - i));
    }
    
    return("");
  }
 


 /*! @function Loops through a text file to count the number of lines.
  * @param file The file to count
  * @param exclude String which, if present in line, causes line to not be counted.
  * @param include String which must be present in the line to be counted.
  * @return Number of valid lines in file
  */
 /*size_t SnowTools::countLines(const std::string &file, const std::string &exclude = "", const std::string &include = "") {
   
   //open the file
   igzstream inFile(file.c_str());
   if (!inFile) 
     return 0;
   
   bool has_include = !include.empty();
   bool has_exclude = !exclude.empty();

   // loop through the file
   size_t count = 0;
   std::string dum;
   while (std::getline(inFile, dum)) {
     if (!include.length() || (dum.find(include) != std::string::npos)) // file must have include
       if (!exclude.length() || (dum.find(exclude) == std::string::npos)) // and not have exclude
	 count++;
     
   }
   
   return count;
   } */
