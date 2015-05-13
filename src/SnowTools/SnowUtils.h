#ifndef SNOWUTILS_H
#define SNOWUTILS_H

#include <string>
#include <time.h>
#include <ctime>
#include <vector>
#include <unistd.h>
#include <sstream>
#include <cmath>
#include <algorithm>

#include "SnowTools/gzstream.h"
#include "SnowTools/GenomicRegion.h"

namespace SnowTools {

  /** Check if a file is readable and exists
   * @param name Name of a file to test
   * @return File is readable and exists
   */
  inline bool read_access_test (const std::string& name) {
    return (access (name.c_str(), R_OK) == 0); 
  }

  /** Format a number to include commas
   * @param data Number to format
   * @param Formatted number
   */
  template <typename T> 
    std::string AddCommas(T data) { 
    std::stringstream ss; 
    ss << data; 
    std::string s = ss.str();
    if (s.length() > 3)
      for (int i = s.length()-3; i > 0; i -= 3)
	s.insert(i,",");
    return s;
  }
  
  inline std::string displayRuntime(const timespec start) {
    
    struct timespec finish;
    clock_gettime(CLOCK_MONOTONIC, &finish);
    double elapsed = (finish.tv_sec - start.tv_sec);
    int t = clock()/CLOCKS_PER_SEC;
    int min = (int)std::floor(elapsed / 60.0);
    int sec = (int)(elapsed-min*60);
    char buffer[100];
    sprintf (buffer, "CPU: %4dm%02ds Wall: %4dm%02ds", 
	     (int)floor( ((double)t) /60.0), t % 60, min, sec);
    buffer[99] = '\0';
    return std::string(buffer);
  }

  /** Deprecated
   */
  inline void rcomplement(std::string &a) {
    
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
  
  // calculate the percentage
  template <typename T> inline int percentCalc(T numer, T denom) {
    if (denom <= 0)
      return 0;
    int perc  = static_cast<int>(floor((float)numer / (float)denom * 100.0));
    return perc;
  }

  // remove the last character from a string
  inline std::string cutLastChar(std::string in) {
    if (in.length() == 0)
      return in;
    else 
      return in.substr(0, in.length() - 1);
  }
  
  // remove substrings from a string
 inline std::string scrubString(std::string toscrub, std::string toremove) {
   std::string::size_type i = toscrub.find(toremove);
   while (i != std::string::npos) {
     toscrub.erase(i, toremove.length());
     i = toscrub.find(toremove);
   }
   return toscrub;
 }

 void genRandomVals(uint32_t &i1, uint32_t &i2, const uint32_t &max);

 void genRandomValue(uint32_t &i, const uint32_t &max);

 /*! @function Loops through a text file to count the number of lines.
  * @param file The file to count
  * @param exclude String which, if present in line, causes line to not be counted.
  * @param include String which must be present in the line to be counted.
  * @return Number of valid lines in file
  */
 inline size_t countLines(const std::string &file, const std::string &exclude = "", const std::string &include = "") {
   
   //open the file
   igzstream inFile(file.c_str());
   if (!inFile) 
     return 0;
   
   // loop through the file
   size_t count = 0;
   std::string dum;
   while (std::getline(inFile, dum)) {
     if (!include.length() || (dum.find(include) != std::string::npos)) // file must have include
       if (!exclude.length() || (dum.find(exclude) == std::string::npos)) // and not have exclude
	 count++;
     
   }
   
   return count;
 } 

}

#endif
