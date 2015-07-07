#ifndef SNOWTOOLS_HISTOGRAM_H__
#define SNOWTOOLS_HISTOGRAM_H__

#include <iostream>
#include <cassert>
#include <string>
#include <utility>
#include <vector>
#include <fstream>

#include "IntervalTree.h"

class Bin;

typedef Interval<Bin> BinInterval;
typedef IntervalTree<Bin> BinIntervalTree;
typedef std::vector<BinInterval> BinIntervalVector;

typedef uint32_t S;

#define INTERCHR 250000000

namespace SnowTools {

  class Matrix;
  class Histogram;

/** Stores one bin in a Histogram
 */
class Bin {

  friend class Histogram;

 public:
  
  /** Construct a new object with 0 count and range [0,1]
   */
  Bin() : m_count(0) 
  { 
    bounds = {0,1}; 
  }
    
    /** Output the bin in format "start range, end range, count"
     */
    friend std::ostream& operator<<(std::ostream &out, const Bin &b) {
      out << b.bounds.first << "," << b.bounds.second << "," << b.m_count;
      return out;
    }

    /** Return the number of counts in this histogram bin 
     */
    size_t getCount() const { return m_count; }
    
    /** Check if a value fits within the range of this bin 
     * @param dist Distance value to check if its in this range
     * @return true if the value is within the range
     */

    /** Check if this bin contains a value
     * @param elem Value to check if its in this range
     * @return true if the value is within the range
     */
    bool contains(const int32_t& elem) const; 

    /** Define bin comparison operator by location of left bound, then right */
    bool operator < (const Bin& b) const;
    
    /** Decrement the histogram bin by one. 
     * Note that this is the prefix version only
     */
    Bin& operator--();

    /** Increment the histogram bin by one. 
     * Note that this is the prefix version only
     */
    Bin& operator++();

 private:
    size_t m_count;
    std::pair<int32_t, int32_t> bounds; //@! was"bin";
};

/** Class to store histogram of numeric values.
 *
 * The bins of the Histogram are not uniformly spaced, and their ranges determined 
 * by partitioning the spans it tablulates into uniform quantiles when initialized
 * by Histogram::initialSpans(). As elements are added and removed this initial bin
 * definition remains constant.
 */
class Histogram { 

 private:

  std::vector<int32_t> m_ind;

 public:

  std::vector<Bin> m_bins;
  /** Construct an empty histogram
   */
  Histogram() {}

  /** Construct a new histogram with bins spaced evenly
   */
  Histogram(const int32_t& start, const int32_t& end, const uint32_t& width);

  std::string toFileString() const;

  friend std::ostream& operator<<(std::ostream &out, const Histogram &h) {
    for (auto& i : h.m_bins)
      out << i << std::endl;
    return out;
  }

  /** Return iterator to the fist bin
   */
  std::vector<int32_t>::iterator begin() { return m_ind.begin(); }
  
  /** Return iterator to the last bin
   */
  std::vector<int32_t>::iterator end() { return m_ind.end(); }
  
  /** initialize histogram from a vector of spans 
   */
  void initialSpans(size_t num_bins, std::vector<S>* pspanv, size_t min_bin_width = 0);

  /** Add an element to the histogram
   * @param elem Length of event to add
   */
  void addElem(const int32_t &elem);

  /** Remove a span from the histogram
   * @param span Length of event to remove
   */
  void removeElem(const int32_t &elem);

  /** Output to CSV file like: bin_start,bin_end,count
   */
  void toCSV(std::ofstream &fs);

  /** Return the total number of elements in the Histogram
   */
  size_t totalCount() const {
    size_t tot = 0;
    for (auto&  i : m_bins)
      tot += i.getCount();
    return tot;
  }
  
  /** Get count for a histogram bin
   * @param i Bin index
   * @return number of events in histogram bin
   */
  int32_t binCount(size_t i) { return m_bins[i].getCount(); }

  /** Get number of bins in histogram
   * @return Number of bins in histogram
   */
  size_t numBins() { return m_bins.size(); }

  /** Find bin corresponding to a span
   * @param elem Event length
   * @return Bin containing event length
   */
  size_t retrieveBinID(const int32_t& elem) const;

};

}
#endif
