#ifndef SWAP_GENOMIC_REGION_COLLECTION_H__
#define SWAP_GENOMIC_REGION_COLLECTION_H__

#include <vector>
#include <string>
#include <cstdlib>
#include <list>
#include <unordered_map>
#include <memory>

//#ifdef BOOST_VERSION
//#include "boost/icl/interval_set.hpp"
//#endif

#include "SeqLib/IntervalTree.h"
#include "SeqLib/GenomicRegionCollection.h"
#include "SeqLib/BamRecord.h"

namespace SeqLib {

/** Class to store vector of intervals on the genome */
typedef TInterval<int32_t> GenomicInterval;
typedef std::unordered_map<int, std::vector<GenomicInterval> > GenomicIntervalMap;
typedef TIntervalTree<int32_t> GenomicIntervalTree;
typedef std::unordered_map<int, GenomicIntervalTree> GenomicIntervalTreeMap;
typedef std::vector<GenomicInterval> GenomicIntervalVector;

  /** @brief Template class for a collection of children of GenomicRegion objects.
   *
   * Contains a set of GenomicRegion objects (or their children), and allows
   * interval queries, reading/writing BED files, etc.
   */
template<typename T=GenomicRegion>
class GenomicRegionCollection {

 public:

  /** Construct an empty GenomicRegionCollection 
   */
 GenomicRegionCollection();

 ~GenomicRegionCollection();
 
  /** Construct from a plain vector of GenomicRegion objects
   */
  GenomicRegionCollection(std::vector<T>& vec);

  /** Construct from a single GenomicRegion
   */
  GenomicRegionCollection(const T& gr);

  /** Construct from a vector of reads
   */
  GenomicRegionCollection(const BamRecordVector& brv);

  /** Construct a GenomicRegionCollection with overlapping intervals
   * 
   * @param width Desired bin width
   * @param ovlp Amount that the bins should overlap
   * @param gr GenomicRegion to divide into smaller overlapping bins
   */
  GenomicRegionCollection(int width, int ovlp, const T &gr);

  /** Read in a MuTect call-stats file and adds to GenomicRegionCollection object.
   *
   * Reads a MuTect call-stats file and imports only
   * lines with KEEP marked. 
   * @param file Path to call-stats file
   * @param pad Amount to pad intervals by
   */
 void readMuTect(const std::string &file, int pad, const SeqLib::BamHeader& hdr);

  /** Read in a BED file and adds to GenomicRegionCollection object
   * @param file Path to BED file
   * @param pad Amount to pad intervals by
   */
 void readBEDfile(const std::string &file, int pad, const SeqLib::BamHeader& hdr);

  /** Read in a VCF file and adds to GenomicRegionCollection object
   * @param file Path to BED file
   * @param pad Amount to pad intervals by
   */
 void readVCFfile(const std::string &file, int pad, const SeqLib::BamHeader& hdr);

  /** Read in a text file (can be gzipped) and construct a GenomicRegionCollection
   *
   * This function will automatically detect which file type is being input:
   * -- ends in .vcf -> readVCFfile
   * -- ends in .bed -> readBEDfile
   * -- header contains "MuTect" -> readMuTect
   * The values are appended to existing vector of GenomicRegion objects
   * @param file Text file to read and store intervals
   * @param hdr BamHeader to serve as dictionary for chromosomes
   */
   GenomicRegionCollection(const std::string &file, const BamHeader& hdr);

  /** Fill in the GenomicIntervalTreeMap stored in this object. 
   *
   * A GenomicIntervalTreeMap is an unordered_map of GenomicIntervalTrees for 
   * each chromosome. A GenomicIntervalTree is an interval tree on the ranges
   * defined by the genomic interval, with cargo set at the same GenomicRegion object.
   */
  void createTreeMap();
  
  /** Send the GenomicRegionCollection to a BED file
  * @param file Output file
  */
  void sendToBED(const std::string file);
  
  /** Reduces the GenomicRegion objects to minimal set
   */
  void mergeOverlappingIntervals();

  /** Return the number of GenomicRegions stored 
   */
  size_t size() const { return m_grv->size(); }

  /** Add a new GenomicRegion to end
   */
 void add(const T& g) { m_grv->push_back(g); /*createTreeMap();*/ }

  /** Is this object empty?
   */
  bool empty() const { return !m_grv->size(); }

  /** Clear out all of the GenomicRegion objects
   */
  void clear() { m_grv->clear(); 
		 m_tree->clear(); 
		 idx = 0;
  }

  /** Retrieve a GenomicRegion at given index. 
   * 
   * Note that this does not move the idx iterator, which is 
   * used to loop through all the regions. Throws an exception
   * if the index is out of bounds.
   * @return GenomicRegion pointed to by index i
   */
  const T& at(size_t i) const;

  /** Find overlaps between this vector and input GenomicRegion.
   *
   * Requires that the GenomicIntervalTreeMap have been created first
   * @param gr Region to test
   * @return Number of overlapping elements in this GenomicRegionCollection
   */
 size_t findOverlapping(const T &gr) const;

 bool overlapSameBin(const T &gr1, const T &gr2) const;

 size_t countContained(const T &gr);
 
 template<class K>
 GenomicRegionCollection<GenomicRegion> findOverlaps(GenomicRegionCollection<K> &subject, std::vector<int32_t>& query_id, std::vector<int32_t>& subject_id, bool ignore_strand = false) const;

 template<class K>
 GenomicRegionCollection<GenomicRegion> findOverlaps(const K& gr, bool ignore_strand) const;

 /** The total amount spanned by this collection
  */
 int width() const { int wid = 0; for (auto& i : m_grv) wid += i.width(); return wid; }

 /** Increase the left and right ends of each contained GenomicRegion by 
  * the pad value.
  * @param v Amount to pad each end by. Result is increase in width by 2*pad.
  */
 void pad(int v);

 /** Set the i'th GenomicRegion */
 T& operator[](size_t i) { return m_grv->at(i); }
 
 /** Retreive the i'th GenomicRegion */
 const T& operator[](size_t i) const { return m_grv->at(i); }
 
  /** Add two GenomicRegionCollection objects together
   */
  void concat(const GenomicRegionCollection<T>& g);

  /** Output the GenomicRegionCollection to a BED format
   *
   * Currently just outputs chr   pos   pos2 with no header.
   * @return BED formated string reprsentation 
   */
  std::string sendToBED() const;

  /** Fill the next GenomicRegion object. 
   * @return false if add end of vector
   */
  bool getNextGenomicRegion(T& gr);

  void gsort();

 /** Expand all the elements so they are sorted and become adjacent by 
     stretching them to the right up to max */
  void SortAndStretchRight(int max);

 /** Expand all the elements so they are sorted and become adjacent by 
     stretching them to the left down to min */
 void SortAndStretchLeft(int min);

  /** Rewind the element pointer to the first GenomicRegion
   */
  void rewind() { idx = 0; }

 GenomicRegionVector asGenomicRegionVector() const;
 
  typename std::vector<T>::iterator begin() { return m_grv->begin(); } 

  typename std::vector<T>::iterator end() { return m_grv->end(); } 
 
  typename std::vector<T>::const_iterator begin() const { return m_grv->begin(); } 

  typename std::vector<T>::const_iterator end() const { return m_grv->end(); } 

  GenomicRegionCollection<GenomicRegion> intersection(GenomicRegionCollection<GenomicRegion>& subject, bool ignore_strand = false);
 
 //#ifdef BOOST_VERSION
 // boost::icl::interval_set<int> m_set;
 //#endif
 
 //#ifdef BOOST_VERSION
 //GenomicRegionCollection<GenomicRegion> complement(GenomicRegionCollection<GenomicRegion>& subject, bool ignore_strand = false);
 //#endif

 private:
 
 // always construct this object any time m_grv is modifed
 std::shared_ptr<GenomicIntervalTreeMap> m_tree;
 
 // hold the genomic regions
 std::shared_ptr<std::vector<T>> m_grv; 
 
 // index for current GenomicRegion
 size_t idx = 0;

 // open the memory
 void __allocate_grc();

};

typedef GenomicRegionCollection<GenomicRegion> GRC;

}

#include "GenomicRegionCollection.cpp"

#endif
