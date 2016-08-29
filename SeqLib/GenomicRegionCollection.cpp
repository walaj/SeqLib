#include "SeqLib/GenomicRegionCollection.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <cassert>
#include <set>
#include <unordered_set>

#ifdef HAVE_BOOST_ICL_INTERVAL_SET_HPP  
#include "boost/icl/interval_set.hpp"
#endif

//#define DEBUG_OVERLAPS 1

namespace SeqLib {

  static bool header_has_chr_string = false;

  template<class T>
  void GenomicRegionCollection<T>::CoordinateSort() {
    
    std::sort(m_grv->begin(), m_grv->end());
  }

  template<class T>
  void GenomicRegionCollection<T>::SortAndStretchRight(int max) {

    if (!m_grv->size())
      return;
    
    std::sort(m_grv->begin(), m_grv->end());

    if (max > 0 && max < m_grv->back().pos2) {
      //std::cerr << "GenomicRegionCollection::SortAndStrech Can't stretch to max, as we are already past max. Max: " << max << " highest GRC " << m_grv->back() << std::endl;
      throw std::out_of_range("GenomicRegionCollection::SortAndStrech Can't stretch to max, as we are already past max.");
    }

    for (size_t i = 0; i < m_grv->size() - 1; ++i) {
      m_grv->at(i).pos2 = m_grv->at(i+1).pos1 - 1;
    }

    if (max > 0)
      m_grv->back().pos2 = max;
    
  }

  template<class T>
  void GenomicRegionCollection<T>::SortAndStretchLeft(int min) {

    if (!m_grv->size())
      return;

    std::sort(m_grv->begin(), m_grv->end());

    if (min >= 0 && min < m_grv->begin()->pos1) {
      std::cerr << "GenomicRegionCollection::SortAndStrechLeft Can't stretch to min, as we are already past max" << std::endl;
      exit(EXIT_FAILURE);
    }

    if (min >= 0)
      m_grv->at(0).pos1 = min;

    for (size_t i = 1; i < m_grv->size(); ++i) {
      m_grv->at(i).pos1 = m_grv->at(i-1).pos2 + 1;
    }
    
  }



  /*template<class T>
bool GenomicRegionCollection<T>::ReadMuTect(const std::string &file, const BamHeader& hdr) {

  std::string curr_chr = "dum";
  
  std::ifstream iss(file.c_str());
  if (!iss || file.length() == 0) { 
    std::cerr << "MuTect call-stats file does not exist: " << file << std::endl;
    return false;
  }

  std::string line;
  while (std::getline(iss, line, '\n')) {
    size_t counter = 0;
      std::string chr, pos, judge;
      std::istringstream iss_line(line);
      std::string val;
      if (line.find("KEEP") != std::string::npos) {
	while(std::getline(iss_line, val, '\t')) {
	  switch (counter) { 
	  case 0 : chr = val; break; 
	  case 1 : pos = val; break;
	  }
	  if (counter >= 1)
	    break;
	  ++counter;
	  
	  if (curr_chr != chr) {
	    std::cerr << "...reading MuTect call-stats -- chr" << chr << std::endl;
	    curr_chr = chr;
	  }

	}

	// parse the strings and send to genomci region
	T gr(chr, pos, pos, hdr);
	if (gr.chr >= 0) {
	  m_grv->push_back(gr);
	}

      } // end "keep" conditional
    } // end main while

  return true;
}
  */

template<class T>
bool GenomicRegionCollection<T>::ReadBED(const std::string & file, const BamHeader& hdr) {

  std::ifstream iss(file.c_str());
  if (!iss || file.length() == 0) { 
    std::cerr << "BED file does not exist: " << file << std::endl;
    return false;
  }

  std::string line;
  std::string curr_chr = "-1";
  while (std::getline(iss, line, '\n')) {

    size_t counter = 0;
    std::string chr, pos1, pos2;
    std::istringstream iss_line(line);
    std::string val;
    
    if (line.find("#") == std::string::npos) {
      while(std::getline(iss_line, val, '\t')) {
	switch (counter) { 
	case 0 : chr = val; break; 
	case 1 : pos1 = val; break;
	case 2 : pos2 = val; break;
	}
	if (counter >= 2)
	  break;
	++counter;
	
	if (chr != curr_chr) {
	  //std::cerr << "...reading from BED - chr" << chr << std::endl;
	  curr_chr = chr;
	}
	
      }

      if (header_has_chr_string) {
	//if (chr == "X" || chr == "Y" || std::stoi(chr) < 23)
	  chr = "chr" + chr;
      }
      
      // construct the GenomicRegion
      T gr(chr, pos1, pos2, hdr);

      if (gr.chr >= 0) {
	m_grv->push_back(gr);
      }
	
	//}
    } // end "keep" conditional
  } // end main while

  return true;
}

template<class T>
bool GenomicRegionCollection<T>::ReadVCF(const std::string & file, const BamHeader& hdr) {

  std::ifstream iss(file.c_str());
  if (!iss || file.length() == 0) { 
    std::cerr << "VCF file does not exist: " << file << std::endl;
    return false;
  }

  std::string line;
  
  while (std::getline(iss, line, '\n')) {
    if (line.length() > 0) {
      if (line.at(0) != '#') { // its a valid line
	std::istringstream iss_this(line);
	int count = 0;
	std::string val, chr, pos;
	
	while (std::getline(iss_this, val, '\t')) {
	  switch (count) {
	  case 0 : chr = val;
	  case 1 : pos = val;
	  }
	  count++;
	}
	if (count < 3) {
	  std::cerr << "Didn't parse VCF line properly: " << line << std::endl;
	  return false;
	} else {

	  // construct the GenomicRegion
	  T gr(chr, pos, pos, hdr);
	  if (gr.chr >= 0) {
	    m_grv->push_back(gr);
	  }
	}
	
      }
    }
  }
  
  return true;
  
}

template<class T>
GenomicRegionCollection<T>::GenomicRegionCollection(const std::string &file, const BamHeader& hdr) {

  std::ifstream iss(file.c_str());
  if (!iss || file.length() == 0) { 
    std::cerr << "Region file does not exist: " << file << std::endl;
    return;
  }

  __allocate_grc();

  // get the header line to check format
  std::string header;
  if (!std::getline(iss, header, '\n')) {
    std::cerr << "Region file is empty: " << file << std::endl;
    return;
  }
  iss.close();

  // MUTECT CALL STATS
  if ((header.find("MuTect") != std::string::npos) ||
      (file.find("call_stats") != std::string::npos) || 
      (file.find("callstats") != std::string::npos))
    std::cerr << "MuTect reading not currently available" << std::endl; //ReadMuTect(file, hdr);
  // BED file
  else if (file.find(".bed") != std::string::npos)
    ReadBED(file, hdr);
  // VCF file
  else if (file.find(".vcf") != std::string::npos) 
    ReadVCF(file, hdr);
  else // default is BED file
    ReadBED(file, hdr);    
}

// reduce a set of GenomicRegions into the minium overlapping set (same as GenomicRanges "reduce")
template <class T>
void GenomicRegionCollection<T>::MergeOverlappingIntervals() {

  // make the list
  std::list<T> intervals(m_grv->begin(), m_grv->end());

  intervals.sort();
  typename std::list<T>::iterator inext(intervals.begin());
  ++inext;
  for (typename std::list<T>::iterator i(intervals.begin()), iend(intervals.end()); inext != iend;) {
    if((i->pos2 > inext->pos1) && (i->chr == inext->chr))
      {
	if(i->pos2 >= inext->pos2) intervals.erase(inext++);
	else if(i->pos2 < inext->pos2)
	  { i->pos2 = inext->pos2; intervals.erase(inext++); }
      }
    else { ++i; ++inext; }
  }

  // move it over to a grv
  m_grv->clear(); // clear the old data 
  std::vector<T> v{ std::make_move_iterator(std::begin(intervals)), 
      std::make_move_iterator(std::end(intervals)) };
  m_grv->insert(m_grv->end(), v.begin(), v.end());
  //m_grv = v;

  // recreate the interval tree
  CreateTreeMap();

}

template <class T>
GenomicRegionVector GenomicRegionCollection<T>::AsGenomicRegionVector() const { 
  GenomicRegionVector gg;
  for (auto& i : *m_grv)
    gg.push_back(GenomicRegion(i.chr, i.pos1, i.pos2, i.strand));
  return gg; 
} 
  
template <class T>
void GenomicRegionCollection<T>::CreateTreeMap() {

  // sort the genomic intervals
  sort(m_grv->begin(), m_grv->end());

  // loop through and make the intervals for each chromosome
  GenomicIntervalMap map;
  for (size_t i = 0; i < m_grv->size(); ++i) {
    map[m_grv->at(i).chr].push_back(GenomicInterval(m_grv->at(i).pos1, m_grv->at(i).pos2, i));
  }

  // for each chr, make the tree from the intervals
  for (auto it : map) {
    GenomicIntervalTreeMap::iterator ff = m_tree->find(it.first);
    if (ff != m_tree->end())
      ff->second = GenomicIntervalTree(it.second);
    else
      m_tree->insert(std::pair<int, GenomicIntervalTree>(it.first, GenomicIntervalTree(it.second)));
    //old //m_tree[it.first] = GenomicIntervalTree(it.second);
  }

}

template<class T>
int GenomicRegionCollection<T>::TotalWidth() const{ 
  int wid = 0; 
  for (auto& i : *m_grv) 
    wid += i.Width(); 
  return wid; 
}

// divide a region into pieces of width and overlaps
template<class T>
GenomicRegionCollection<T>::GenomicRegionCollection(int width, int ovlp, const T &gr) {

  __allocate_grc();

  // undefined otherwise
  assert(width > ovlp);
  if (width >= gr.Width()) {
    m_grv->push_back(gr);
    return;
  }

  int32_t start = gr.pos1;
  int32_t end = gr.pos1 + width;

  // region is smaller than width
  if ( end > gr.pos2 ) {
    std::cerr << "GenomicRegionCollection constructor: GenomicRegion is smaller than bin width" << std::endl;
    return; 
  }

  // loop through the sizes until done
  while (end <= gr.pos2) {
    m_grv->push_back(T(gr.chr, start, end));
    end += width - ovlp; // make the new one
    start += width - ovlp;
  }
  assert(m_grv->size() > 0);
  
  // finish the last one if we need to
  if (m_grv->back().pos2 != gr.pos2) {
    start = m_grv->back().pos2 - ovlp; //width;
    end = gr.pos2;
    m_grv->push_back(T(gr.chr, start, end));
  }

}

template<class T>
size_t GenomicRegionCollection<T>::CountOverlaps(const T &gr) const {

  if (m_tree->size() == 0 && m_grv->size() != 0) 
    {
      std::cerr << "!!!!!! WARNING: Trying to find overlaps on empty tree. Need to run this->createTreeMap() somewhere " << std::endl;
      return 0;
    }

  GenomicIntervalVector giv;

  GenomicIntervalTreeMap::const_iterator ff = m_tree->find(gr.chr);
  if (ff == m_tree->end())
    return 0;
  ff->second.findOverlapping(gr.pos1, gr.pos2, giv);
  return (giv.size());
}

  template<class T>
  bool GenomicRegionCollection<T>::OverlapSameInterval(const T &gr1, const T &gr2) const {

    if (m_tree->size() == 0 && m_grv->size() != 0) 
     {
     std::cerr << "!!!!!! WARNING: Trying to find overlaps on empty tree. Need to run this->createTreeMap() somewhere " << std::endl;
     return 0;
    }
  
  // events on diff chr do not overlap same bin
  if (gr1.chr != gr2.chr)
    return false;

  GenomicIntervalVector giv1, giv2;
  GenomicIntervalTreeMap::const_iterator ff1 = m_tree->find(gr1.chr);
  GenomicIntervalTreeMap::const_iterator ff2 = m_tree->find(gr2.chr);
  if (ff1 == m_tree->end() || ff2 == m_tree->end())
    return 0;

  ff1->second.findOverlapping(gr1.pos1, gr1.pos2, giv1);
  ff2->second.findOverlapping(gr2.pos1, gr2.pos2, giv2);

  if (giv1.size() == 0 || giv2.size() == 0)
    return false;
  
  // each one only overlapped one element
  if (giv1.size() == 1 && giv2.size() == 1)
    if (giv1[0].value == giv2[0].value)
      return true;

  // make a set of the possible starts
  std::unordered_set<int> vals;
  for (auto& i : giv1)
    vals.insert(i.value);
  
  // loop the other side and see if they mix
  for (auto& j : giv2)
    if (vals.count(j.value))
      return true;

  return false;
  
  }

template<class T>
std::string GenomicRegionCollection<T>::AsBEDString() const {
  
  if (m_grv->size() ==  0)
    return ""; 

  std::stringstream ss;
  for (auto& i : *m_grv)
    ss << i.chr << "\t" << i.pos1 << "\t" << i.pos2 << std::endl;

  return ss.str();

}

template<class T>
void GenomicRegionCollection<T>::Concat(const GenomicRegionCollection<T>& g)
{
  m_grv->insert(m_grv->end(), g.m_grv->begin(), g.m_grv->end());
}

template<class T>
bool GenomicRegionCollection<T>::GetNextGenomicRegion(T& gr)
{
  if (idx >= m_grv->size())
    return false;

  gr = m_grv->at(++idx);
  return true;
  
}

template<class T>
GenomicRegionCollection<T>::GenomicRegionCollection() {
  __allocate_grc();
}

template<class T>
GenomicRegionCollection<T>::~GenomicRegionCollection() {
}


template<class T>
void GenomicRegionCollection<T>::__allocate_grc() {
  m_grv = std::shared_ptr<std::vector<T>>(new std::vector<T>()) ;
  m_tree = std::shared_ptr<GenomicIntervalTreeMap>(new GenomicIntervalTreeMap()) ;
}

template<class T>
GenomicRegionCollection<T>::GenomicRegionCollection(const BamRecordVector& brv) {

  __allocate_grc();

  for (auto& i : brv) 
    m_grv->push_back(GenomicRegion(i.ChrID(), i.Position(), i.PositionEnd()));

}

template<class T>
const T& GenomicRegionCollection<T>::at(size_t i) const
{ 
  if (i >= m_grv->size()) 
    throw 20;
  return m_grv->at(i); 
}  


// this is query
template<class T>
template<class K>
GenomicRegionCollection<GenomicRegion> GenomicRegionCollection<T>::FindOverlaps(const K& gr, bool ignore_strand) const
{  

  GenomicRegionCollection<GenomicRegion> output;

  if (m_tree->size() == 0 && m_grv->size() != 0) {
    std::cerr << "!!!!!! findOverlaps with GRanges: WARNING: Trying to find overlaps on empty tree. Need to run this->createTreeMap() somewhere " << std::endl;
    return output;
  }
  
  // which chr (if any) are common between query and subject
  GenomicIntervalTreeMap::const_iterator ff = m_tree->find(gr.chr);

  //must as least share a chromosome  
  if (ff == m_tree->end())
    return output;

  // get the subject hits
  GenomicIntervalVector giv;
  ff->second.findOverlapping(gr.pos1, gr.pos2, giv);
  
#ifdef DEBUG_OVERLAPS
  std::cerr << "ff->second.intervals.size() " << ff->second.intervals.size() << std::endl;
  for (auto& k : ff->second.intervals)
    std::cerr << " intervals " << k.start << " to " << k.stop << " value " << k.value << std::endl;
  std::cerr << "GIV NUMBER OF HITS " << giv.size() << " for query " << gr << std::endl;
#endif

  // loop through the hits and define the GenomicRegion
  for (auto& j : giv) { // giv points to positions on subject
    if (ignore_strand || (m_grv->at(j.value).strand == gr.strand) ) {
#ifdef DEBUG_OVERLAPS
	std::cerr << "find overlaps hit " << j.start << " " << j.stop << " -- " << j.value << std::endl;
#endif
	output.add(GenomicRegion(gr.chr, std::max(static_cast<int32_t>(j.start), gr.pos1), std::min(static_cast<int32_t>(j.stop), gr.pos2)));
      }
  }

  return output;
  
}

  // this is query
  template<class T>
  template<class K>
GenomicRegionCollection<GenomicRegion> GenomicRegionCollection<T>::FindOverlaps(GenomicRegionCollection<K>& subject, std::vector<int32_t>& query_id, std::vector<int32_t>& subject_id, bool ignore_strand) const
{  

  GenomicRegionCollection<GenomicRegion> output;
  if (subject.NumTree() == 0 && subject.size() != 0) {
    std::cerr << "!!!!!! findOverlaps: WARNING: Trying to find overlaps on empty tree. Need to run this->createTreeMap() somewhere " << std::endl;
    return output;
  }

  // we loop through query, so want it to be smaller
  if (subject.size() < m_grv->size() && m_grv->size() - subject.size() > 20) 
    std::cerr << "findOverlaps warning: Suggest switching query and subject for efficiency." << std::endl;

#ifdef DEBUG_OVERLAPS
  std::cerr << "OVERLAP SUBJECT: " << std::endl;
  for (auto& i : subject)
    std::cerr << i << std::endl;
#endif

  // loop through the query GRanges (this) and overlap with subject
  for (size_t i = 0; i < m_grv->size(); ++i) 
    {
      // which chr (if any) are common between query and subject
      GenomicIntervalTreeMap::const_iterator ff = subject.GetTree()->find(m_grv->at(i).chr);

      GenomicIntervalVector giv;

#ifdef DEBUG_OVERLAPS
      std::cerr << "TRYING OVERLAP ON QUERY " << m_grv->at(i) << std::endl;
#endif
      //must as least share a chromosome
      if (ff != m_tree->end())
	{
	  // get the subject hits
	  ff->second.findOverlapping(m_grv->at(i).pos1, m_grv->at(i).pos2, giv);

#ifdef DEBUG_OVERLAPS
	  std::cerr << "ff->second.intervals.size() " << ff->second.intervals.size() << std::endl;
	  for (auto& k : ff->second.intervals)
	    std::cerr << " intervals " << k.start << " to " << k.stop << " value " << k.value << std::endl;
	  std::cerr << "GIV NUMBER OF HITS " << giv.size() << " for query " << m_grv->at(i) << std::endl;
#endif
	  // loop through the hits and define the GenomicRegion
	  for (auto& j : giv) { // giv points to positions on subject
	    if (ignore_strand || (subject.at(j.value).strand == m_grv->at(i).strand) )
	      {
		query_id.push_back(i);
		subject_id.push_back(j.value);
#ifdef DEBUG_OVERLAPS
		std::cerr << "find overlaps hit " << j.start << " " << j.stop << " -- " << j.value << std::endl;
#endif
		output.add(GenomicRegion(m_grv->at(i).chr, std::max(static_cast<int32_t>(j.start), m_grv->at(i).pos1), std::min(static_cast<int32_t>(j.stop), m_grv->at(i).pos2)));
	      }
	  }
	}
    }

  return output;
  
}


template<class T>
GenomicRegionCollection<T>::GenomicRegionCollection(const T& gr)
{
  __allocate_grc();
  m_grv->push_back(gr);
}

template<class T>
GRC GenomicRegionCollection<T>::intersection(GRC& subject, bool ignore_strand /* false */)
{
  std::vector<int32_t> sub, que;
  GRC out = this->FindOverlaps(subject, que, sub, ignore_strand);
  return out;
}

template<class T>
void GenomicRegionCollection<T>::Pad(int v)
{
  for (auto& i : *m_grv)
    i.Pad(v);
}

}

