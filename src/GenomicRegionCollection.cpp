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
  void GenomicRegionCollection<T>::gsort() {
    
    std::sort(m_grv->begin(), m_grv->end());
  }

  template<class T>
  void GenomicRegionCollection<T>::SortAndStretchRight(int max) {

    if (!m_grv->size())
      return;
    
    std::sort(m_grv->begin(), m_grv->end());

    if (max < m_grv->back().pos2) {
      std::cerr << "GenomicRegionCollection::SortAndStrech Can't stretch to max, as we are already past max. Max: " << max << " highest GRC " << m_grv->back() << std::endl;
      exit(EXIT_FAILURE);
    }

    for (size_t i = 0; i < m_grv->size() - 1; ++i) {
      m_grv->at(i).pos2 = m_grv->at(i+1).pos1 - 1;
    }
    m_grv->back().pos2 = max;
    
  }

  template<class T>
  void GenomicRegionCollection<T>::SortAndStretchLeft(int min) {

    if (!m_grv->size())
      return;

    std::sort(m_grv->begin(), m_grv->end());

    if (min < m_grv->begin()->pos1) {
      std::cerr << "GenomicRegionCollection::SortAndStrechLeft Can't stretch to min, as we are already past max" << std::endl;
      exit(EXIT_FAILURE);
    }

    m_grv->at(0).pos1 = min;
    for (size_t i = 1; i < m_grv->size(); ++i) {
      m_grv->at(i).pos1 = m_grv->at(i-1).pos2 + 1;
    }
    
  }



template<class T>
void GenomicRegionCollection<T>::readMuTect(const std::string &file, int pad, const BamHeader& hdr) {
  
  assert(pad >= 0);

  std::cerr << "Reading MuTect CallStats"  << std::endl;
  std::string curr_chr = "dum";
  
  std::ifstream iss(file.c_str());
  if (!iss || file.length() == 0) { 
    std::cerr << "MuTect call-stats file does not exist: " << file << std::endl;
    exit(EXIT_FAILURE);
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
	  counter++;
	  
	  if (curr_chr != chr) {
	    std::cerr << "...reading MuTect call-stats -- chr" << chr << std::endl;
	    curr_chr = chr;
	  }

	}

	// parse the strings and send to genomci region
	T gr(chr, pos, pos, hdr);
	if (gr.valid()) {
	  gr.pad(pad);
	  m_grv->push_back(gr);
	}

      } // end "keep" conditional
    } // end main while
}

template<class T>
void GenomicRegionCollection<T>::readBEDfile(const std::string & file, int pad, const BamHeader& hdr) {

  assert(pad >= 0);

  std::ifstream iss(file.c_str());
  if (!iss || file.length() == 0) { 
    std::cerr << "BED file does not exist: " << file << std::endl;
    exit(EXIT_FAILURE);
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
	counter++;
	
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

      if (gr.valid()) {
	gr.pad(pad);
	m_grv->push_back(gr);
      }
	
	//}
    } // end "keep" conditional
  } // end main while

  
}

template<class T>
void GenomicRegionCollection<T>::readVCFfile(const std::string & file, int pad, const BamHeader& hdr) {

  assert(pad >= 0);

  std::ifstream iss(file.c_str());
  if (!iss || file.length() == 0) { 
    std::cerr << "VCF file does not exist: " << file << std::endl;
    exit(EXIT_FAILURE);
  }

  std::string line;
  
  std::cerr << "Parsing VCF file "  << std::endl;
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
	} else {

	  // construct the GenomicRegion
	  T gr(chr, pos, pos, hdr);
	  if (gr.valid()) {
	    gr.pad(pad);
	    m_grv->push_back(gr);
	  }
	}
	
      }
    }
  }
  
  std::cerr << "...parsed " << m_grv->size() << " VCF lines " << std::endl;
  
}

template<class T>
GenomicRegionCollection<T>::GenomicRegionCollection(const std::string &file, const BamHeader& hdr) {

  std::ifstream iss(file.c_str());
  if (!iss || file.length() == 0) { 
    std::cerr << "Region file does not exist: " << file << std::endl;
    exit(EXIT_FAILURE);
  }

  // get the header line to check format
  std::string header;
  if (!std::getline(iss, header, '\n'))
    std::cerr << "Region file is empty: " << file << std::endl;
  iss.close();

  const int pad = 0;
  // MUTECT CALL STATS
  if ((header.find("MuTect") != std::string::npos) ||
      (file.find("call_stats") != std::string::npos) || 
      (file.find("callstats") != std::string::npos))
    readMuTect(file, pad, hdr);
  // BED file
  else if (file.find(".bed") != std::string::npos)
    readBEDfile(file, pad, hdr);
  // VCF file
  else if (file.find(".vcf") != std::string::npos) 
    readVCFfile(file, pad, hdr);
  else // default is BED file
    readBEDfile(file, pad, hdr);    
}

// reduce a set of GenomicRegions into the minium overlapping set (same as GenomicRanges "reduce")
template <class T>
void GenomicRegionCollection<T>::mergeOverlappingIntervals() {

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
  m_grv = nullptr; // clear the old data 
  std::vector<T> v{ std::make_move_iterator(std::begin(intervals)), 
      std::make_move_iterator(std::end(intervals)) };
  m_grv->insert(m_grv->end(), v.begin(), v.end());
  //m_grv = v;

}

template <class T>
GenomicRegionVector GenomicRegionCollection<T>::asGenomicRegionVector() const { 
  GenomicRegionVector gg;
  for (auto& i : *m_grv)
    gg.push_back(GenomicRegion(i.chr, i.pos1, i.pos2, i.strand));
  return gg; 
} 
  
template <class T>
void GenomicRegionCollection<T>::createTreeMap() {

  // sort the 
  sort(m_grv->begin(), m_grv->end());

  GenomicIntervalMap map;
  //for (auto it : *m_grv) {
  //  map[it.chr].push_back(GenomicInterval(it.pos1, it.pos2, it));
  //}
  for (size_t i = 0; i < m_grv->size(); ++i) {
    map[m_grv->at(i).chr].push_back(GenomicInterval(m_grv->at(i).pos1, m_grv->at(i).pos2, i));
  }

  for (auto it : map) {
    GenomicIntervalTreeMap::iterator ff = m_tree->find(it.first);
    ff->second = it.second;
    //old //m_tree[it.first] = GenomicIntervalTree(it.second);
  }

}

template <class T>
void GenomicRegionCollection<T>::sendToBED(const std::string file) {
  
  if (m_grv->size() ==  0) {
    std::cerr << "sendToBED: GenomicRegionCollection is empty" << std::endl;
    return; 
  }

  std::ofstream ofile(file.c_str(), std::ios::out);
  for (auto it : *m_grv)
    ofile << GenomicRegion::chrToString(it.chr) << "\t" << it.pos1 << "\t" << it.pos2 << "\t" << it.strand << std::endl;
  ofile.close();

}

// divide a region into pieces of width and overlaps
template<class T>
GenomicRegionCollection<T>::GenomicRegionCollection(int width, int ovlp, const T &gr) {

  __allocate_grc();

  // undefined otherwise
  assert(width > ovlp);
  if (width >= gr.width()) {
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
size_t GenomicRegionCollection<T>::findOverlapping(const T &gr) const {

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
  bool GenomicRegionCollection<T>::overlapSameBin(const T &gr1, const T &gr2) const {

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
std::string GenomicRegionCollection<T>::sendToBED() const {
  
  if (m_grv->size() ==  0)
    return ""; 

  std::stringstream ss;
  for (auto& i : *m_grv)
    ss << i.chr << "\t" << i.pos1 << "\t" << i.pos2 << std::endl;

  return ss.str();

}

template<class T>
void GenomicRegionCollection<T>::concat(const GenomicRegionCollection<T>& g)
{
  m_grv->insert(m_grv->end(), g.m_grv->begin(), g.m_grv->end());
}

template<class T>
bool GenomicRegionCollection<T>::getNextGenomicRegion(T& gr)
{
  if (idx >= m_grv->size())
    return false;

  gr = m_grv->at(++idx);
  return true;
  
}

template<class T>
GenomicRegionCollection<T>::GenomicRegionCollection() {
  __allocate_grc();
};

template<class T>
GenomicRegionCollection<T>::~GenomicRegionCollection() {
};


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
GenomicRegionCollection<GenomicRegion> GenomicRegionCollection<T>::findOverlaps(const K& gr, bool ignore_strand) const
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
    if (ignore_strand || (m_grv->at(j.value).strand == gr.strand) )
      {
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
GenomicRegionCollection<GenomicRegion> GenomicRegionCollection<T>::findOverlaps(GenomicRegionCollection<K>& subject, std::vector<int32_t>& query_id, std::vector<int32_t>& subject_id, bool ignore_strand) const
{  

  GenomicRegionCollection<GenomicRegion> output;
  if (subject.m_tree->size() == 0 && subject.m_grv->size() != 0) {
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
      GenomicIntervalTreeMap::const_iterator ff = subject.m_tree->find(m_grv->at(i).chr);

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
	    if (ignore_strand || (subject.m_gr->at(j.value).strand == m_grv->at(i).strand) )
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
  GRC out = this->findOverlaps(subject, que, sub, ignore_strand);
  return out;
}

  /*
template<class T>
GRC GenomicRegionCollection<T>::complement(GRC& subject, bool ignore_strand)
{

  GRC out;
  // get this - that
  
#ifdef BOOST_CONFIG_HPP_dum
  using namespace boost::icl;
  typedef interval_set<int> TIntervalSet; 
  
  // flatten each
  subject.mergeOverlappingIntervals();
  this->mergeOverlappingIntervals();

  // get list of THAT chr to go through
  std::set<int> chr;
  for (auto& i : subject)
    if (!chr.count(i.chr))
      chr.insert(i.chr);

  // the THIS events that have chr not in THAT, automatically get added
  for (auto& i : m_grv)
    if (!chr.count(i.chr))
      out.add(i);

  // loop through the chromosomes of THAT to subtract
  for (auto& c : chr) 
    {
      TIntervalSet intSetA;
      TIntervalSet intSetB;
  
      // build up the intervals for THIS
      for (auto& i : m_grv)
	if (i.chr == c)
	  intSetA += discrete_interval<int>::closed(i.pos1, i.pos2);
      // build up the intervals for THAT
      for (auto& i : subject)
	if (i.chr == c)
	  intSetB += discrete_interval<int>::closed(i.pos1, i.pos2);
      
      TIntervalSet intSetD = intSetA - intSetB;

      for (auto& i : intSetD)
	out.add(GenomicRegion(c, first(i), upper(i)));
    }

  return out;

#else
  std::cerr << "No Boost library detected. GenomicRegionCollection::complement will not function. Returning Empty GRC" << std::endl;

  // sort this
  //gsort();

  // sort that
  //subject.gsort();

  subject.createTreeMap();
  createTreeMap();

  std::vector<int32_t> sub, que;
  GRC tt = findOverlaps(subject, sub, que, true);
  for (auto& i : tt)
    std::cerr << i << std::endl;
  
  
  return out;
#endif

}
  */

template<class T>
void GenomicRegionCollection<T>::pad(int v)
{

  for (auto& i : *m_grv)
    i.pad(v);

}

}

