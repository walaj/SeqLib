//#include "SnowTools/GenomicRegionCollection.h"

#include <iostream>
#include <sstream>
#include <cassert>
#include "SnowTools/gzstream.h"

namespace SnowTools {

template<class T>
void GenomicRegionCollection<T>::readMuTect(const std::string &file, int pad) {
  
  assert(pad >= 0);

  std::cout << "Reading MuTect CallStats"  << std::endl;
  std::string curr_chr = "dum";
  
  igzstream iss(file.c_str());
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
	    std::cout << "...reading MuTect call-stats -- chr" << chr << std::endl;
	    curr_chr = chr;
	  }

	}
	if (GenomicRegion::chrToNumber(chr) >= 0) {
	  T gr(chr, pos, pos);
	  gr.pad(pad);
	  m_grv.push_back(gr);
	}
      } // end "keep" conditional
    } // end main while

}

template<class T>
void GenomicRegionCollection<T>::readBEDfile(const std::string & file, int pad) {

  assert(pad >= 0);

  igzstream iss(file.c_str());
  if (!iss || file.length() == 0) { 
    std::cerr << "BED file does not exist: " << file << std::endl;
    exit(EXIT_FAILURE);
  }

  std::string line;
  std::cout << "Reading normal BED" << std::endl;
  std::string curr_chr = "-1";
  while (std::getline(iss, line, '\n')) {

    size_t counter = 0;
    std::string chr, pos1, pos2;
    std::istringstream iss_line(line);
    std::string val;
    
    if (line.find("#") == std::string::npos) {
      while(std::getline(iss_line, val, '\t')) {
	switch (counter) { 
	case 0 : chr = SnowTools::scrubString(val, "chr"); break; 
	case 1 : pos1 = val; break;
	case 2 : pos2 = val; break;
	}
	if (counter >= 2)
	  break;
	counter++;
	
	if (chr != curr_chr) {
	  //std::cout << "...reading from BED - chr" << chr << std::endl;
	  curr_chr = chr;
	}
	
      }
      if (GenomicRegion::chrToNumber(chr) >= 0) {
	T gr(chr, pos1, pos2);
	gr.pad(pad);
	m_grv.push_back(gr);
      }
    } // end "keep" conditional
  } // end main while

  
}

template<class T>
void GenomicRegionCollection<T>::readVCFfile(const std::string & file, int pad) {

  assert(pad >= 0);

  igzstream iss(file.c_str());
  if (!iss || file.length() == 0) { 
    std::cerr << "VCF file does not exist: " << file << std::endl;
    exit(EXIT_FAILURE);
  }

  std::string line;
  
  std::cout << "Parsing VCF file "  << std::endl;
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
	  T gr(chr, pos, pos);
	  gr.pad(pad);
	  m_grv.push_back(gr);
	}
	
      }
    }
  }

  
}

template<class T>
void GenomicRegionCollection<T>::regionFileToGRV(const std::string &file, int pad) {

  igzstream iss(file.c_str());
  if (!iss || file.length() == 0) { 
    std::cerr << "Region file does not exist: " << file << std::endl;
    exit(EXIT_FAILURE);
  }

  // get the header line to check format
  std::string header;
  if (!std::getline(iss, header, '\n'))
    std::cerr << "Region file is empty: " << file << std::endl;
  iss.close();

  GenomicRegionCollection<T> grv;

  // MUTECT CALL STATS
  if (header.find("MuTect") != std::string::npos)
    readMuTect(file, pad);
  // BED file
  else if (file.find(".bed") != std::string::npos)
    readBEDfile(file, pad);
  // VCF file
  else if (file.find(".vcf") != std::string::npos) 
    readVCFfile(file, pad);
  else // default is BED file
    readBEDfile(file, pad);    
}

// reduce a set of GenomicRegions into the minium overlapping set (same as GenomicRanges "reduce")
template <class T>
void GenomicRegionCollection<T>::mergeOverlappingIntervals() {

  // make the list
  std::list<T> intervals(m_grv.begin(), m_grv.end());

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
  std::vector<T> v{ std::make_move_iterator(std::begin(intervals)), 
      std::make_move_iterator(std::end(intervals)) };
  m_grv = v;

}

template <class T>
void GenomicRegionCollection<T>::createTreeMap() {

  // sort the 
  sort(m_grv.begin(), m_grv.end());

  GenomicIntervalMap map;
  //for (auto it : m_grv) {
  //  map[it.chr].push_back(GenomicInterval(it.pos1, it.pos2, it));
  //}
  for (size_t i = 0; i < m_grv.size(); ++i)
    map[m_grv[i].chr].push_back(GenomicInterval(m_grv[i].pos1, m_grv[i].pos2, i));
    
  for (auto it : map) 
    m_tree[it.first] = GenomicIntervalTree(it.second);

}

template <class T>
void GenomicRegionCollection<T>::sendToBED(const std::string file) {
  
  if (m_grv.size() ==  0) {
    std::cerr << "sendToBED: GenomicRegionCollection is empty" << std::endl;
    return; 
  }

  std::ofstream ofile(file.c_str(), ios::out);
  for (auto it : m_grv)
    ofile << GenomicRegion::chrToString(it.chr) << "\t" << it.pos1 << "\t" << it.pos2 << std::endl;
  ofile.close();

}

// divide a region into pieces of width and overlaps
template<class T>
GenomicRegionCollection<T>::GenomicRegionCollection(int width, int ovlp, const T &gr) {

  // undefined otherwise
  assert(width > ovlp);

  uint32_t start = gr.pos1;
  uint32_t end = gr.pos1 + width;

  // region is smaller than width
  if ( end >= gr.pos2 ) {
    std::cerr << "GenomicRegionCollection constructor: GenomicRegion is smaller than bin width" << std::endl;
    return; 
  }

  // loop through the sizes until done
  while (end <= gr.pos2) {
    m_grv.push_back(T(gr.chr, start, end));
    end += width - ovlp; // make the new one
    start += width - ovlp;
  }
  assert(m_grv.size() > 0);
  
  // finish the last one
  start = m_grv.back().pos2 - width;
  end = gr.pos2;
  m_grv.push_back(T(gr.chr, start, end));

}

template<class T>
size_t GenomicRegionCollection<T>::findOverlapping(const T &gr) {

  if (m_tree.size() == 0 && m_grv.size() != 0) 
    {
      std::cerr << "!!!!!! WARNING: Trying to find overlaps on empty tree. Need to run this->createTreeMap() somewhere " << std::endl;
      return 0;
    }

  GenomicIntervalVector giv;
  GenomicIntervalTreeMap::iterator ff = m_tree.find(gr.chr);
  if (ff == m_tree.end())
    return 0;
  ff->second.findOverlapping(gr.pos1, gr.pos2, giv);
  return (giv.size());


}

template<class T>
string GenomicRegionCollection<T>::sendToBED() const {
  
  if (m_grv.size() ==  0)
    return ""; 

  stringstream ss;
  for (auto& i : m_grv)
    ss << i.chr << "\t" << i.pos1 << "\t" << i.pos2 << endl;

  return ss.str();

}

template<class T>
void GenomicRegionCollection<T>::concat(const GenomicRegionCollection<T>& g)
{
  m_grv.insert(m_grv.begin(), g.m_grv.begin(), g.m_grv.end());
}

template<class T>
bool GenomicRegionCollection<T>::getNextGenomicRegion(T& gr)
{
  if (idx >= m_grv.size())
    return false;

  gr = m_grv[++idx];
  return true;
  
}

template<class T>
GenomicRegionCollection<T>::GenomicRegionCollection(std::vector<T>& vec) : m_grv(vec)
{
}

template<class T>
const T& GenomicRegionCollection<T>::at(size_t i) const
{ 
  if (i >= m_grv.size()) 
    throw 20;
  return m_grv[i]; 
}

template<class T>
GenomicRegionCollection<GenomicRegion> GenomicRegionCollection<T>::findOverlaps(GenomicRegionCollection<GenomicRegion>& subject, std::vector<size_t>& query_id, std::vector<size_t>& subject_id, bool ignore_strand)
{  

  GenomicRegionCollection<GenomicRegion> output;

  // loop through the query GRanges (this) and overlap with subject
  for (size_t i = 0; i < m_grv.size(); ++i) 
    {
      GenomicIntervalVector giv;
      GenomicIntervalTreeMap::iterator ff = subject.m_tree.find(m_grv[i].chr);

      std::cout << "TRYING OVERLAP ON QUERY " << m_grv[i] << std::endl;

      //must as least share a chromosome
      if (ff != m_tree.end())
	{
	  // get the subject hits
	  ff->second.findOverlapping(m_grv[i].pos1, m_grv[i].pos2, giv);

	  //debug
	  for (auto& k : ff->second.intervals)
	    std::cout << " intervals " << k.start << " to " << k.stop << " value " << k.value << std::endl;

	  std::cout << "GIV NUMBER OF HITS " << giv.size() << " for query " << m_grv[i] << std::endl;

	  // loop through the hits and define the GenomicRegion
	  for (auto& j : giv) { // giv points to positions on subject
	    if (!ignore_strand || (subject.m_grv[j.value].strand == m_grv[i].strand) )
	      {
		query_id.push_back(i);
		subject_id.push_back(j.value);
		std::cout << "find overlaps hit " << j.start << " " << j.stop << " -- " << j.value << std::endl;
		output.add(T(m_grv[i].chr, max(j.start, (int)m_grv[i].pos1), min(j.stop, (int)m_grv[i].pos2)));
	      }
	  }
	}
    }
  
  return output;
  
}

template<class T>
GenomicRegionCollection<T>::GenomicRegionCollection(const T& gr)
{
  m_grv.push_back(gr);
}

template<class T>
GRC GenomicRegionCollection<T>::intersection(GRC& subject)
{
  std::vector<size_t> sub, que;
  GRC out = this->findOverlaps(subject, que, sub);
  return out;
}


}
