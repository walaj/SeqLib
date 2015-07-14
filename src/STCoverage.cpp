#include "SnowTools/STCoverage.h"
#include "SnowTools/SnowToolsCommon.h"
#include <stdexcept>

namespace SnowTools {

  STCoverage::STCoverage(const SnowTools::GenomicRegion& gr) {
    std::cerr << "*********** making cov *************** " << std::endl;//debug
    m_gr = gr;
    v = uint16_sp(new std::vector<uint16_t>(gr.width(),0));
  }

  void STCoverage::addRead(const BamRead &r) {
    
    int p = r.Position() - m_gr.pos1;
    int e = r.PositionEnd() - m_gr.pos1;
    if (p < 0 || e < 0)
      return;
    
    try {
      while (p <= e) {
	if (v->at(p) < 60000) // 60000 is roughly int16 lim
	  v->at(p)++;
	++p;
	
      }
    } catch (std::out_of_range &oor) {
      //std::cerr << "Position " << p << " on tid " << r.ChrID()
      //      << " is greater than expected max of " << v->size() << " -- skipping" << std::endl;
      
    }
    
  }

  std::ostream& operator<<(std::ostream &out, const STCoverage &c) {
    out << "Region " << c.m_gr << " v.size() " << c.v->size() << std::endl;
    return out;
  }

  void STCoverage::ToBedgraph(ogzstream * o, const bam_hdr_t * h) const {

    // unitialized so nothing to do
    if (m_gr.chr == -1 || v->size() == 0)
      return;

    size_t curr_start = 0;
    size_t curr_val = v->at(0);
    for (size_t i = 0; i < v->size(); ++i) {
      if (v->at(i) != curr_val) {
	(*o) << m_gr.ChrName(h) << "\t" << (curr_start + m_gr.pos1) << "\t" << (i+m_gr.pos1) << "\t" << curr_val << std::endl;
	curr_start = i;
	curr_val = v->at(i);
      }
    }
    
    // need to dump last one
    if ( (curr_start+1) != v->size()) 
      (*o) << m_gr.ChrName(h) << "\t" << (curr_start + m_gr.pos1) << "\t" << (v->size()+m_gr.pos1-1) << "\t" << curr_val << std::endl;
  }
  
  uint16_t STCoverage::getCoverageAtPosition(size_t pos) const {

  if (pos < m_gr.pos1 || pos > m_gr.pos2) {
    //std::cerr << "Coverage query out of bounds for location " << m_gr.chr << ":" << pos << std::endl;
    return 0;
  }
  
  size_t q = pos - m_gr.pos1;
  if (q >= v->size()) {
    std::cerr << "Coverage query out of bounds for location " << m_gr.chr << ":" << pos << " with pos-start of " << q << " attempt on v of size " << v->size() << std::endl;
    return 0;
  }

  return (v->at(q));
    

}

}
