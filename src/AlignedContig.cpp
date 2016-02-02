#include "SnowTools/AlignedContig.h"

#include <unordered_map>

#define MAX_CONTIG_SIZE 5000000

namespace SnowTools {

  AlignedContig::AlignedContig(const BamReadVector& bav) {
    
    if (!bav.size())
      return;
    
    // set the sequence. Convention is store as it came off assembler for first alignment
    m_seq = bav.begin()->Sequence();
    if (bav.begin()->ReverseFlag()) {
      SnowTools::rcomplement(m_seq);
    }

    // zero the coverage
    for (auto& i : cov)
      i.second = std::vector<int>(m_seq.length(), 0);
    
    // find the number of primary alignments
    size_t num_align = 0;
    for (auto& i : bav)
      if (!i.SecondaryFlag())
	++num_align;

    // make the individual alignments and add
    for (auto& i : bav) {
      if (!i.SecondaryFlag()) {
	bool flip = (m_seq != i.Sequence()); // if the seq was flipped, need to flip the AlignmentFragment
	m_frag_v.push_back(AlignmentFragment(i, flip, prefixes));
	m_frag_v.back().num_align = num_align;
      } else {
	bool flip = (m_seq != i.Sequence()); // if the seq was flipped, need to flip the AlignmentFragment
	if (m_frag_v.size())
	  m_frag_v.back().secondaries.push_back(AlignmentFragment(i, flip, prefixes));
	//m_frag_v_secondary.push_back(AlignmentFragment(i, flip));
	//m_frag_v_secondary.back().num_align = bav.size();
      }      
    }

    // sort fragments by order on fwd-strand contig
    if (m_frag_v.size() > 1)
      std::sort(m_frag_v.begin(), m_frag_v.end());

    if (!m_frag_v.size())
      return;

    // get breaks out of it
    setMultiMapBreakPairs();
  }

  void AlignedContig::printContigFasta(std::ofstream& os) const {
    os << ">" << getContigName() << std::endl;
    os << getSequence() << std::endl;
  }
  
  std::string AlignedContig::print() const {

    std::stringstream out;
    
    // print the AlignmentFragments alignments
    for (auto& i : m_frag_v) 
      out << i.print() << getContigName() << std::endl;
    bool draw_divider = true;
    for (auto& i : m_frag_v) 
      for (auto& j : i.secondaries) {
	if (draw_divider) {
	  out << std::string(m_seq.length(), BAM_CSOFT_CLIP) << std::endl;
	  draw_divider = false;
	}
      }
    
    // print the break locations for indel deletions
    for (auto& i : m_frag_v) {
      for (auto& j : i.getIndelBreaks()) {
	if (j.num_align == 1 && j.insertion == "") // deletion
	  out << std::string(j.b1.cpos, ' ') << "|" << std::string(j.b2.cpos-j.b1.cpos-1, ' ') << '|' << "   " << getContigName() << std::endl;	
      }
    }
    
    out << getSequence() << "    " << getContigName() << std::endl; 
    PlottedReadVector plot_vec;
    
    // print out the individual reads
    for (auto& i : m_bamreads) {
      
      int pos = -1;
      int aln = -1;
      int rc = 0;
      std::string this_cig;
      //bool was_trimmed = false;
      //std::string seq = i.QualityTrimmedSequence(4, dum, was_trimmed);
      //std::string seq = ""; // dummy
      std::string seq = i.Sequence(); //i.QualitySequence(); // they don't have qual scores, since they are artifically created from alignSingleSequence
      std::string sr = i.GetZTag("SR");
      
      // get the more complex tags (since there can be multiple annotations per tag)
      std::vector<int> posvec = i.GetSmartIntTag("SL"); // start positions ON CONTIG
      std::vector<int> alnvec = i.GetSmartIntTag("TS"); // start positions ON READ
      std::vector<int> rcvec = i.GetSmartIntTag("RC"); // read reverse complmented relative to contig
      std::vector<std::string> cigvec = i.GetSmartStringTag("SC"); // read against contig CIGAR
      std::vector<std::string> cnvec = i.GetSmartStringTag("CN");
      std::vector<std::string> hcvec = i.GetSmartStringTag("HC"); // hard-clipped seq on CONTIG      

      if (posvec.size() != alnvec.size() ||
	  posvec.size() != rcvec.size() ||
	  cigvec.size() != posvec.size() ||
	  cnvec.size() != posvec.size())
	continue;
      
      std::string hc;
      assert(cnvec.size() == posvec.size());
      size_t kk = 0;
      for (; kk < cnvec.size(); kk++) 
	if (cnvec[kk] == getContigName()) {
	  pos = posvec[kk];
	  aln = alnvec[kk];
	  rc = rcvec[kk];
	  this_cig = cigvec[kk];
	  if (hcvec.size()) {
	    hc = hcvec[kk];
	    // bump up the alignments if starts with hardclip
	    boost::regex reg("([0-9]+)H.*"); 
	    boost::cmatch match;
	    if (boost::regex_match(this_cig.c_str(), match, reg))
	      pos += std::stoi(match[1].str());
	  }
	  break;
	}
      
      // if we have harclip of read to contig, use that
      if (!hc.empty())
	seq = hc;

      // reverse complement if need be
      if (rc)
	SnowTools::rcomplement(seq);      
      
      if (aln > 0)
	try {
	  seq = seq.substr(aln, seq.length() - aln);
	} catch (...) {
	  std::cerr << "AlignedContig::operator<< error: substring out of bounds. seqlen " << 
	    seq.length() << " start " << aln << " length " << (seq.length() - aln) << std::endl;
	}
      
      if ( (pos + seq.length() ) > getSequence().length()) 
	try { 
	  seq = seq.substr(0, getSequence().length() - pos);
	} catch (...) {
	  std::cerr << "AlignedContig::operator<< (2) error: substring out of bounds. seqlen " << 
	    seq.length() << " start " << 0 << " pos " << pos << " getSequence().length() " << 
	    getSequence().length() << std::endl;

	}

      assert(kk != cnvec.size()); // assure that we found something
      pos = abs(pos);
      int padlen = getSequence().size() - pos - seq.size() + 5;
      padlen = std::max(5, padlen);
      
      std::stringstream rstream;
      assert(pos < MAX_CONTIG_SIZE && padlen < MAX_CONTIG_SIZE); // bug, need to check
      rstream << sr << "--" << (i.ChrID()+1) << ":" << i.Position() << " r2c CIGAR: " << this_cig;
      
      plot_vec.push_back({pos, seq, rstream.str()});
    }
    
    std::sort(plot_vec.begin(), plot_vec.end());
    
    PlottedReadLineVector line_vec;
    
    // plot the reads from the ReadPlot vector
    for (auto& i : plot_vec) {
      bool found = false;
      for (auto& j : line_vec) {
	if (j.readFits(i)) { // it fits here
	  j.addRead(&i);
	  found = true;
	  break;
	}
      }
      if (!found) { // didn't fit anywhere, so make a new line
	PlottedReadLine prl;
	prl.contig_len = getSequence().length();
	prl.addRead(&i);
	line_vec.push_back(prl);
      }
    }
    
    // plot the lines. Add contig identifier to each
    for (auto& i : line_vec) 
      out << i << " " << getContigName() << std::endl;
    
    return out.str();
  }
  
  void AlignedContig::setMultiMapBreakPairs() {
    
    // if single mapped contig, nothing to do here
    if (m_frag_v.size() == 1)
      return;
    
    // if single mapped contig, nothing to do here  
    // initialize the breakpoint, fill with basic info
    BreakPoint bp;

    bp.seq = getSequence();
    bp.num_align = m_frag_v.size();
    assert(bp.num_align > 0);
    
    bp.cname = getContigName(); 
    assert(bp.cname.length());
    
    // walk along the ordered contig list and make the breakpoint pairs  
    for (AlignmentFragmentVector::const_iterator it = m_frag_v.begin(); it != m_frag_v.end() - 1; it++) {
      
      AlignmentFragmentVector bwa_hits_1, bwa_hits_2;
      bwa_hits_1.push_back(*it);
      bwa_hits_2.push_back(*(it+1));
      bwa_hits_1.insert(bwa_hits_1.end(), it->secondaries.begin(), it->secondaries.end());
      bwa_hits_2.insert(bwa_hits_2.end(), (it+1)->secondaries.begin(), (it+1)->secondaries.end());
      
      // make all of the local breakpoints
      for (auto& a : bwa_hits_1) {
	for (auto& b : bwa_hits_2) {
	  
	  bp.b1 = a.makeBreakEnd(true);
	  bp.b2 = b.makeBreakEnd(false); 
	  
	  // order the breakpoint
	  bp.order();

	  if (a.m_align.SecondaryFlag() || b.m_align.SecondaryFlag())
	    bp.secondary = true;
	  
	  assert(bp.valid());

	  // add the the vector of breakpoints
	  if (!bp.secondary) 
	    m_local_breaks.push_back(bp);
	  else
	    m_local_breaks_secondaries.push_back(bp);	  
	}
      }
    } // end frag iterator loop
    
    // if this is a double mapping, we are done
    if (m_frag_v.size() == 2) {
      m_global_bp = m_local_breaks[0];
      m_local_breaks.clear();
      for (auto& i : m_local_breaks_secondaries) {
	m_global_bp_secondaries.push_back(i);
      }
      m_local_breaks_secondaries.clear();
      return;
    }
    
    // 3+ mappings. If all good, then don't make the "global"
    // Actually, do make the global
    //bool make_locals = true;
    //for (size_t i = 1; i < m_frag_v.size() - 1; ++i)
    //  if (m_frag_v[i].m_align.MapQuality() < 50)
    //	make_locals = false;
    
    //if (make_locals) { // intermediates are good, so just leave locals as-is
    //  return;
    //}
    
    // TODO support 3+ mappings that contain secondary
    
    // go through alignments and find start and end that reach mapq 
    size_t bstart = MAX_CONTIG_SIZE; //1000 is a dummy
    size_t bend = m_frag_v.size() - 1;

    for (size_t i = 0; i < m_frag_v.size(); i++)
      if (m_frag_v[i].m_align.MapQuality() >= 60) {
	bend = i;
	if (bstart == MAX_CONTIG_SIZE)
	  bstart = i;
      }
    if (bstart == bend || bstart==MAX_CONTIG_SIZE) {
      bstart = 0;
      bend = m_frag_v.size() -1 ;
    }

    assert(bend <= m_frag_v.size());
    assert(bstart <= m_frag_v.size());
    assert(bstart != MAX_CONTIG_SIZE);
    
    // there are 3+ mappings, and middle is not great. Set a global break
    m_global_bp = bp;
    m_global_bp.b1 = m_frag_v[bstart].makeBreakEnd(true);
    m_global_bp.b2 = m_frag_v[bend].makeBreakEnd(false);
    
    // set the strands
    //m_global_bp.gr1.strand = m_frag_v[bstart].align.IsReverseStrand() ? '-' : '+';
    //m_global_bp.gr2.strand = m_frag_v[bend].align.IsReverseStrand()   ? '+' : '-';
    ///////m_global_bp.b1.gr.strand = !m_frag_v[bstart].m_align.ReverseFlag() ? '+' : '-';
    ///////m_global_bp.b2.gr.strand =  m_frag_v[bend].m_align.ReverseFlag() ? '+' : '-';
    
    // order the breakpoint
    m_global_bp.order();

  // clear out the locals
  m_local_breaks.clear();
  //std::cerr << m_global_bp << " " << (*this) << std::endl;
  assert(m_global_bp.valid());
  }
  
  //std::pair<int,int> AlignedContig::getCoverageAtPosition(int pos) const {
    //if (pos <= 0 || pos >= (int)tum_cov.size())
    //  return std::pair<int,int>(0,0);
    
    //return std::pair<int,int>(tum_cov[pos], norm_cov[pos]);
  //}

  //std::unordered_map<std::string, int> AlignedContig::getCoverageAtPosition(int pos) const {
  //for (auto& i : allele) 
      //   i.second.cov[pos];
  //    return std::pair<int,int>(tum_cov[pos], norm_cov[pos]);
  // }

  
  AlignmentFragment::AlignmentFragment(const BamRead &talign, bool flip, const std::unordered_set<std::string>& prefixes) {

    m_align = talign;

    sub_n = talign.GetIntTag("SB");

    // orient cigar so it is on the contig orientation. 
    // need to do this to get the right ordering of the contig fragments below
    // We only flip if we flipped the sequence, and that was determined
    // by the convention set in AlignedContig::AlignedContig, so we need
    // to pass that information explicitly
    if (flip/*m_align.ReverseFlag()*/) {
      m_cigar = m_align.GetReverseCigar();
    } else { 
      m_cigar = m_align.GetCigar();
    }

    // find the start position of alignment ON CONTIG
    start = 0; 
    for (auto& i : /*m_align.GetCigar()*/ m_cigar) {
      if (i.RawType() != BAM_CMATCH)
	start += i.Length();
      else
	break;
    }

    // set the left-right breaks
    unsigned currlen  = 0; 
    
    // CPOS is zero based

    // cigar is oriented to as is from aligner
    for (auto& i : m_cigar/*m_align.GetCigar()*/ /*align.CigarData*/) { //CigarOpVec::const_iterator j = align.cigar.begin(); j != align.cigar.end(); j++) {
      
      // SET THE CONTIG BREAK (treats deletions and leading S differently)
      // the first M gets the break1, pos on the left
      if (i.RawType() == BAM_CMATCH && break1 == -1)
	break1 = currlen;
      if (i.RawType() != BAM_CDEL) // m_skip deletions but not leading S, but otherwise update
	currlen += i.Length();
      if (i.RawType() == BAM_CMATCH) // keeps triggering every M, with pos at the right
	break2 = currlen;
    }

    // assign the genomic coordinates of the break
    if (m_align.ReverseFlag()) {
      gbreak2 = m_align.Position() + 1;
      gbreak1 = m_align.PositionEnd();
    } else {
      gbreak1 = m_align.Position() + 1;
      gbreak2 = m_align.PositionEnd();
    }

    if (break1 >= MAX_CONTIG_SIZE || break2 >= MAX_CONTIG_SIZE || break1 < 0 || break2 < 0) 
      std::cerr << " break1 " << break1 << " break2 " << break2 << " " << this->print() << std::endl;
    assert(break1 < MAX_CONTIG_SIZE);
    assert(break2 < MAX_CONTIG_SIZE);

    if (break1 < 0 || break2 < 0) {
      std::cout << this->print() << std::endl;
    }
    assert(break1 >= 0);
    assert(break2 >= 0);

    // parse right away to see if there are indels on this alignment
    BreakPoint bp;
    size_t fail_safe_count = 0;
    while (parseIndelBreak(bp) && fail_safe_count++ < 100 && !m_align.SecondaryFlag()) {
      //std::cerr << bp << " " << (*this) << std::endl;
      assert(bp.valid());
      //for (auto& i : prefixes)
      //bp.allele[i].indel = true;
      m_indel_breaks.push_back(bp);
      assert(bp.num_align == 1);
    }

    assert(fail_safe_count != 100);

}

  std::string AlignmentFragment::print() const {
    
    std::stringstream out;

    // sets the direction to print
    char jsign = '>'; 
    if (m_align.ReverseFlag())
      jsign = '<';
    
    // print the cigar value per base
    for (auto& j : /*m_align.GetCigar()*/ m_cigar) { //align.CigarData) { // print releative to forward strand
      if (j.RawType() == BAM_CMATCH)
	out << std::string(j.Length(), jsign);
      else if (j.RawType() == BAM_CINS) 
	out << std::string(j.Length(), BAM_CINS);
      else if (j.RawType() == BAM_CSOFT_CLIP || j.RawType() == BAM_CHARD_CLIP)
	out << std::string(j.Length(), '.');
    }
    
    // print contig and genome breaks
    out << "\tC[" << break1 << "," << break2 << "] G[" << gbreak1 << "," << gbreak2 << "]";
    
    // add local info
    std::string chr_name = m_align.GetZTag("MC");
    if (!chr_name.length())
      chr_name = std::to_string(m_align.ChrID()+1);
    out << "\tAligned to: " << chr_name << ":" << m_align.Position() << "(" << (m_align.ReverseFlag() ? "-" : "+") << ") CIG: " << m_align.CigarString() << " MAPQ: " << m_align.MapQuality() << " SUBN " << sub_n << " ";
  
    return out.str();
  }
  
  bool AlignmentFragment::parseIndelBreak(BreakPoint &bp) {
    
     // make sure we have a non-zero cigar
    if (m_cigar.size() == 0) {
      std::cerr << "CIGAR of length 0 on " << this->print() << std::endl;
      return false;
    }
 
    // reject if too many mismatches
    //size_t di_count = 0;
    for (auto& i : m_cigar)
      if (i.RawType() == BAM_CDEL || i.RawType() == BAM_CINS)
    	++di_count;
    if (di_count > 2 && m_align.Qname().substr(0,2) == "c_") // only trim for snowman assembled contigs
      return false;

    // reject if it has small matches, could get confused. Fix later
    //for (auto& i : m_cigar) 
    //  if (i.RawType() == BAM_CMATCH && i.Length() < 5)
    //    return false;
    
     // reject if first alignment is I or D or start with too few M
    if (m_cigar.begin()->RawType() == BAM_CINS || m_cigar.begin()->RawType() == BAM_CDEL || m_cigar.back().RawType() == BAM_CDEL || m_cigar.back().RawType() == BAM_CINS) {
      //std::cerr << "rejecting cigar for starting in I or D or < 5 M" << std::endl;
      return false;
    }
    
    // use next available D / I by looping until can increment idx
    size_t loc = 0; // keep track of which cigar field
    for (auto& i : m_cigar) {
      ++loc;
      if (i.RawType() == BAM_CDEL || i.RawType() == BAM_CINS) {
	assert (loc != 1 && loc != m_cigar.size()); // shuldn't start with I or D
	if (loc > idx) {
	  idx = loc;
	  break;
	}
      }
    }
    

    // we made it to the end, no more indels to report
    if (loc == m_cigar.size())
      return false;

    // clear out the old bp just in case
    bp = BreakPoint();
    bp.num_align = 1;
    
    int curr = 0;
    int gcurrlen = 0; 
    
    // make break ends with dummy positions
    bp.b1 = BreakEnd(m_align);
    bp.b2 = BreakEnd(m_align);

    // assign the contig-wide properties
    bp.cname = m_align.Qname();
    bp.seq = m_align.Sequence();
    assert(bp.cname.length());
    
    size_t count = 0; // count to make sure we are reporting the right indel
    for (auto& i : m_cigar) { // want breaks in CONTIG coordinates, so use oriented cigar
      count++;
      
      // set the contig breakpoint
      if (i.RawType() == BAM_CMATCH || i.RawType() == BAM_CINS || i.RawType() == BAM_CSOFT_CLIP) 
	curr += i.Length();
      if (i.RawType() == BAM_CDEL && bp.b1.cpos == -1 && count == idx) {
	
	bp.b1.cpos = curr-1;
	bp.b2.cpos = curr;
      } 
      if (i.RawType() == BAM_CINS && bp.b1.cpos == -1 && count == idx) {
	bp.b1.cpos = curr - i.Length() - 1; // -1 because cpos is last MATCH
	bp.b2.cpos = curr/* - 1*/; 
	bp.insertion = m_align.Sequence().substr(bp.b1.cpos+1, i.Length()); // +1 because cpos is last MATCH.
      }
      
      // set the genome breakpoint
      if (bp.b1.cpos > 0) {
	if (i.RawType() == BAM_CDEL) {
	  if (!m_align.ReverseFlag()) {
	    bp.b1.gr.pos1 =  m_align.Position() + gcurrlen; // dont count this one//bp.b1.cpos + align.Position; //gcurrlen + align.Position;
	    bp.b2.gr.pos1 = bp.b1.gr.pos1 + i.Length() + 1;
	  } else {
	    bp.b2.gr.pos1 =  (m_align.PositionEnd()-1) - gcurrlen; //bp.b1.cpos + align.Position; //gcurrlen + align.Position;
	    bp.b1.gr.pos1 =  bp.b2.gr.pos1 - i.Length() - 1;
	  }
	} else if (i.RawType() == BAM_CINS) {
	  if (!m_align.ReverseFlag()) {
	    bp.b1.gr.pos1 = m_align.Position() + gcurrlen; //gcurrlen + align.Position;
	    bp.b2.gr.pos1 = bp.b1.gr.pos1 + 1;	
	  } else {
	    // GetEndPosition is 1 too high
	    bp.b2.gr.pos1 = (m_align.PositionEnd()-1) - gcurrlen; //gcurrlen + align.Position;
	    bp.b1.gr.pos1 = bp.b2.gr.pos1 - 1;	
	  }
	}
	break; // already got it, so quit cigar loop
      }
      
      // update the position on the genome
      if (i.RawType() == BAM_CMATCH || i.RawType() == BAM_CDEL) {
	gcurrlen += i.Length();
      } 
      
      
    } // end cigar loop
    
    // set the dummy other end
    bp.b1.gr.pos2 = bp.b1.gr.pos1; 
    bp.b2.gr.pos2 = bp.b2.gr.pos1;
   
    // should have been explicitly ordered in the creation above
    if (!(bp.b1.gr < bp.b2.gr)) {
      //std::cerr << "Warning: something went wrong in indelParseBreaks. Can't order breaks" << std::endl;
      //std::cerr << bp.b1.gr << " " << bp.b2.gr << std::endl;
      //std::cerr << (*this) << std::endl;
      return false;
    }
    //bp.order();

    bp.b1.gr.strand = '+';
    bp.b2.gr.strand = '-';

    // find if it has repeat seq near break
    int rstart = std::max(0, bp.b1.cpos - 7);
    int rlen = std::min(rstart + 7,(int)bp.seq.length() - rstart);
    std::string rr;
    try {
      rr = bp.seq.substr(rstart, rlen);
    } catch (...) { 
      std::cerr << "Caught substring error for string: " << bp.seq << " start " << rstart << " len " << rlen << std::endl;
    }
    assert(bp.valid());
    return true;
  }
  
  std::vector<BreakPoint> AlignedContig::getAllBreakPoints(bool local_restrict) const {
    
    std::vector<BreakPoint> out;
    for (auto& i : m_frag_v) {
      if (i.local || !local_restrict) // only output if alignment is local for indels
	for (auto& k : i.m_indel_breaks)
	  out.push_back(k);
    }
 
    if (!m_global_bp.isEmpty())
      out.push_back(m_global_bp);
    
    out.insert(out.end(), m_local_breaks.begin(), m_local_breaks.end());
    out.insert(out.end(), m_global_bp_secondaries.begin(), m_global_bp_secondaries.end());
    
    return out;
  }
  
  std::vector<BreakPoint> AlignedContig::getAllBreakPointsSecondary() const {
    
    std::vector<BreakPoint> out;
    
    for (auto& i : m_global_bp_secondaries)
      if (!i.isEmpty())
	out.push_back(i);
    
    return out;
  }
  
  
  bool AlignedContig::hasVariant() const { 
    
    if (!m_global_bp.isEmpty())
      return true;

    if (m_local_breaks.size())
      return true; 

    if (m_global_bp_secondaries.size())
      return true; 

    if (m_local_breaks_secondaries.size())
      return true; 

    for (auto& i : m_frag_v)
      if (i.local && i.m_indel_breaks.size())
	return true;
    
    return false;
    
  }
  
  BreakEnd AlignmentFragment::makeBreakEnd(bool left) {
    
    BreakEnd b;

    b.chr_name = m_align.GetZTag("MC");
    assert(b.chr_name.length());
    b.mapq = m_align.MapQuality();
    b.matchlen = m_align.NumMatchBases();
    b.nm = m_align.GetIntTag("NM");
	
    // if alignment is too short, zero the mapq
    // TODO get rid of this by integrating matchlen later
    if (b.matchlen < 30)
      b.mapq = 0;

    if (left) {
      b.gr = SnowTools::GenomicRegion(m_align.ChrID(), gbreak2, gbreak2);
      b.gr.strand = m_align.ReverseFlag() ? '-' : '+'; 
      b.cpos = break2; // take the right-most breakpoint as the first
    } else {
      b.gr = SnowTools::GenomicRegion(m_align.ChrID(), gbreak1, gbreak1);
      b.gr.strand = m_align.ReverseFlag() ? '+' : '-';
      b.cpos = break1;  // take the left-most of the next one
    }

    assert(b.cpos < MAX_CONTIG_SIZE);
    
    return b;
  }
  
  void AlignmentFragment::writeToBAM(BamWalker& bw) { 
      bw.writeAlignment(m_align); 
    } 

  void AlignedContig::writeToBAM(BamWalker& bw) { 
    for (auto& i : m_frag_v) 
      i.writeToBAM(bw);
    } 


  void AlignedContig::writeAlignedReadsToBAM(BamWalker& bw) {
    for (auto& i : m_bamreads)
      bw.writeAlignment(i);
  }

}
