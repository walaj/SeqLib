#include "SnowTools/AlignedContig.h"
#include "SnowTools/BamWalker.h"

#include <regex>
#include <unordered_map>

namespace SnowTools {

  /*void AlignedContig::addAlignment(const BamRead &align) {

  AlignmentFragment tal(align); 
  m_frag_v.push_back(tal);
  
  if (m_frag_v.size() > 1)
  sort(m_frag_v.begin(), m_frag_v.end());

  }*/

  AlignedContig::AlignedContig(const BamReadVector& bav) 
  {
    if (!bav.size())
      return;
    
    // set the sequence. Convention is store as it came off assembler for first alignment
    m_seq = bav.begin()->Sequence();
    if (bav.begin()->ReverseFlag()) {
      SnowTools::rcomplement(m_seq);
    }
    
    // zero the coverage
    tum_cov = std::vector<int>(m_seq.length(), 0);
    norm_cov = std::vector<int>(m_seq.length(), 0);

    // make the individual alignments and add
    for (auto& i : bav) {
      bool flip = (m_seq != i.Sequence()); // if the seq was flipped, need to flip the AlignmentFragment
      m_frag_v.push_back(AlignmentFragment(i, flip));
      m_frag_v.back().num_align = bav.size();
    }

    // sort fragments by order on fwd-strand contig
    if (m_frag_v.size() > 1)
      std::sort(m_frag_v.begin(), m_frag_v.end());

    // get breaks out of it
    setMultiMapBreakPairs();
  }

void AlignedContig::printContigFasta(std::ofstream& os) const {
  os << ">" << getContigName() << std::endl;
  os << getSequence() << std::endl;
}
  
void AlignedContig::blacklist(GRC &grv) {

  //if (m_skip)
  //  return;

 // loop through the indel breaks and blacklist
 for (auto& i : m_frag_v) 
   for (auto& j : i.m_indel_breaks) 
     j.checkBlacklist(grv);

}

void AlignedContig::splitCoverage() { 
  
  //if (m_skip)
  //  return;
  
  for (auto& i : m_frag_v) {
    for (auto& j : i.m_indel_breaks) {
      //std::cout << "calling indel split coverage " << std::endl;
      j.splitCoverage(m_bamreads);
    }
  }
  
  for (auto& i : m_local_breaks) {
    //std::cout << "calling local split coverage " << std::endl;
    i.splitCoverage(m_bamreads);
  }

  if (!m_global_bp.isEmpty()) {
    //std::cout << "calling global split coverage " << std::endl;
    m_global_bp.splitCoverage(m_bamreads);
  }

}

std::ostream& operator<<(std::ostream& out, const AlignedContig &ac) {

  // print the global breakpoint
  if (!ac.m_global_bp.isEmpty())
    out << "Global BP: " << ac.m_global_bp << 
      " ins_aginst_contig " << ac.insertion_against_contig_read_count << 
      " del_against_contig " << ac.deletion_against_contig_read_count  << std::endl;       

  // print the multi-map breakpoints
  for (auto& i : ac.m_local_breaks)
    if (!i.isEmpty())
      out << "Multi-map BP: " << i << " -- " << ac.getContigName() << std::endl;       

  // print the indel breakpoints
  for (auto& i : ac.m_frag_v)
    for (auto& j : i.getIndelBreaks()) 
      if (!j.isEmpty())
	out << "Indel: " << j << " -- " << ac.getContigName() << " ins_a_contig " << ac.insertion_against_contig_read_count << 
	  " del_a_contig " << ac.deletion_against_contig_read_count << std::endl;       

  // print the AlignmentFragments alignments
  for (auto& i : ac.m_frag_v) 
    out << i << " Disc: " << ac.printDiscordantClusters() << " -- " << ac.getContigName() << std::endl;

  // print the break locations for indel deletions
  for (auto& i : ac.m_frag_v) {
    for (auto& j : i.getIndelBreaks()) {
      if (j.isindel && j.insertion == "") // deletion
	out << std::string(j.cpos1, ' ') << "|" << std::string(j.cpos2-j.cpos1-1, ' ') << '|' << "   " << ac.getContigName() << std::endl;	
    }
  }

  // print the contig base-pairs
  out << ac.getSequence() << "    " << ac.getContigName() << std::endl; 
  
  PlottedReadVector plot_vec;

  // print out the individual reads
  for (auto& i : ac.m_bamreads) {
    
    int pos = -1;
    int aln = -1;
    int dum= 0;
    int rc = 0;
    std::string this_cig;
    std::string seq = i.QualityTrimmedSequence(4, dum);
    std::string sr = i.GetZTag("SR");

    // reverse complement if need be
    //int32_t rc = i.GetIntTag("RC");
    //if (rc/* && false*/)
    //  SnowTools::rcomplement(seq);
    
    // get the more complex tags (since there can be multiple annotations per tag)
    std::vector<int> posvec = i.GetSmartIntTag("SL"); // start positions ON CONTIG
    std::vector<int> alnvec = i.GetSmartIntTag("TS"); // start positions ON READ
    std::vector<int> rcvec = i.GetSmartIntTag("RC"); // read reverse complmented relative to contig
    std::vector<std::string> cigvec = i.GetSmartStringTag("SC"); // read against contig CIGAR
    std::vector<std::string> cnvec = i.GetSmartStringTag("CN");

    assert(cnvec.size() == posvec.size());
    size_t kk = 0;
    for (; kk < cnvec.size(); kk++) 
      if (cnvec[kk] == ac.getContigName()) {
	pos = posvec[kk];
	aln = alnvec[kk];
        rc = rcvec[kk];
	this_cig = cigvec[kk];
	break;
      }
    
    // reverse complement if need be
    if (rc)
      SnowTools::rcomplement(seq);      
    
    // trim the sequence if it hangs off the end
    //if (i.GetZTag("SR") == "t147_D0BK6ACXX111110:6:2101:5352:156091") {
    //  std::cerr << "i.PositionEnd() " << i.PositionEnd() << " i.Position " << i.Position() << " len " << ac.getSequence().length() << std::endl;
    //  exit(1);
    //	 }
    if (aln > 0)
      seq = seq.substr(aln, seq.length() - aln);
    
    if ( (pos + seq.length() ) > ac.getSequence().length()) 
      seq = seq.substr(0, ac.getSequence().length() - pos);
    
    
    assert(kk != cnvec.size()); // assure that we found something
    pos = abs(pos);
    int padlen = ac.getSequence().size() - pos - seq.size() + 5;
    padlen = std::max(5, padlen);

     std::stringstream rstream;
     assert(pos < 1e4 && padlen < 1e4); // bug, need to check
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
       prl.contig_len = ac.getSequence().length();
       prl.addRead(&i);
       line_vec.push_back(prl);
     }
   }

   // plot the lines
   for (auto& i : line_vec) 
     out << i << " contig " << ac.getContigName() << std::endl;
   
   return out;
}

void AlignedContig::setMultiMapBreakPairs() {
   
  //if (m_skip)
  //  return;

  // if single mapped contig, nothing to do here
  if (m_frag_v.size() == 1)
    return;

  // initialize the breakpoint, fill with basic info
  BreakPoint bp;
  bp.seq = getSequence();
  bp.num_align = m_frag_v.size();
  assert(bp.num_align > 0);

  //bp.num_frag_v = m_frag_v.size();
  bp.cname = getContigName(); 
  assert(bp.cname.length());
  
  // set the discovar support if it's there
  //parseDiscovarName(bp.disco_norm, bp.disco_tum);

  // walk along the ordered contig list and make the breakpoint pairs  
  for (auto it = m_frag_v.begin(); it != m_frag_v.end() - 1; it++) {
    
    bp.gr1 = SnowTools::GenomicRegion(it->m_align.ChrID(), it->gbreak2, it->gbreak2);
    bp.gr2 = SnowTools::GenomicRegion((it+1)->m_align.ChrID(), (it+1)->gbreak1, (it+1)->gbreak1);
    
    bp.gr1.strand = it->m_align.ReverseFlag() ? '-' : '+'; 
    bp.gr2.strand = (it+1)->m_align.ReverseFlag() ? '+' : '-';

    bp.cpos1 = it->break2; // take the right-most breakpoint as the first
    bp.cpos2 = (it+1)->break1;  // take the left-most of the next one
    
    assert(bp.cpos1 < 10000);
    assert(bp.cpos2 < 10000);

    // set the mapq
    bp.mapq1 = it->m_align.MapQuality();
    bp.mapq2 = (it+1)->m_align.MapQuality();
    
    assert(bp.mapq1 < 1000 && bp.mapq2 < 1000);
    
    // set the local
    bp.local1 = it->local;
    bp.local2 = (it+1)->local;
    
    // set the match length
    //for (CigarOpVec::const_iterator cc = it->align.CigarData.begin(); cc != it->align.CigarData.end(); cc++)
    //  if (cc->Type == 'M')
    //	bp.matchlen1 += cc->Length;
    //for (CigarOpVec::const_iterator cc = (it+1)->align.CigarData.begin(); cc != (it+1)->align.CigarData.end(); cc++)
    //  if (cc->Type == 'M')
    //bp.matchlen2 += cc->Length;
    
    // set the NM
    /*int nmtag;
    if (!it->align.GetTag("NM", nmtag))
      nmtag = 0;
    bp.nm1 = nmtag;
    if (!(it+1)->align.GetTag("NM", nmtag))
      nmtag = 0;
    bp.nm2 = nmtag;
    */

    // set the insertion / homology
    try {
      if (bp.cpos1 >= bp.cpos2)
	bp.homology = m_seq.substr(bp.cpos2, bp.cpos1-bp.cpos2);
      else if (bp.cpos2 >= bp.cpos1)
	bp.insertion = m_seq.substr(bp.cpos1, bp.cpos2-bp.cpos1);
      if (bp.insertion.length() == 0)
	bp.insertion = "";
      if (bp.homology.length() == 0)
	bp.homology = "";
    } catch (...) {
      unsigned hom = abs(bp.cpos1 - bp.cpos2);
      std::cout << "Caught error with contig on fine-getBreakPairs: " << it->m_align.Qname() << std::endl; 
      std::cout << "m_seq length: " << m_seq.length() << " bp.cpos1: " << bp.cpos1 << " bp.cpos2: " << bp.cpos2 << " bp.cpos1-bp.cpos2: " << hom << " m_seq: " << m_seq << std::endl;
    }

    m_local_breaks.push_back(bp);
    
  }
  
  // if this is a double mapping, we are done
  if (m_frag_v.size() == 2) {
    m_global_bp = bp;
    m_local_breaks.clear();
    return;
  }
  
  // go through alignments and find start and end that reach mapq 
  size_t bstart = 1000; //1000 is a dummy
  size_t bend = m_frag_v.size() - 1;
  for (size_t i = 0; i < m_frag_v.size(); i++)
    if (m_frag_v[i].m_align.MapQuality() >= 60) {
      bend = i;
      if (bstart == 1000)
	bstart = i;
    }
  if (bstart == bend || bstart==1000) {
    bstart = 0;
    bend = m_frag_v.size() -1 ;
  }
  assert(bend <= m_frag_v.size());
  assert(bstart <= m_frag_v.size());
  assert(bstart != 1000);
  
  // there are 3+ mappings
  m_global_bp = bp;
  m_global_bp.cpos1 = m_frag_v[bstart].break2; // first mapping
  m_global_bp.gr1.pos1 = m_frag_v[bstart].gbreak2;
  m_global_bp.gr1.pos2 = m_global_bp.gr1.pos1;
  m_global_bp.gr2.pos1 = m_frag_v[bend].gbreak1;
  m_global_bp.gr2.pos2 = m_global_bp.gr2.pos1;
  m_global_bp.gr1.chr = m_frag_v[bstart].m_align.ChrID();
  m_global_bp.gr2.chr = m_frag_v[bend].m_align.ChrID();
  //m_global_bp.pos1  = m_frag_v[bstart].gbreak2;
  m_global_bp.cpos2 = m_frag_v[bend].break1; // last mapping
   //m_global_bp.pos2  = m_frag_v[bend].gbreak1;
   //m_global_bp.refID1 = m_frag_v[bstart].align.RefID; 
   //m_global_bp.refID2 = m_frag_v[bend].align.RefID; 

   // set the strands
   //m_global_bp.gr1.strand = m_frag_v[bstart].align.IsReverseStrand() ? '-' : '+';
   //m_global_bp.gr2.strand = m_frag_v[bend].align.IsReverseStrand()   ? '+' : '-';
  m_global_bp.gr1.strand = !m_frag_v[bstart].m_align.ReverseFlag() ? '+' : '-';
  m_global_bp.gr2.strand = m_frag_v[bend].m_align.ReverseFlag() ? '+' : '-';


   // set the splits
   //m_global_bp.nsplit1 = m_frag_v[bstart].nsplit2;
   //m_global_bp.tsplit1 = m_frag_v[bstart].tsplit2;
   //m_global_bp.nsplit2 = m_frag_v[bend].nsplit1;
   //m_global_bp.tsplit2 = m_frag_v[bend].tsplit1;

   //m_global_bp.nsplit = std::min(m_global_bp.nsplit1, m_global_bp.nsplit2);
   //m_global_bp.tsplit = std::min(m_global_bp.tsplit1, m_global_bp.tsplit2);

   // set the mapping quality
  m_global_bp.mapq1 = m_frag_v[bstart].m_align.MapQuality();
  m_global_bp.mapq2 = m_frag_v[bend].m_align.MapQuality();

   if (m_global_bp.mapq1 > 60 || m_global_bp.mapq2 > 60) {
     std::cerr << "bad mapq GLOBAL" << std::endl;
     std::cerr << m_global_bp.toString() << std::endl;
     std::cerr << " m_frag_v size " << m_frag_v.size() << " bstart " << bstart << " bend " << bend << std::endl;
     exit(EXIT_FAILURE); 
   }

   // set the homologies
   try {
     if (m_global_bp.cpos1 >= m_global_bp.cpos2)
       m_global_bp.homology = m_seq.substr(m_global_bp.cpos2, m_global_bp.cpos1-m_global_bp.cpos2);
     else if (m_global_bp.cpos2 >= m_global_bp.cpos1)
       m_global_bp.insertion = m_seq.substr(m_global_bp.cpos1, m_global_bp.cpos2-m_global_bp.cpos1);
     if (m_global_bp.insertion.length() == 0)
       m_global_bp.insertion = "N";
     if (m_global_bp.homology.length() == 0)
       m_global_bp.homology = "N";
     m_global_bp.order();
   } catch (...) {
     std::cout << "Caught error with contig on global-getBreakPairs: " << getContigName() << std::endl;
   }

}

  void AlignedContig::alignReads(BamReadVector& bav)
  {
    
    // construc the index
    SnowTools::BWAWrapper bw;
    SnowTools::USeqVector v = {   {getContigName(), getSequence()}   };
    bw.constructIndex(v);
    
    // align the reads
    for (auto& i : bav)
      {
	BamReadVector brv;
	int dum = 0;
	std::string seqr = i.QualityTrimmedSequence(4, dum);
	bw.alignSingleSequence(seqr, i.Qname(), brv, false);

	if (brv.size() > 1) {
	  continue;
	  //std::cerr << "Mulitple alignments for " << brv[0].Qname() << std::endl;
	}
	
	if (brv.size() == 0) {
	  continue;
	}
	
	bool length_pass = (brv[0].PositionEnd() - brv[0].Position()) >= (seqr.length() * 0.95);
	bool mapq_pass = brv[0].MapQuality();
	int ins_bases = brv[0].MaxInsertionBases();
	int del_bases = brv[0].MaxDeletionBases();
	
	//std::cerr << "just aligned " << brv[0] << " size " << bav.size() << std::endl;

	// store read2contig alignment info in this read
	if ( length_pass && mapq_pass && ins_bases == 0 && del_bases == 0)
	  {
	    if (brv[0].ReverseFlag())
	      i.SmartAddTag("RC","1");
	    else 
	      i.SmartAddTag("RC","0");
	    i.SmartAddTag("SL", std::to_string(brv[0].Position()));
	    i.SmartAddTag("SE", std::to_string(brv[0].PositionEnd()));
	    i.SmartAddTag("TS", std::to_string(brv[0].AlignmentPosition()));
	    i.SmartAddTag("TE", std::to_string(brv[0].AlignmentEndPosition()));
	    i.SmartAddTag("SC", brv[0].CigarString());
	    i.SmartAddTag("CN", getContigName());
	    
	    m_bamreads.push_back(i);

	    // add the coverage
	    int cc = brv[0].Position();
	    std::string srr = i.GetZTag("SR");
	    if (srr.length()) {
	      if (srr.at(0) == 't')
		while (cc <= brv[0].PositionEnd() && cc < (int)tum_cov.size())
		  ++tum_cov[cc++];
	      else
		while (cc <= brv[0].PositionEnd() && cc < (int)norm_cov.size())
		  ++norm_cov[cc++];
	    }
	    
	  }

	else if (ins_bases && mapq_pass && length_pass)
	  ++insertion_against_contig_read_count;
	else if (del_bases && mapq_pass && length_pass)
	  ++deletion_against_contig_read_count;
      }


    assignSupportCoverage();

  }

  std::pair<int,int> AlignedContig::getCoverageAtPosition(int pos) const {
    if (pos <= 0 || pos >= (int)tum_cov.size())
      return std::pair<int,int>(0,0);
    
    return std::pair<int,int>(tum_cov[pos], norm_cov[pos]);
  }
 
  std::string AlignedContig::printDiscordantClusters() const {
    
    std::stringstream out;
    if (m_dc.size() == 0)
      return "none";
    
    for (std::vector<DiscordantCluster>::const_iterator it = m_dc.begin(); it != m_dc.end(); it++)
      out << *it << " ";
    return out.str();
    
  }/*
     
  // reject if not at last partially in window
  m_skip = true;
  for (auto& ca : m_frag_v)
    if (ca.local)
      m_skip = false;

  // set break pairs for multi mapped reads
  setMultiMapBreakPairs();

  // set the std::min alignment end length
  //if (m_frag_v.size() > 1) {
  //  m_std::min_align_len = std::min(m_frag_v[0].align.Length, m_frag_v.back().align.Length);
  //}

}
*/

  bool AlignedContig::checkLocal(const GenomicRegion& window)
  {
    bool has_loc = false;
    for (auto& i : m_frag_v) 
      if (i.checkLocal(window))
	has_loc = true;
    return has_loc;

  }

  bool AlignmentFragment::checkLocal(const GenomicRegion& window)
  {
    // make a regino for this frag
    GenomicRegion gfrag(m_align.ChrID(), m_align.Position(), m_align.PositionEnd());
    
    if (window.getOverlap(gfrag)) {
      local = true;
      return true;
    }
      
    return false;
  }

  AlignmentFragment::AlignmentFragment(const BamRead &talign, bool flip) {
    
    m_align = talign;
    
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
      if (i.Type != 'M')
	start += i.Length;
      else
	break;
    }
    
    // set the left-right breaks
    unsigned currlen  = 0; 
    
    // cigar is oriented to as is from aligner
    for (auto& i : m_cigar/*m_align.GetCigar()*/ /*align.CigarData*/) { //CigarOpVec::const_iterator j = align.cigar.begin(); j != align.cigar.end(); j++) {
      
      // SET THE CONTIG BREAK (treats deletions and leading S differently)
      // the first M gets the break1, pos on the left
      if (i.Type == 'M' && break1 == -1)
	break1 = currlen;
      if (i.Type != 'D') // m_skip deletions but not leading S, but otherwise update
	currlen += i.Length;
      if (i.Type == 'M') // keeps triggering every M, with pos at the right
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

    assert(break1 < 10000);
    assert(break2 < 10000);

    assert(break1 >= 0);
    assert(break2 >= 0);
    
    // parse right away to see if there are indels on this alignment
    BreakPoint bp;
    size_t fail_safe_count = 0;
    while (parseIndelBreak(bp) && fail_safe_count++ < 100) {
      m_indel_breaks.push_back(bp);
      assert(bp.num_align == 1);
    }
    
    assert(fail_safe_count != 100);
    
    // set the cigar matches
    //indelCigarMatches(nmap, tmap);
    
}

std::ostream& operator<<(std::ostream &out, const AlignmentFragment &c) {
  
  // sets the direction to print
  char jsign = '>'; 
  if (c.m_align.ReverseFlag())
    jsign = '<';
  
  // print the cigar value per base
  for (auto& j : /*c.m_align.GetCigar()*/ c.m_cigar) { //c.align.CigarData) { // print releative to forward strand
    if (j.Type == 'M')
      out << std::string(j.Length, jsign);
    else if (j.Type == 'I') 
      out << std::string(j.Length, 'I');
    else if (j.Type == 'S' || j.Type == 'H')
      out << std::string(j.Length, '.');
  }

  // print contig and genome breaks
  out << "\tC[" << c.break1 << "," << c.break2 << "] G[" << c.gbreak1 << "," << c.gbreak2 << "]";
  
  // add local info
  out << "\tLocal: " << c.local << "\tAligned to: " << (c.m_align.ChrID()+1) << ":" << c.m_align.Position() << "(" << (c.m_align.ReverseFlag() ? "-" : "+") << ") CIG: " << c.m_align.CigarString() << " MAPQ: " << c.m_align.MapQuality();

  // print the del if there is one TODO
  /*  std::string del_string = std::string(c.m_align.Length(), ' ');
  int dpos = 0;
  bool has_del = false;
  for (auto& i : c.m_align.GetCigar()) {
    if (i.Type != 'D')
      dpos += i.Length;
    if (i.Type == 'D' && dpos != 0) {
      del_string[dpos-1] = '|';
      del_string[dpos] = '|';
      has_del = true;
    }
  }
  if (has_del)
    out << del_string << std::endl;
  */

  // print the info
  /*  out << "    " << c.m_align.RefID + 1 << ":" << c.align.Position 
      << " MAPQ: " << c.align.MapQuality << " OrientedCigar: " << BamToolsUtils::cigarToString(c.align.CigarData)
	  << " Start: " << c.start
	  << " Breaks: " << c.break1 << " " << c.break2 
      << " GenomeBreaks: " << c.gbreak1 << " " << c.gbreak2 << " cname " << c.m_name;
  */
    //	  << " Tsplits: " << c.tsplit1 << " " << c.tsplit2 
    //    << " Nsplits: " << c.nsplit1 << " " << c.nsplit2;
  
  return out;
}


  void AlignedContig::checkAgainstCigarMatches(const CigarMap& nmap, const CigarMap& tmap) {
    
    for (auto& i : m_frag_v)
      i.indelCigarMatches(nmap, tmap);
    
  }

void AlignmentFragment::indelCigarMatches(const CigarMap &nmap, const CigarMap &tmap) { 

  // loop through the indel breakpoints
  for (auto& i : m_indel_breaks) {

    assert(i.isindel); // make sure we only call on indels
    if (i.getSpan() <= 0) {
      std::cerr << "weird span detected " << i.getSpan();
      std::cerr << i << std::endl;
      std::cerr << *this << std::endl;
    }

    // get the hash string in same formate as cigar map (eg. pos_3D)
    std::string st = i.getHashString();

    // check if this breakpoint from assembly is in the cigarmap
    CigarMap::const_iterator ffn = nmap.find(st);
    CigarMap::const_iterator fft = tmap.find(st);

    // if it is, add it
    if (ffn != nmap.end())
      i.ncigar = ffn->second;
    if (fft != tmap.end())
      i.tcigar = fft->second;

  }
}

bool AlignedContig::isWorse(const AlignedContig &ac) const {

  // TODO.  right now is TRUE if any more than 1 bp per contig
  std::vector<const BreakPoint*> bpv1 = getAllBreakPointPointers();
  std::vector<const BreakPoint*> bpv2 = ac.getAllBreakPointPointers();
  if (bpv1.size() != 1 || bpv2.size() != 1)
    return false;
  
  bool same = (*(bpv2.back())) == (*(bpv1.back()));
  if (!same)
    return false;
    
  int mapq_this = std::max(bpv1.back()->mapq1, bpv1.back()->mapq2);
  int mapq_that = std::max(bpv2.back()->mapq1, bpv2.back()->mapq2);

  if (mapq_this < mapq_that) {
    return true;
  } 

  if ( (m_seq.length() < ac.m_seq.length()) && (mapq_this == mapq_that)) {
    return true;
  }

  return false;
  

}

bool AlignmentFragment::parseIndelBreak(BreakPoint &bp) {

  // make sure we have a non-zero cigar
  if (m_cigar.size() == 0) {
    std::cerr << "CIGAR of length 0 on " << *this << std::endl;
    return false;
  }

  // reject if it has small matches, could get confused. Fix later
  for (auto& i : m_cigar) 
    if (i.Type == 'M' && i.Length < 5)
      return false;

  // reject if first alignment is I or D
  if (m_cigar[0].Type == 'I' || m_cigar[0].Type == 'D' || m_cigar[m_cigar.size()-1].Type == 'D' || m_cigar[m_cigar.size()-1].Type == 'I') {
    std::cerr << "rejcting cigar for starting in I or D" << std::endl;
    return false;
  }

  // use next available largest D / I
  size_t loc = 0; // keep track of which cigar field
  for (auto& i : m_cigar) {
    loc++;
    if (i.Type == 'D' || i.Type == 'I') {
      
      // if it starts / stop with I, D, reject
      if (loc == 1 || loc == m_cigar.size()) 
	return false;
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
  bp.isindel = true;

  // set the number of alignments for this bp
  bp.num_align = 1;

  //bp.insertion = "";
  //bp.homology = "";

  int curr = 0;
  int gcurrlen = -1;

  bp.gr1.pos1 = -1;
  bp.gr1.pos2 = -1;
  bp.gr2.pos1 = -1;
  bp.gr2.pos2 = -1;
  bp.gr1.chr = m_align.ChrID();
  bp.gr2.chr = m_align.ChrID();

  bp.cpos1 = -1;
  bp.cpos2 = -1;

  bp.cname = m_align.Qname();
  bp.mapq1 = m_align.MapQuality();
  bp.mapq2 = m_align.MapQuality();
  bp.gr1.strand = '+';
  bp.gr2.strand = '-';

  assert(bp.cname.length());

  size_t count = 0; // count to make sure we are reporting the right indel

  for (auto& i : m_cigar) { // want breaks in CONTIG coordaintes, so use oriented cigar
    count++;

    // set the contig breakpoint
    if (i.Type == 'M' || i.Type == 'I') 
      curr += i.Length;
    if (i.Type == 'D' && bp.cpos1 == -1 && count == idx) {

      bp.cpos1 = curr-1;
      bp.cpos2 = curr;
    } 
    if (i.Type == 'I' && bp.cpos1 == -1 && count == idx) {
      bp.cpos1 = curr - i.Length - 1;
      bp.cpos2 = curr - 1;
      bp.insertion = m_align.Sequence().substr(curr, i.Length);
    }

    // set the genome breakpoint
    if (bp.cpos1 > 0) {
      if (i.Type == 'D') {
	if (!m_align.ReverseFlag()) {
	  bp.gr1.pos1 =  m_align.Position() + gcurrlen; // dont count this one//bp.cpos1 + align.Position; //gcurrlen + align.Position;
	  bp.gr2.pos1 =  bp.gr1.pos1 + i.Length + 1;
	} else {
	  bp.gr2.pos1 =  (m_align.PositionEnd()-1) - gcurrlen; //bp.cpos1 + align.Position; //gcurrlen + align.Position;
	  bp.gr1.pos1 =  bp.gr2.pos1 - i.Length - 1;
	}
      } else if (i.Type == 'I') {
	if (!m_align.ReverseFlag()) {
	  bp.gr1.pos1 = m_align.Position() + gcurrlen; //gcurrlen + align.Position;
	  bp.gr2.pos1 = bp.gr1.pos1 + 1;	
	} else {
	  // GetEndPosition is 1 too high
	  bp.gr2.pos1 = (m_align.PositionEnd()-1) - gcurrlen; //gcurrlen + align.Position;
	  bp.gr1.pos1 = bp.gr2.pos1 - 1;	
	}
      }
      break; // already got it, so quit cigar loop
    }
    
    // update the position on the genome
    if (i.Type == 'M' || i.Type == 'D') {
      gcurrlen += i.Length;
    } 


  } // end cigar loop

  // set the dummy other end
  bp.gr1.pos2 = bp.gr1.pos1;
  bp.gr2.pos2 = bp.gr2.pos1;

  // error check the length
  //int seq_len = m_seq.length();
  /*  if (bp.cpos1 > seq_len || bp.cpos2 > seq_len) {
    cerr << "bp " << bp << std::endl;
    cerr << "CA " << *this << std::endl;
    cerr << "seq len " << seq_len << std::endl;
    }*/
  //assert(bp.cpos1 <= seq_len);
  //assert(bp.cpos2 <= seq_len);

  bp.order();
  
  return true;
}

  //bool AlignedContig::parseDiscovarName(size_t &tumor, size_t &normal) {
  /*
  std::string s = "";
  bool valid = m_frag_v[0].m_align.GetZTag("TN", s);
  if (!valid)
    return false;

  // set the tumor support
  std::regex reg("^t([0-9]+)_*");
  std::smatch match;
  //if (std::regex_search(s.begin(), s.end(), match, reg))
  if (std::regex_search(s, match, reg))
    tumor = std::stoi(match[1].str());
  
  // set the normal support
  std::regex reg2("^t[0-9]+n([0-9]+)");
  std::smatch match2;
  if (std::regex_search(s, match2, reg2))
    normal = std::stoi(match2[1].str());
  
  return false;
  */
  //}

std::vector<const BreakPoint*> AlignedContig::getAllBreakPointPointers() const  {

  std::vector<const BreakPoint*> out;
  for (auto& i : m_frag_v) {
    for (auto& k : i.m_indel_breaks)
      out.push_back(&k);
  }
  
  if (!m_global_bp.isEmpty())
    out.push_back(&m_global_bp);
  
  return out;
}

std::vector<BreakPoint> AlignedContig::getAllBreakPoints() const {

  std::vector<BreakPoint> out;
  for (auto& i : m_frag_v) {
    for (auto& k : i.m_indel_breaks)
      out.push_back(k);
  }
  
  if (!m_global_bp.isEmpty())
    out.push_back(m_global_bp);
  
  return out;
}


bool AlignedContig::hasVariant() const { 
  
  if (!m_global_bp.isEmpty())
    return true;

  for (auto& i : m_frag_v)
    if (i.local && i.m_indel_breaks.size())
      return true;

  return false;
  
}

  void AlignedContig::addDiscordantCluster(DiscordantClusterMap& dmap)
  {
    
    // loop through the breaks and compare with the map
    for (auto& i : m_local_breaks)
      i.__combine_with_discordant_cluster(dmap);
    m_global_bp.__combine_with_discordant_cluster(dmap);

  }

  void AlignedContig::assignSupportCoverage() {
    
    // go through the indel breaks and assign support coverage
    for (auto& j : m_frag_v) {
      for (auto& i : j.m_indel_breaks) {
	std::pair<int,int> p1 = getCoverageAtPosition(i.cpos1);
	std::pair<int,int> p2 = getCoverageAtPosition(i.cpos2);
	i.tcov_support = std::min(p1.first, p2.first);
	i.ncov_support = std::min(p1.second, p2.second);
      }
    }
    
    // go through the SV breaks and assign support coverage
    for (auto& i : m_local_breaks) {
      std::pair<int,int> p1 = getCoverageAtPosition(i.cpos1);
      std::pair<int,int> p2 = getCoverageAtPosition(i.cpos2);
      i.tcov_support = std::min(p1.first, p2.first);
      i.ncov_support = std::min(p1.second, p2.second);
    }
    
    std::pair<int,int> p1 = getCoverageAtPosition(m_global_bp.cpos1);
    std::pair<int,int> p2 = getCoverageAtPosition(m_global_bp.cpos2);
    m_global_bp.tcov_support = std::min(p1.first, p2.first);
    m_global_bp.ncov_support = std::min(p1.second, p2.second);

  }
  
  
}
