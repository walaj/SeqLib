#include "SnowTools/BreakPoint.h"

#include <getopt.h>
#include "SnowTools/gzstream.h"

#define SPLIT_BUFF 10

// define repeats
static std::vector<std::string> repr = {"AAAAAAAA", "TTTTTTTT", "CCCCCCCC", "GGGGGGGG", 
				 "TATATATATATATATA", "ATATATATATATATAT", 
				 "GCGCGCGCGCGCGCGC", "CGCGCGCGCGCGCGCG", 
				 "TGTGTGTGTGTGTGTG", "GTGTGTGTGTGTGTGT", 
				 "TCTCTCTCTCTCTCTC", "CTCTCTCTCTCTCTCT", 
				 "CACACACACACACACA", "ACACACACACACACAC", 
				 "GAGAGAGAGAGAGAGA", "AGAGAGAGAGAGAGAG"};

namespace SnowTools {

  // send breakpoint to a string
  std::ostream& operator<<(std::ostream& out, const BreakPoint& b) {
    
    out << b.b1.chr_name << ":" << b.b1.gr.pos1 << "(" << b.b1.gr.strand << ")" << "-" 
	<< b.b2.chr_name << ":" << b.b2.gr.pos1 << "(" << b.b2.gr.strand << ")" <<
      " SPAN: " << b.getSpan() << " MAPQ: " << 
      b.b1.mapq << "/" << b.b2.mapq << " HOM: " << 
      b.homology << " INS: " << b.insertion << " NS: " << 
      b.nsplit << " TS: " << b.tsplit << " TD: " << b.dc.tcount << " ND: " << b.dc.ncount 
	<< " NC " << b.ncigar << " TC " << b.tcigar << " TCOV " << b.tcov << " NCOV " << b.ncov << " -- " << b.cname; 
    return out;
  }

  // make the file string
  std::string BreakPoint::toFileString(bool noreads) {
    
    // make sure we already ran scoring
    assert(evidence.length());
    assert(confidence.length());
    
    std::string sep = "\t";
    std::stringstream ss;
    
    // put the read names into a string
    if (!noreads) 
      __format_readname_string();
    
    ss << b1.chr_name << sep << b1.gr.pos1 << sep << b1.gr.strand << sep 
       << b2.chr_name << sep << b2.gr.pos1 << sep << b2.gr.strand << sep 
       << ref << sep << alt << sep 
       << getSpan() << sep
       << b1.mapq << sep << b2.mapq << sep 
       << nsplit << sep << tsplit << sep
       << b1.sub_n << sep << b2.sub_n << sep
       << dc.ncount << sep << dc.tcount << sep
       << dc.mapq1 << sep << dc.mapq2 << sep
       << ncigar << sep << tcigar << sep
       << (homology.length() ? homology : "x") << sep 
       << (insertion.length() ? insertion : "x") << sep 
       << cname << sep
       << num_align << sep 
       << confidence << sep << evidence << sep
       << quality << sep
       << secondary << sep << somatic_score << sep 
       << pon << sep << (repeat_seq.length() ? repeat_seq : "x") << sep 
       << ncov << sep << tcov << sep << af_n << sep << af_t << sep
       << blacklist << sep << (rs.length() ? rs : "x") << sep 
       << (read_names.length() ? read_names : "x");
    
    return ss.str();
    
  }
  
  BreakEnd::BreakEnd(const GenomicRegion& g, int mq, const std::string& chr_n) {
    gr = g;
    mapq = mq;
    chr_name = chr_n;
    cpos = -1; nm = -1; matchlen = -1; tsplit = 0; nsplit = 0; sub_n = 0; local = false;
  }
  
  ReducedBreakEnd::ReducedBreakEnd(const GenomicRegion& g, int mq, const std::string& chr_n) {
    gr = g;
    mapq = mq;
    chr_name = chr_n;
    //chr_name = (char*)malloc(chr_n.length() + 1);
    //strcpy(chr_name, chr_n.c_str());
  }


  // make a breakpoint from a discordant cluster 
  BreakPoint::BreakPoint(const DiscordantCluster& tdc, const BWAWrapper * bwa) {
    
    num_align = 0;
    dc = tdc;
    
    assert(tdc.reads.size());
    std::string chr_name = bwa->ChrIDToName(tdc.reads.begin()->second.ChrID());

    int pos1 = (dc.m_reg1.strand == '+') ? dc.m_reg1.pos2 : dc.m_reg1.pos1;
    int pos2 = (dc.m_reg2.strand == '+') ? dc.m_reg2.pos2 : dc.m_reg2.pos1;
    b1 = BreakEnd(GenomicRegion(dc.m_reg1.chr, pos1, pos1), dc.mapq1, chr_name);
    b2 = BreakEnd(GenomicRegion(dc.m_reg2.chr, pos2, pos2), dc.mapq2, chr_name);
    b1.gr.strand = dc.m_reg1.strand;
    b2.gr.strand = dc.m_reg2.strand;

    cname = dc.toRegionString();
    
    //assert(b1.gr < b2.gr);

  }
  
  bool BreakPoint::hasDiscordant() const {
    return (dc.ncount || dc.tcount);
  }
  
  /** 
   * Has at least two supporting reads
   */
  bool BreakPoint::hasMinimal() const {
    if (tsplit + nsplit + dc.tcount + dc.ncount >= 2)
      return true;
    return false;
  }
  
  bool BreakPoint::operator==(const BreakPoint &bp) const {
    return (b1.gr == bp.b1.gr && b2.gr == bp.b2.gr); 
  }
    
  /*void BreakPoint::repeatFilter(faidx_t * f) {

    if (!f) {
      std::cerr << "Need to open reference for reading. BreakPoint::repeatFilter" << std::endl;
      return;
    }
   
    if (gr1.chr >= 24)
      return;


    GenomicRegion gr = gr1;
    gr.pad(40);

    int len;    
    std::string chrstring = GenomicRegion::chrToString(gr.chr);
    char * seq = faidx_fetch_seq(f, const_cast<char*>(chrstring.c_str()), gr.pos1-1, gr.pos2-1, &len);
    std::string seqr = std::string(seq);

    for (auto& i : repr)
      if (seqr.find(i) != std::string::npos)
	repeat_seq = i;
    
	}*/

  BreakPoint::BreakPoint(const std::string &line, bam_hdr_t* h) {
    
    if (!h) {
      std::cerr << "BreakPoint::BreakPoint - Must supply non-empty header" << std::endl;
      exit(EXIT_FAILURE);
    }

    std::istringstream iss(line);
    std::string val;
    size_t count = 0;

    std::string chr1, pos1, chr2, pos2, chr_name1, chr_name2; 
    char strand1 = '*', strand2 = '*';
    while (std::getline(iss, val, '\t')) {

      try{
	switch(++count) {
	case 1: chr1 = val; chr_name1 = val; break;
	case 2: pos1 = val; break; 
	case 3: assert(val.length()); strand1 = val.at(0); break;
	case 4: chr2 = val; chr_name2 = val; break;
	case 5: pos2 = val; break; 
	case 6: assert(val.length()); strand2 = val.at(0); break;
	case 7: ref = val; break;
	case 8: alt = val; break;
	  //case 9: span = stoi(val); break; // automatically calculated
	case 10: 
	  b1 = BreakEnd(GenomicRegion(chr1, pos1, pos1, h), std::stoi(val), chr_name1); b1.gr.strand = strand1; break;
	case 11:
	  b2 = BreakEnd(GenomicRegion(chr2, pos2, pos2, h), std::stoi(val), chr_name2); b2.gr.strand = strand2; break;
	case 12: nsplit = std::stoi(val); break;
	case 13: tsplit = std::stoi(val); break;
	case 14: b1.sub_n = std::stoi(val); break;
	case 15: b2.sub_n = std::stoi(val); break;
	case 16: dc.ncount = std::stoi(val); break;
	case 17: dc.tcount = std::stoi(val); break;
	case 18: dc.mapq1 = std::stoi(val); break;  
	case 19: dc.mapq2 = std::stoi(val); break;  
	case 20: ncigar = std::stoi(val); break;
	case 21: tcigar = std::stoi(val); break;
	case 22: homology = (val == "x" ? "" : val); break;
	case 23: insertion = (val == "x" ? "" : val); break;
	case 24: cname = val; break;
	case 25: num_align = std::stoi(val); break;
	case 26: confidence = val; break;
	case 27: evidence = val; break;
	case 28: quality = std::stoi(val); break;
	case 29: secondary = val == "1";
	case 30: somatic_score = std::stod(val); break;
	case 31: pon = std::stoi(val); break;
	case 32: repeat_seq = val; break;
	case 33: ncov = std::stoi(val); break;
	case 34: tcov = std::stoi(val); break;
	case 35: af_n = stod(val); break;
	case 36: af_t = stod(val); break;
	case 37: blacklist = (val=="1"); break;
	case 38: rs = val; break;
	case 39: read_names = val; break;
	}
      } catch(...) {
	std::cerr << "caught stoi error on: " << val << std::endl;
	std::cerr << line << std::endl;
	exit(1);
      }
    }
  }
  
  void BreakPoint::splitCoverage(BamReadVector &bav) {
    
    // zero the counts
    tsplit = 0; nsplit = 0;
    b1.tsplit = 0; b1.nsplit = 0; b2.tsplit = 0; b2.nsplit = 0;
    
    int rightbreak1= b1.cpos + SPLIT_BUFF; // read must extend this far right of break1
    int leftbreak1 = b1.cpos - SPLIT_BUFF; // read must extend this far left break2
    int rightbreak2= b2.cpos + SPLIT_BUFF;
    int leftbreak2 = b2.cpos - SPLIT_BUFF;

    for (auto& j : bav) {
      
      // for indels, reads with ins alignments to del contig dont count, vice versa
      bool read_should_be_skipped = false;
      if (num_align == 1) {
	std::vector<std::string> cigvec = j.GetSmartStringTag("SC"); // read against contig CIGAR
	std::vector<std::string> cnvec = j.GetSmartStringTag("CN");
	assert(cigvec.size() == cnvec.size());
	size_t kk = 0;
	for (; kk < cnvec.size(); kk++) 
	  if (cnvec[kk] == cname) 
	    if ( (cigvec[kk].find("D") != std::string::npos/* && insertion.length()*/) || 
	    	 (cigvec[kk].find("I") != std::string::npos/* && !insertion.length()*/)) {
	      read_should_be_skipped = true;
	      break;
	    }
      }
      
      if (read_should_be_skipped)
	continue;
      
      // get read ID
      std::string sr = j.GetZTag("SR");
      assert(sr.at(0) == 't' || sr.at(0) == 'n');
      bool tumor_read = (sr.at(0) == 't');

      // get the contig qname
      std::string qname = j.GetSmartStringTag("CN").back(); //GetStringTag(j, "CN").back();
      int pos = j.GetSmartIntTag("SL").back();
      int te = 0;
      try {
	te = j.GetSmartIntTag("SE").back();
      } catch (...) {
	std::cerr << "error grabbing SE tag for tag " << j.GetZTag("SE") << std::endl;
      }
      assert(qname == cname);
    
      int rightend = te; //seq.length();
      int leftend  = pos;
      bool issplit1 = (leftend <= leftbreak1) && (rightend >= rightbreak1);
      bool issplit2 = (leftend <= leftbreak2) && (rightend >= rightbreak2);
      
      // add the split reads for each end of the break
      if (issplit1 || issplit2)
	reads.push_back(j);
      
      // update the counters for each break end
      if (issplit1 && tumor_read)
	++b1.tsplit;
      if (issplit1 && !tumor_read)
	++b1.nsplit;
      if (issplit2 && tumor_read)
	++b2.tsplit;
      if (issplit2 && !tumor_read)
	++b2.nsplit;
      
      // read spans both ends
      if ((issplit1 || issplit2) && tumor_read)
	++tsplit;
      if ((issplit1 || issplit2) && !tumor_read)
	++nsplit;
      
    }
  }
  
  std::string BreakPoint::getHashString() const {
    
    bool isdel = insertion.length() == 0;
    //if (isdel) // del breaks are stored as last non-deleted base. CigarMap stores as THE deleted base
    //  pos1++;
    std::string st = std::to_string(b1.gr.chr) + "_" + std::to_string(b1.gr.pos1) + "_" + std::to_string(this->getSpan()) + (isdel ? "D" : "I");
    return st;
  }
  
  int BreakPoint::checkPon(const PONFilter * p) {
    
    // only built for indels for now
    if (num_align != 1)
      return 0;
    
    std::string chr = std::to_string(b1.gr.chr);
    
    std::string key1, key2, key3, key4, key5;
    int span = getSpan();
    std::string type = (insertion != "") ? "I" : "D";

    key1 = chr + "_" + std::to_string(b1.gr.pos1-1);// + "_" + std::to_string(span) + type;
    key2 = chr + "_" + std::to_string(b1.gr.pos1  );// + "_" + std::to_string(span) + type;
    key3 = chr + "_" + std::to_string(b1.gr.pos1+1);// + "_" + std::to_string(span) + type;
    key4 = chr + "_" + std::to_string(b1.gr.pos1-2);// + "_" + std::to_string(span) + type;
    key5 = chr + "_" + std::to_string(b1.gr.pos1+2);//+ "_" + std::to_string(span) + type;
    
    pon = std::max(p->NSamps(key1), pon);
    pon = std::max(p->NSamps(key2), pon);
    pon = std::max(p->NSamps(key3), pon);
    pon = std::max(p->NSamps(key4), pon);
    pon = std::max(p->NSamps(key5), pon);
    
    return pon;
    
  }
  
  void BreakPoint::addAllelicFraction(STCoverage * t_cov, STCoverage * n_cov) {
    tcov = t_cov->getCoverageAtPosition(b1.gr.chr, b1.gr.pos1); 
    ncov = n_cov->getCoverageAtPosition(b1.gr.chr, b1.gr.pos1); 
  }
  
  std::string BreakPoint::toPrintString() const {
    
    std::stringstream ss;
    
    std::string s_af_t = std::to_string(af_t);
    if (s_af_t.length() > 4)
      s_af_t = s_af_t.substr(0, 4);      
    std::string s_af_n = std::to_string(af_n);
    if (s_af_n.length() > 4)
      s_af_n = s_af_n.substr(0, 4);
    
    size_t rs_t = std::count(rs.begin(), rs.end(), 'r');
    
    if (isindel)
      ss << ">" << (insertion.size() ? "INS: " : "DEL: ") << getSpan() << " " << b1.gr /*<< " "  << cname */
	 << " T/N split: " << tsplit << "/" << nsplit << " T/N cigar: " 
	 << tcigar << "/" << ncigar << " T/N AF " << s_af_t << "/" << s_af_n << " T/N Cov " << tcov << "/" << ncov << " DBSNP: " << rs_t;
    else
      ss << "> SV: " << b1.gr.pointString() << " to " << b2.gr.pointString() << " SPAN " << getSpan() /*<< " "  << cname */
	 << " T/N split: " << tsplit << "/" << nsplit << " T/N disc: " 
	 << dc.tcount << "/" << dc.ncount << " " << evidence;
    
    return ss.str();
    
  }
  
  double BreakPoint::__sv_is_somatic() const {

    double somatic_ratio = 100;
    size_t ncount = std::max(nsplit, dc.ncount);
    if (ncount > 0)
      somatic_ratio = (std::max(tsplit,dc.tcount)) / ncount;
    
    if (somatic_ratio >= 12 && ncount < 2) 
      return somatic_ratio;
    return 0;
  }

  double BreakPoint::__indel_is_somatic() const {
    
    double somatic_ratio = 100;
    size_t ncount = std::max(nsplit, ncigar);
    if (ncount > 0)
      somatic_ratio = (std::max(tsplit, tcigar)) / ncount;

    // delete if crosses DBSnp
    if (!rs.empty()) 
      return 0;
    if ( (pon * (1-af_t)) >= 2 && tcov >= 5)
      return 0;
    if (pon > 1 && tcov < 5)
      return 0;

    if (somatic_ratio >= 10 && ncount < 2 && af_n < 0.05) 
      return somatic_ratio;

    return 0;

  }

  void BreakPoint::checkBlacklist(GRC &grv) {
    
    // only check for indels
    if (num_align != 1)
      return;
    
    if (grv.findOverlapping(b1.gr))
      blacklist = true;
  }
  
  void BreakPoint::__set_homologies_insertions() {
    try { 
      if (b1.cpos > b2.cpos)
	homology = seq.substr(b2.cpos, b1.cpos-b2.cpos);
      else if (b2.cpos > b1.cpos)
	insertion = seq.substr(b1.cpos, b2.cpos-b1.cpos);
      if (insertion.length() == 0)
	insertion = "x";
      if (homology.length() == 0)
	homology = "x";
    } catch (...) {
      std::cerr << "Caught error with contig on global-getBreakPairs: " << cname << std::endl;
      std::cerr << b1.cpos << " " << b2.cpos << " seq.length() " << seq.length() << " num_align " << num_align << std::endl;
    }
  }
  
  BreakEnd::BreakEnd(const BamRead& b) {
    gr.chr = b.ChrID(); 
    gr.pos1 = -1;
    gr.pos2 = -1;
    cpos = -1;
    mapq = b.MapQuality();
    chr_name = b.GetZTag("MC"); 
    assert(chr_name.length());
    assert(chr_name != "23");
  }

  void BreakPoint::__combine_with_discordant_cluster(DiscordantClusterMap& dmap)
  {
    int PAD = 400;
    GenomicRegion bp1 = b1.gr;
    GenomicRegion bp2 = b2.gr;
    bp1.pad(PAD);
    bp2.pad(PAD);

    for (auto& d : dmap)
      {
	bool bp1reg1 = bp1.getOverlap(d.second.m_reg1) > 0;
	bool bp2reg2 = bp2.getOverlap(d.second.m_reg2) > 0;

	bool pass = bp1reg1 && bp2reg2;

	/*	std::cerr << " gr1 " << gr1 << " gr2 " << gr2 << std::endl << 
	  " m_reg1 " << d.second.m_reg1 << " m_reg2 " << 
	  d.second.m_reg2 << std::endl << " bp1 " << bp1 << 
	  " bp2 " << bp2 << " pass " << pass << 
	  " bp1reg1 " << bp1reg1 << " bp2reg2 " << bp2reg2 << std::endl << 
	  " bp1reg2 " << bp1reg2 << " bp2reg1 " << bp2reg1 << std::endl;
	*/

	if (pass)
	  // check that we haven't already added a cluster to this breakpoint
	  // if so, chose the one with more tumor support
	  if (dc.isEmpty() || dc.tcount < d.second.tcount) {
	      dc = d.second;
	      d.second.m_contig = cname;
	  } 
	
      }
    
  }

  void BreakPoint::__set_evidence() {

    bool isdisc = (dc.tcount + dc.ncount) != 0;
    bool issplit = (tsplit + nsplit) != 0;//dc.m_contig != ""; 

    if (num_align == 1)
      evidence = "INDEL";
    else if ( isdisc && issplit )
      evidence = "ASDIS";
    else if ( isdisc )
      evidence = "DSCRD";
    else if (num_align == 2)
      evidence = "ASSMB";
    else if (num_align > 2)
      evidence = "COMPL";
    else 
      std::cerr << "num_align" << num_align << " isdisc " << " issplit "  << std::endl;

    assert(evidence.length());
  }

  void BreakPoint::__set_allelic_fraction() {

    if (tcov > 0) 
      af_t = static_cast<double>(std::max(tsplit, tcigar)) / static_cast<double>(tcov);
    else
      af_t = 0;
    if ( ncov > 0) 
      af_n = static_cast<double>(std::max(nsplit, ncigar)) / static_cast<double>(ncov);
    else
      af_n = 0;
  }

  void BreakPoint::__score_assembly_only() {

    int split_count = tsplit + nsplit;
    //int this_mapq1 = b1.local ? 60 : b1.mapq;
    //int this_mapq2 = b2.local ? 60 : b2.mapq;
    int span = getSpan();
    //bool germ = dc.ncount > 0 || nsplit > 0;

    if (seq.length() < 101 + 30)
      confidence = "TOOSHORT";
    else if (split_count < 6 && (span > 1500 || span == -1))  // large and inter chrom need 7+
      confidence = "NODISC";
    else if (std::max(b1.mapq, b2.mapq) <= 50 || std::min(b1.mapq, b2.mapq) <= 10) 
      confidence = "LOWMAPQ";
    else if ( split_count <= 3 && (span <= 1500 && span != -1) ) // small with little split
      confidence = "WEAKASSEMBLY";
    //else if ( (germ && span == -1) || (germ && span > 1000000) ) // super short alignemtns are not to be trusted. Also big germline events
    //  confidence = "WEAKASSEMBLY";
    else if ((b1.sub_n && b1.mapq < 30) || (b2.sub_n && b2.mapq < 30))
      confidence = "MULTIMATCH";
    else if (secondary)
      confidence = "SECONDARY";
    else if ( (insertion.length() >= 100 || span < 1000) && (tsplit <= 8 || nsplit) ) // be extra strict for huge insertions or really low spans
      confidence = "WEAKASSEMBLY";
    else
      confidence = "PASS";
    
  }

  void BreakPoint::__score_assembly_dscrd() {

    int this_mapq1 = /*b1.local ? 60 : */b1.mapq;
    int this_mapq2 = /*b2.local ? 60 : */b2.mapq;
    int span = getSpan();
    bool germ = dc.ncount > 0 || nsplit > 0;

    int max_a_mapq = std::max(this_mapq1, dc.mapq1);
    int max_b_mapq = std::max(this_mapq2, dc.mapq2);

    int min_assm_mapq = std::min(this_mapq1, this_mapq2);
    int max_assm_mapq = std::max(this_mapq1, this_mapq2);
    
    int total_count = nsplit + tsplit + dc.ncount + dc.tcount;

    double min_disc_mapq = std::min(dc.mapq1, dc.mapq2);

    if (max_a_mapq <= 10 || max_b_mapq <= 10)
      confidence = "LOWMAPQ";
    else if (std::max(max_a_mapq, max_b_mapq) <= 30)
      confidence = "LOWMAPQ";
    //if ( (min_disc_mapq <= 10 && min_assm_mapq < 30) || (max_assm_mapq < 40))
    //  confidence = "LOWMAPQ";
    else if ( std::max(tsplit, nsplit) == 0 || total_count < 4 || (germ && (total_count <= 6) )) // stricter about germline
      confidence = "WEAKASSEMBLY";
    else if ( total_count < 15 && germ && span == -1) // be super strict about germline interchrom
      confidence = "WEAKASSEMBLY";
    else if ((b1.sub_n && dc.mapq1 < 1) || (b2.sub_n && dc.mapq2 < 1))
      confidence = "MULTIMATCH";
    else if (secondary)
      confidence = "SECONDARY";
    else
      confidence = "PASS";
  }

  void BreakPoint::__score_dscrd() {
    int disc_count = dc.ncount + dc.tcount;
    if (std::min(dc.mapq1, dc.mapq2) < 20 || std::max(dc.mapq1, dc.mapq2) <= 30) // mapq here is READ mapq (37 std::max)
      confidence = "LOWMAPQ";
    else if (std::min(dc.mapq1, dc.mapq2) == 37 && disc_count >= 6)
      confidence = "PASS";
    else if ( disc_count < 8 || (dc.ncount > 0 && disc_count < 12) )  // be stricter about germline disc only
      confidence = "WEAKDISC";
    else 
      confidence = "PASS";
  }

  void BreakPoint::__score_indel() {

    assert(b1.mapq == b2.mapq);
    
    __set_allelic_fraction();

    double max_af = std::max(af_t, af_n);
    int cigar_count = ncigar+tcigar; 
    int max_count = std::max(tsplit+nsplit, cigar_count);
    bool blacklist_and_low_count = blacklist && (tsplit + nsplit) < 5 && (tcigar + ncigar) < 5;
    bool blacklist_and_low_AF = (max_af < 0.2 && max_count < 8) && blacklist;

    if (rs.length())
      confidence="DBSNP";
    if (blacklist && pon > 3)
      confidence="GRAYLISTANDPON";
    else if (blacklist_and_low_count || blacklist_and_low_AF)
      confidence="LOWAF";
    else if ( (max_count < 4 && max_af < 0.2) || (max_count < 2 && b1.mapq < 60) || (max_count < 5 && b1.mapq < 30))
      confidence="WEAKASSEMBLY";
    else if (b1.mapq < 10)
      confidence="LOWMAPQ";
    else if ( (max_af < 0.2 && (nsplit+ncigar) > 0) || (max_af < 0.30 && nsplit == 0 && tsplit < 3)) // more strict for germline bc purity is not issue
      confidence = "LOWAF";
    else if (ncov <= 5)
      confidence = "LOWNORMCOV";
    else
      confidence="PASS";

  }

  void BreakPoint::scoreBreakpoint() {
    
    // set the evidence (INDEL, DSCRD, etc)
    __set_evidence();
    
    // ensure that discordant cluster is oriented
    //if (hasDiscordant()) 
    //  assert(dc.m_reg1 < dc.m_reg2);
    
    // ensure that breakpoint is oriented
    assert(valid()); 
    
    // do the scoring
    if (evidence == "ASSMB" || (evidence == "COMPL" && (dc.ncount + dc.tcount)==0))
      __score_assembly_only();
    else if (evidence == "ASDIS" || (evidence == "COMPL" && (dc.ncount + dc.tcount)))
      __score_assembly_dscrd();
    else if (evidence == "DSCRD")
      __score_dscrd();
    else if (evidence == "INDEL") 
      __score_indel();
    else {
      std::cerr << "evidence " << evidence << std::endl;
      std::cerr << "BreakPoint not classified. Exiting" << std::endl;
      exit(EXIT_FAILURE);
    }

    somatic_score = (evidence == "INDEL") ? __indel_is_somatic() : __sv_is_somatic();

    if (confidence == "PASS")
      quality = 99;
    else
      quality = 0;
    
    assert(confidence.length());
    assert(getSpan() > -2);

  }

  void BreakPoint::order() {

    if (b1.gr < b2.gr)
      return;
    
    flip(b1, b2);
    
  }

  std::string BreakPoint::__format_readname_string() {
    
    std::string supporting_reads = "";
    std::unordered_map<std::string, bool> supp_reads;
    
    //add the discordant reads
    for (auto& r : dc.reads) {
      std::string tmp = r.second.GetZTag("SR");
      supp_reads[tmp] = true;
    }
    for (auto& r : dc.mates) {
      std::string tmp = r.second.GetZTag("SR");
      supp_reads[tmp] = true;
    }
    
    //add the reads from the breakpoint
    for (auto& r : reads) {
      std::string tmp = r.GetZTag("SR");
      supp_reads[tmp];
    }
    
    // print reads to a string, delimit with a ,
    size_t lim = 0;
    for (auto& i : supp_reads) {
      if (++lim > 50)
	break;
      supporting_reads = supporting_reads + "," + i.first;
    }
    if (supporting_reads.size() > 0)
      supporting_reads = supporting_reads.substr(1, supporting_reads.size() - 1); // remove first _
    
    if (read_names.length() == 0)
      read_names = supporting_reads;
    
    return read_names;
  }

  bool BreakPoint::valid() const {

    // debug
    return true;
    
    if (!(b1.gr.strand == '+' || b1.gr.strand == '-') || !(b2.gr.strand == '+' || b2.gr.strand == '-')) {
      std::cerr << "b1.strand " << b1.gr.strand << " b2.strand " << b2.gr.strand << std::endl;
      return false;
    }

    // b1 is less than b2, or the same but signifies inverted connection
    if ((b1.gr < b2.gr) || (b1.gr.chr == b2.gr.chr && b1.gr.pos1 == b2.gr.pos1)) 
      return true;
    
    std::cerr << b1.gr << " " << b2.gr << std::endl;
    return false;
  }

  void BreakPoint::setRefAlt(faidx_t * main_findex, faidx_t * viral_findex) {
    
    int len;

    assert(main_findex);

    if (evidence != "INDEL") {
      
      // get the reference for BP1
      char * ref1 = faidx_fetch_seq(main_findex, const_cast<char*>(b1.chr_name.c_str()), b1.gr.pos1-1, b1.gr.pos1-1, &len);
      if (!ref1) {
	if (viral_findex)
	  ref1 = faidx_fetch_seq(viral_findex, const_cast<char*>(b1.chr_name.c_str()), b1.gr.pos1-1, b1.gr.pos1-1, &len);
      }
      if (!ref1) {
	std::cerr << "couldn't find reference on BP1 for ref " << b1.chr_name << " in either viral or human" << std::endl;
      }
      
      char * ref2 = faidx_fetch_seq(main_findex, const_cast<char*>(b2.chr_name.c_str()), b2.gr.pos1-1, b2.gr.pos1-1, &len);
      if (!ref2) {
	if (viral_findex)
	  ref2 = faidx_fetch_seq(viral_findex, const_cast<char*>(b2.chr_name.c_str()), b2.gr.pos1-1, b2.gr.pos1-1, &len);
      } 
      if (!ref2) {
	std::cerr << "couldn't find reference on BP2 for ref " << b2.chr_name << " in either viral or human" << std::endl;
      }

      // by convention, set ref to 1 and alt to 2. Gets sorted in VCF creation
      ref = std::string(ref1);
      alt = std::string(ref2);

      if (!ref.length())
	ref = "N";
      if (!alt.length())
	alt = "N";
      
    } else {

      if (insertion.length() && insertion != "x") {

	// reference
	char * refi = faidx_fetch_seq(main_findex, const_cast<char*>(b1.chr_name.c_str()), b1.gr.pos1-1, b1.gr.pos1-1, &len);
	if (!refi) {
	  std::cerr << "couldn't find reference sequence for ref " << b1.chr_name << " for indel, on human reference " << std::endl;
	  return;
	}
	ref = std::string(refi);
	if (!ref.length()) {
	  ref = "N";
	}

	// alt 
	alt = ref + insertion;
      
      // deletion
      } else {	

	// reference
	assert(b2.gr.pos1 - b1.gr.pos1 - 1 >= 0);
	char * refi = faidx_fetch_seq(main_findex, const_cast<char*>(b1.chr_name.c_str()), b1.gr.pos1-1, b2.gr.pos1-2, &len);
	if (!refi) {
	  std::cerr << "couldn't find reference sequence for ref " << b1.chr_name << " for indel, on human reference " << std::endl;
	  return;
	}
	ref = std::string(refi);
	if (!ref.length())
	  ref = "N";
	alt = ref.substr(0,1);
      }
    }
  }

  int BreakPoint::getSpan() const { 
    if (num_align == 1 && insertion == "") {// deletion
      return (abs((int)b1.gr.pos1 - (int)b2.gr.pos1) - 1);
    }
    if (num_align == 1) 
      return (insertion.length()); // insertion
    if (b1.gr.chr == b2.gr.chr)
      return abs((int)b1.gr.pos1-(int)b2.gr.pos1);
    else
      return -1;
  }

  ReducedBreakPoint::ReducedBreakPoint(const std::string &line, bam_hdr_t* h) {
    
    if (!h) {
      std::cerr << "ReducedBreakPoint::ReducedBreakPoint - Must supply non-empty header" << std::endl;
      exit(EXIT_FAILURE);
    }

    //ref = (char*) malloc(line.length() + 1);
    //strcpy(ref, line.c_str());
    //exit(0);

    std::istringstream iss(line);
    std::string val;
    size_t count = 0;

    ref = nullptr;
    alt = nullptr;
    cname = nullptr;
    evidence = nullptr;
    confidence = nullptr;
    insertion = nullptr;
    homology = nullptr;

    float afn, aft;
    std::string ref_s, alt_s, cname_s, insertion_s, homology_s, evidence_s, confidence_s;
    
    std::string chr1, pos1, chr2, pos2, chr_name1, chr_name2; 
    char strand1 = '*', strand2 = '*';
    while (std::getline(iss, val, '\t')) {
      try{
	switch(++count) {
	case 1: chr1 = val; chr_name1 = val; break;
	case 2: pos1 = val; break; 
	case 3: assert(val.length()); strand1 = val.at(0); break;
	case 4: chr2 = val; chr_name2 = val; break;
	case 5: pos2 = val; break; 
	case 6: assert(val.length()); strand2 = val.at(0); break;
	case 7: 
	  ref_s = val;
	  break; 
	case 8: 
	  alt_s = val;
	  break;
	  //case 9: span = stoi(val); break; // automatically calculated
	case 10: 
	  b1 = ReducedBreakEnd(GenomicRegion(chr1, pos1, pos1, h), std::stoi(val), chr_name1); b1.gr.strand = strand1; break;
	case 11:
	  b2 = ReducedBreakEnd(GenomicRegion(chr2, pos2, pos2, h), std::stoi(val), chr_name2); b2.gr.strand = strand2; break;
	case 12: nsplit = std::min((int)255,std::stoi(val)); break;
	case 13: tsplit = std::min((int)255,std::stoi(val)); break;
	case 14: b1.sub_n = std::min((int)255, std::stoi(val)); break;
	case 15: b2.sub_n = std::min((int)255, std::stoi(val)); break;
	case 16: dc.ncount = std::min((int)255, std::stoi(val)); break;
	case 17: dc.tcount = std::min((int)255,std::stoi(val)); break;
	case 18: dc.mapq1 = std::stoi(val); break;  
	case 19: dc.mapq2 = std::stoi(val); break;  
	case 20: ncigar = std::min((int)255, std::stoi(val)); break;
	case 21: tcigar = std::min((int)255, std::stoi(val)); break;
	case 22: 
	  homology_s = val;
	  break; 
	case 23: 
	  //insertion_size = (val == "x") ? 0 : val.length();
	  insertion_s = val;
	  break; 
	case 24: cname_s = val; break;
	case 25: num_align = std::min((int)31, std::stoi(val)); break;
	case 26: 
	  pass = val == "PASS";
	  confidence_s = val;
	  break;
	case 27: 
	  evidence_s = val;
	  indel = val == "INDEL"; 
	  imprecise = val == "DSCRD"; 
	  break; 
	case 28: quality = std::min((int)255,std::stoi(val)); break;
	case 29: secondary = val == "1" ? 1 : 0;
	case 30: somatic_score = std::stof(val); break;
	case 31: pon = std::min(255,std::stoi(val)); break;
	  //case 32: repeat_seq = val; break;
	case 33: ncov = std::min((int)255,std::stoi(val)); break;
	case 34: tcov = std::min((int)255,std::stoi(val)); break;
	case 35: 
	  afn = std::max((float)0, stof(val)); 
	  if (!(afn >= 0 && afn <= 1))
	    afn = 0;
	  af_n = afn * 100; // between 0 and 100
	  break;
	case 36: 
	  aft = std::max((float)0, stof(val)); 
	  if (!(aft >= 0 && aft <= 1))
	    aft = 0;
	  af_t = aft * 100; // between 0 and 100
	  break;
	case 37: blacklist = (val=="1" ? 1 : 0); break;
	case 38: 
	  dbsnp = val != "x";
	  break;
	  //case 39: read_names = val; break;
	}

      } catch(...) {
	std::cerr << "caught stoi error on: " << val << " for count " << count << std::endl;
	std::cerr << line << std::endl;
	exit(1);
      }
    }

    confidence = __string_alloc2char(confidence_s, confidence);
    evidence   = __string_alloc2char(evidence_s, evidence);
    insertion  = __string_alloc2char(insertion_s, insertion);
    homology   = __string_alloc2char(homology_s, homology);
    cname      = __string_alloc2char(cname_s, cname);
    ref        = __string_alloc2char(ref_s, ref);
    alt        = __string_alloc2char(alt_s, alt);

  }


}
