#include "SnowTools/BreakPoint2.h"

#include <getopt.h>
#include <iomanip>
#include <cassert>
#include "SnowTools/gzstream.h"

#include <boost/math/distributions/binomial.hpp>

#define SPLIT_BUFF 8
#define LOD_CUTOFF 8
#define DBCUTOFF 12
#define NODBCUTOFF 8

// define repeats
static std::vector<std::string> repr = {"AAAAAAAA", "TTTTTTTT", "CCCCCCCC", "GGGGGGGG", 
				 "TATATATATATATATA", "ATATATATATATATAT", 
				 "GCGCGCGCGCGCGCGC", "CGCGCGCGCGCGCGCG", 
				 "TGTGTGTGTGTGTGTG", "GTGTGTGTGTGTGTGT", 
				 "TCTCTCTCTCTCTCTC", "CTCTCTCTCTCTCTCT", 
				 "CACACACACACACACA", "ACACACACACACACAC", 
				 "GAGAGAGAGAGAGAGA", "AGAGAGAGAGAGAGAG"};

static std::unordered_map<int, double> ERROR_RATES = {{0,1e-4}, {1,1e-4}, {2, 1e-4}, {3, 1e-4}, {4, 1e-4}, {5, 2e-4}, {6, 5e-4}, {7, 1e-3},
						      {8,2e-3}, {9,3e-3}, {10, 1e-2}, {11, 2e-2}, {12, 3e-5}};

namespace SnowTools {

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
    
    double max_lod = 0;
    for (auto& s : allele) {
      max_lod = std::max(max_lod, s.second.LO);
    }


    ss << b1.chr_name << sep << b1.gr.pos1 << sep << b1.gr.strand << sep 
       << b2.chr_name << sep << b2.gr.pos1 << sep << b2.gr.strand << sep 
       << ref << sep << alt << sep 
       << getSpan() << sep
       << b1.mapq << sep << b2.mapq << sep 
       << b1.nm << sep << b2.nm << sep 
       << dc.mapq1 << sep << dc.mapq2 << sep
       //<< dc.ncount << sep << dc.tcount << sep
       << b1.sub_n << sep << b2.sub_n << sep      
       << (homology.length() ? homology : "x") << sep 
       << (insertion.length() ? insertion : "x") << sep 
       << cname << sep
       << num_align << sep 
       << confidence << sep << evidence << sep
       << quality << sep
       << secondary << sep << somatic_score << sep << somatic_lod << sep 
       << max_lod << sep 
       << pon << sep << (repeat_seq.length() ? repeat_seq : "x") << sep 
       << blacklist << sep << (rs.length() ? rs : "x") << sep 
       << (read_names.length() ? read_names : "x");

    for (auto& a : allele)
      ss << sep << a.second.toFileString();

    return ss.str();

  }
  
  BreakEnd::BreakEnd(const GenomicRegion& g, int mq, const std::string& chr_n) {
    gr = g;
    mapq = mq;
    chr_name = chr_n;
    cpos = -1; nm = -1; matchlen = -1; sub_n = 0; local = false;
  }
  
  ReducedBreakEnd::ReducedBreakEnd(const GenomicRegion& g, int mq, const std::string& chr_n) {
    gr = g;
    mapq = mq;
    chr_name = chr_n;
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

  }
  
  bool BreakPoint::hasDiscordant() const {
    return (dc.ncount || dc.tcount);
  }
  
  bool BreakPoint::hasMinimal() const {
    int count;
    count = dc.tcount + dc.ncount + t.split + n.split;
    if (count >= 2)
      return true;
    return false;
  }
  
  bool BreakPoint::operator==(const BreakPoint &bp) const {
    return (b1.gr == bp.b1.gr && b2.gr == bp.b2.gr && bp.insertion == insertion); 
  }
    
  void BreakPoint::repeatFilter() {

    int span = getSpan();

    if (seq.length() == 0 || span <= 0)
      return;

    std::string r1, r2; //, r3, r4;
    __rep(1, r1);
    repeat_seq = r1.length() > repeat_seq.length() ? r1 : repeat_seq;
    __rep(2, r2);
    repeat_seq = r2.length() > repeat_seq.length() ? r2 : repeat_seq;
    //__rep(3, r3);
    //repeat_seq = r3.length() > repeat_seq.length() ? r3 : repeat_seq;
    // __rep(4, r4);
    //repeat_seq = r4.length() > repeat_seq.length() ? r4 : repeat_seq;

    //std::cerr << seq << " cpos1 " << b1.cpos << " r1 " << r1 << " r2 " << r2 << " r3 " << r3 << " r4 " << r4 << std::endl;

  }


  BreakPoint::BreakPoint(const std::string &line, bam_hdr_t* h) {

    if (!h) {
      std::cerr << "BreakPoint::BreakPoint - Must supply non-empty header" << std::endl;
      exit(EXIT_FAILURE);
    }

    std::istringstream iss(line);
    std::string val;
    size_t count = 0;

    SampleInfo aaa;

    std::string chr1, pos1, chr2, pos2, chr_name1, chr_name2, id; 
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
	case 9: break; //span = stoi(val); break; // automatically calculated
	case 10: 
	  b1 = BreakEnd(GenomicRegion(chr1, pos1, pos1, h), std::stoi(val), chr_name1); b1.gr.strand = strand1; break;
	case 11:
	  b2 = BreakEnd(GenomicRegion(chr2, pos2, pos2, h), std::stoi(val), chr_name2); b2.gr.strand = strand2; break;
	case 12: b1.nm = std::stoi(val); break;
	case 13: b2.nm = std::stoi(val); break;
	  //case 14: dc.ncount = std::stoi(val); break;
	  //case 15: dc.tcount = std::stoi(val); break;
	case 14: dc.mapq1 = std::stoi(val); break;  
	case 15: dc.mapq2 = std::stoi(val); break;  
	case 16: b1.sub_n = std::stoi(val); break;
	case 17: b2.sub_n = std::stoi(val); break;
	case 18: homology = (val == "x" ? "" : val); break;
	case 19: insertion = (val == "x" ? "" : val); break;
	case 20: cname = val; break;
	case 21: num_align = std::stoi(val); break;
	case 22: confidence = val; break;
	case 23: evidence = val; break;
	case 24: quality = std::stoi(val); break;
	case 25: secondary = val == "1";
	case 26: somatic_score = std::stod(val); break;
	case 27: somatic_lod = std::stod(val); break;
	case 28: a.LO = std::stod(val); break;
	case 29: pon = std::stoi(val); break;
	case 30: repeat_seq = val; break;
	case 31: blacklist = (val=="1"); break;
	case 32: rs = val; break;
	case 33: read_names = val; break;
        default: 
	  aaa.indel = evidence == "INDEL";
	  aaa.fromString(val);
	  id += "A";
	  allele[id] = aaa; //id is dummy. just keep in order as came in;
	  break;
	}
      } catch(...) {
	std::cerr << "caught stoi error on: " << val << std::endl;
	std::cerr << line << std::endl;
	exit(1);
      }
    }

  }
  
  void BreakPoint::splitCoverage(BamReadVector &bav) {
    
    int rightbreak1= b1.cpos + SPLIT_BUFF; // read must extend this far right of break1
    int leftbreak1 = b1.cpos - SPLIT_BUFF; // read must extend this far left break2
    int rightbreak2= b2.cpos + SPLIT_BUFF;
    int leftbreak2 = b2.cpos - SPLIT_BUFF;

    int span = getSpan();

    for (auto& j : bav) {

      std::vector<std::string> cnvec = j.GetSmartStringTag("CN");
      
      // for indels, reads with ins alignments to del contig dont count, vice versa
      bool read_should_be_skipped = false;
      if (num_align == 1) {
	std::vector<std::string> cigvec = j.GetSmartStringTag("SC"); // read against contig CIGAR
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
      std::string sample_id = sr.substr(0,4);

      std::string qname;
      // get the alignment position on contig
      int pos = 0, te = 0;
      try {
	std::vector<int> posvec = j.GetSmartIntTag("SL");
	std::vector<int> tevec  = j.GetSmartIntTag("SE");
	
	for (size_t y = 0; y < cnvec.size(); ++y) {
	  if (cnvec[y] == cname) {
	    pos = posvec[y];
	    te  = tevec[y];
	    qname = cname;
	  }
	}
      } catch (...) {
	std::cerr << "error grabbing SE tag for tag " << j.GetZTag("SE") << std::endl;
      }
      assert(qname == cname);
    
      int rightend = te; //seq.length();
      int leftend  = pos;
      bool issplit1 = (leftend <= leftbreak1) && (rightend >= rightbreak1);
      bool issplit2 = (leftend <= leftbreak2) && (rightend >= rightbreak2);

      // add the split reads for each end of the break
      if ( (issplit1 && issplit2 && num_align > 1 ) ||  ((issplit1 || issplit2) && num_align == 1 )) {
	reads.push_back(j);
	std::string rname = j.GetZTag("SR");
	split_reads.insert(rname);
	std::string qn = j.Qname();
	if (qnames.count(qn) && span > 1000)
	  continue; // don't count support if already added and not a short event
	qnames.insert(qn);
	allele[sample_id].supporting_reads.insert(rname);
      }

      // update the counters for each break end
      if (issplit1)
	++b1.split[sample_id];
      if (issplit2)
	++b2.split[sample_id];	
      
      // read spans both ends
      if (issplit1 && issplit2) {
	++allele[sample_id].split;
      }
      
    } // end read loop

      // adjust the alt count
    for (auto& i : allele) {
      i.second.indel = num_align == 1;
      i.second.alt = std::max((int)i.second.supporting_reads.size(), i.second.cigar);
    }
    
  }
  
  std::string BreakPoint::getHashString() const {
    
    bool isdel = insertion.length() == 0;
    //if (isdel) // del breaks are stored as last non-deleted base. CigarMap stores as THE deleted base
    //  pos1++;
    std::string st = std::to_string(b1.gr.chr) + "_" + std::to_string(b1.gr.pos1) + "_" + std::to_string(this->getSpan()) + (isdel ? "D" : "I");
    //std::string st = std::to_string(b1.gr.chr) + "_" + std::to_string(b1.gr.pos1); // + "_" + std::to_string(this->getSpan()) + (isdel ? "D" : "I");
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
  
    std::ostream& operator<<(std::ostream& out, const BreakPoint& b) {
      
      size_t rs_t = std::count(b.rs.begin(), b.rs.end(), 'r');
      
      if (b.isindel) {
	out << ">" << (b.insertion.size() ? "INS: " : "DEL: ") << b.getSpan() << " " << b.b1.gr;
	  //<< " T/N split: " << b.t.split << "/" << b.n.split << " T/N cigar: " 
          //  << b.t.cigar << "/" << b.n.cigar << " T/N Cov " << b.t.cov << "/" << b.n.cov << " DBSNP: " << rs_t;
	for (auto& i : b.allele)
	  out << " " << i.first << ":" << i.second.split;  
      } else {
	out << ": " << b.b1.gr.pointString() << " to " << b.b2.gr.pointString() << " SPAN " << b.getSpan();
	  //<< " T/N split: " << b.t.split << "/" << b.n.split << " T/N disc: " 
	  //  << b.dc.tcount << "/" << b.dc.ncount << " " << b.evidence;
	for (auto& i : b.allele)
	  out << " " << i.first << ":" << i.second.split;  
      }
      
      return out;
      
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
    nm = std::max(b.GetIntTag("NM") - (int)b.MaxInsertionBases() - (int)b.MaxDeletionBases(), 0);
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

	      // add the counts
	      for (auto& c : d.second.counts) {
		allele[c.first].disc = c.second;
		allele[c.first].indel = num_align == 1;

		// add the discordant reads names to supporting reads for each sampleinfo
		for (auto& rr : d.second.reads) {
		  std::string nn = rr.second.GetZTag("SR");
		  allele[nn.substr(0,4)].supporting_reads.insert(nn);
		}
		for (auto& rr : d.second.mates){
		  std::string nn = rr.second.GetZTag("SR");
		  allele[nn.substr(0,4)].supporting_reads.insert(nn);
		}
		
		// adjust the alt counts
		for (auto& aa : allele)
		  aa.second.alt = aa.second.supporting_reads.size();
	      }

	  } 
	
      }
    
  }

  void BreakPoint::__set_evidence() {

    bool isdisc = (dc.tcount + dc.ncount) != 0;
    //bool issplit = (t.split + n.split) != 0;

    if (num_align == 1)
      evidence = "INDEL";
    else if ( isdisc && num_align == 2 )
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

  void BreakPoint::__score_assembly_only() {

    int span = getSpan();
    int this_mapq1 = b1.local ? 60 : b1.mapq;
    int this_mapq2 = b2.local ? 60 : b2.mapq;

    if (seq.length() < 101 + 30)
      confidence = "TOOSHORT";
    else if (a.split < 6 && (span > 1500 || span == -1) && a.split < 7)  // large and inter chrom need 7+
      confidence = "NODISC";
    else if (std::max(b1.mapq, b2.mapq) <= 50 || std::min(b1.mapq, b2.mapq) <= 10) 
      confidence = "LOWMAPQ";
    else if ( std::min(b1.mapq, b2.mapq) <= 30 && a.split < 12 ) 
      confidence = "LOWMAPQ";
    else if (a.split <= 3 && span <= 1500 && span != -1) // small with little split
      confidence = "WEAKASSEMBLY";
    //else if ( (germ && span == -1) || (germ && span > 1000000) ) // super short alignemtns are not to be trusted. Also big germline events
    //  confidence = "WEAKASSEMBLY";
    else if ((b1.sub_n && b1.mapq < 50) || (b2.sub_n && b2.mapq < 50))
      confidence = "MULTIMATCH";
    //else if (n.cov > 500 || t.cov==0 || n.cov == 0)
    //  confidence = "BADREGION";
    else if (secondary)
      confidence = "SECONDARY";
    //else if ( (insertion.length() >= 100 || span < 1000) && (t.split <= 8 || n.split) ) // be extra strict for huge insertions or really low spans
    //  confidence = "WEAKASSEMBLY";
    else
      confidence = "PASS";

  }

  void BreakPoint::__score_somatic() {

    if (evidence == "INDEL") {

      double ratio = n.alt > 0 ? t.alt / n.alt : 100;
      bool immediate_reject = ratio <= 5 || n.alt > 2; // can't call somatic with 3+ reads or <5x more tum than norm ALT

      somatic_lod = n.LO_n;
      if (confidence != "DBSNP")
	somatic_score = somatic_lod > NODBCUTOFF && !immediate_reject;
      else
	somatic_score = somatic_lod > DBCUTOFF && !immediate_reject;	

    } else {
    
      // old model
      double somatic_ratio = 100;
      size_t ncount = std::max(n.split, dc.ncount);
      if (ncount > 0)
	somatic_ratio = (std::max(t.split,dc.tcount)) / ncount;
      
      //double fraction_nclip = ncov == 0 ? 0 : nclip_cov / ncov;
      
      if (somatic_ratio >= 20 && n.split < 2) { // && fraction_nclip < 0.3) // 5 or more nclips is bad
	somatic_score = 1; //somatic_ratio;
	somatic_lod = somatic_ratio;
      } else {
	somatic_lod = 0;
	somatic_score = 0;
      }
    }

    // kill all if too many normal support
    if (n.alt > 2)
      somatic_score = 0;

    // kill if bad ratio
    int r = n.alt > 0 ? t.alt / n.alt : 100;
    if (r < 5)
      somatic_score = 0;

  }

  void BreakPoint::__set_total_reads() {

    // total unique read support
    
    if (evidence == "ASDIS") {
      std::unordered_set<std::string> this_reads_t;
      std::unordered_set<std::string> this_reads_n;
      for (auto& i : split_reads)
	if (!this_reads_t.count(i) && i.at(0) == 't')
	  this_reads_t.insert(i);
      for (auto& i : dc.reads)
	if (!this_reads_t.count(i.first) && i.first.at(0) == 't')
	  this_reads_t.insert(i.first);
      for (auto& i : split_reads)
	if (!this_reads_n.count(i) && i.at(0) == 'n')
	  this_reads_n.insert(i);
      for (auto& i : dc.reads)
	if (!this_reads_n.count(i.first) && i.first.at(0) == 'n')
	  this_reads_n.insert(i.first);
      
      t_reads = this_reads_t.size();
      n_reads = this_reads_n.size();
    } else if (evidence == "INDEL" || evidence == "ASSMB") {
      t_reads = t.alt;
      n_reads = n.alt;
    } else if (evidence == "DSCRD") {
      t_reads = dc.tcount;
      n_reads = dc.ncount;
    } 

  }

  void BreakPoint::__score_assembly_dscrd() {

    int this_mapq1 = b1.mapq;
    int this_mapq2 = b2.mapq;
    int span = getSpan();
    bool germ = dc.ncount > 0 || n.split > 0;


    int max_a_mapq = std::max(this_mapq1, dc.mapq1);
    int max_b_mapq = std::max(this_mapq2, dc.mapq2);

    int min_assm_mapq = std::min(this_mapq1, this_mapq2);
    int max_assm_mapq = std::max(this_mapq1, this_mapq2);
    
    int total_count = t_reads + n_reads; //n.split + t.split + dc.ncount + dc.tcount;

    double min_disc_mapq = std::min(dc.mapq1, dc.mapq2);

    if (max_a_mapq <= 10 || max_b_mapq <= 10)
      confidence = "LOWMAPQ";
    else if (std::max(max_a_mapq, max_b_mapq) <= 30)
      confidence = "LOWMAPQ";
    //if ( (min_disc_mapq <= 10 && min_assm_mapq < 30) || (max_assm_mapq < 40))
    //  confidence = "LOWMAPQ";
    else if ( std::max(t.split, n.split) == 0 || total_count < 4 || (germ && (total_count <= 6) )) // stricter about germline
      confidence = "WEAKASSEMBLY";
    //else if ( std::min(b1.t.split, b2.t.split) >= 5  // take care of really long homoologies (longer than read length)
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

    int min_span = 10e3;
    int disc_count = dc.ncount + dc.tcount;
    if (std::min(dc.mapq1, dc.mapq2) < 20 || std::max(dc.mapq1, dc.mapq2) <= 30) // mapq here is READ mapq (37 std::max)
      confidence = "LOWMAPQ";
    else if (getSpan() >= min_span && std::min(dc.mapq1, dc.mapq2) >= 35 && disc_count >= 7 && n.cov >= 10 & n.cov <= 200 && t.cov >= 10 && t.cov <= 1000)
      confidence = "PASS";
    else if ( disc_count < 8 || (dc.ncount > 0 && disc_count < 15) )  // be stricter about germline disc only
      confidence = "WEAKDISC";
    else if ( getSpan() < min_span)
      confidence = "LOWSPAN";
    else 
      confidence = "PASS";
  }

  void BreakPoint::__score_indel() {

    assert(b1.mapq == b2.mapq);
    
    double ratio = n.alt > 0 ? t.alt / n.alt : 100;

    double max_lod = 0;
    for (auto& s : allele) {
      max_lod = std::max(max_lod, s.second.SLO);
    }

    //bool blacklist_and_low_count = blacklist && (t.split + n.split) < 5 && (t.cigar + n.cigar) < 5;
    //bool blacklist_and_low_AF = (max_af < 0.2 && max_count < 8) && blacklist;

    //if (rs.length() && (af_t < 0.1 || t.cov < 10 || std::max(n.split, n.cigar)))
    //  confidence="DBSNP";
    //if (blacklist && pon > 3)
    //  confidence="GRAYLISTANDPON";
    //else if (blacklist_and_low_count || blacklist_and_low_AF)
    //  confidence="LOWAF";
    //else if (max_af * (double)t.cov < 5)
    //  confidence = "LOWAF";
    //else if ( (max_count < 4 && max_af < 0.2) || (max_count < 2 && b1.mapq < 60) || (max_count < 5 && b1.mapq < 30))
    //  confidence="WEAKASSEMBLY";
    if (b1.mapq < 10)
      confidence="LOWMAPQ";
    else if (max_lod <= LOD_CUTOFF) 
      confidence = "LOWAF";
    //else if ( (max_af < 0.1 && (n.split+n.cigar) > 0)) // || (max_af < 0.30 && n.split == 0 && t.split < 3)) // more strict for germline bc purity is not issue
    //  confidence = "LOWAF";
    //else if (max_af < 0.1)
    //  confidence = "LOWAF";
    //else if (!repeat_seq.empty() && max_af < 0.2)
    //  confidence = "LOWAF";
    //else if (n.cov <= 5)
    //  confidence = "LOWNORMCOV";
    else
      confidence="PASS";
    
  }

  void BreakPoint::scoreBreakpoint() {
    
    // set the evidence (INDEL, DSCRD, etc)
    __set_evidence();
    
    // ensure that breakpoint is oriented
    assert(valid()); 

    // 
    double er = repeat_seq.length() > 10 ? 0.04 : ERROR_RATES[repeat_seq.length()] * 5;
    if (evidence == "INDEL" || true)
      for (auto& i : allele)
	i.second.modelSelection(er);
    
    __combine_alleles();

    // scale the LOD for MAPQ    
    int mapqr1 = b1.local ? std::max(30, b1.mapq) : b1.mapq; // if local, don't drop below 30
    int mapqr2 = b2.local ? std::max(30, b2.mapq) : b2.mapq; // if local (aligns to within window), don't drop below 30
    double scale = (double)( std::min(b1.mapq, b2.mapq) - 2 * b1.nm) / (double)60;
    for (auto& i : allele) 
      i.second.SLO = i.second.LO * scale;
    t.LO = t.SLO * scale;
    n.LO = n.SLO * scale;
   a.LO = a.SLO * scale;
    
    // sanity check
    int split =0;
    for (auto& i : allele) {
      split += i.second.split;
    }
    assert( (split == 0 && t.split == 0 && n.split==0) || (split > 0 && (t.split + n.split > 0)));

    __set_total_reads();

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

    __score_somatic();

    if (confidence == "PASS")
      quality = 99;
    else
      quality = 0;
    
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

      if (ref1)
	free(ref1);
      if (ref2)
	free(ref2);
      
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
      
	if (refi)
	  free(refi);

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
	if (refi)
	  free(refi);
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
    
    std::string chr1, pos1, chr2, pos2, chr_name1, chr_name2, repeat_s; 
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
	case 9: break; //span = stoi(val); break; // automatically calculated
	case 10: 
	  b1 = ReducedBreakEnd(GenomicRegion(chr1, pos1, pos1, h), std::stoi(val), chr_name1); b1.gr.strand = strand1; break;
	case 11:
	  b2 = ReducedBreakEnd(GenomicRegion(chr2, pos2, pos2, h), std::stoi(val), chr_name2); b2.gr.strand = strand2; break;
	case 12: b1.nm = std::stoi(val); break;
	case 13: b2.nm = std::stoi(val); break;
	case 16: b1.sub_n = std::min((int)255, std::stoi(val)); break;
	case 17: b2.sub_n = std::min((int)255, std::stoi(val)); break;
	  //case 14: dc.ncount = std::min((int)255, std::stoi(val)); break;
	  //case 15: dc.tcount = std::min((int)255,std::stoi(val)); break;
	case 14: dc.mapq1 = std::stoi(val); break;  
	case 15: dc.mapq2 = std::stoi(val); break;  
	case 18: 
	  homology_s = val;
	  break; 
	case 19: 
	  //insertion_size = (val == "x") ? 0 : val.length();
	  insertion_s = val;
	  break; 
	case 20: cname_s = val; break;
	case 21: num_align = std::min((int)31, std::stoi(val)); break;
	case 22: 
	  pass = val == "PASS";
	  confidence_s = val;
	  break;
	case 23: 
	  evidence_s = val;
	  indel = val == "INDEL"; 
	  imprecise = val == "DSCRD"; 
	  break; 
	case 24: quality = std::min((int)255,std::stoi(val)); break;
	case 25: secondary = val == "1" ? 1 : 0;
	case 26: somatic_score = std::stod(val); break;
	case 27: somatic_lod = std::stod(val); break;
	case 28: true_lod = std::stod(val); break;
	case 29: pon = std::min(255,std::stoi(val)); break;
	case 30: repeat_s = val; break; // repeat_seq
	case 31: blacklist = (val=="1" ? 1 : 0); break;
	case 32: dbsnp = val != "x"; break;
	case 33: break; //reads
	default:
	  format_s.push_back(val);
	}

      } catch(...) {
	std::cerr << "caught stoi/stod/stof error on: " << val << " for count " << count << std::endl;
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
    repeat     = repeat_s.empty() ? nullptr : __string_alloc2char(repeat_s, repeat);

  }

  double SampleInfo::__log_likelihood(int ref, int alt, double f, double e) {
    
    // mutect log liklihood against error
    double ll = 0;
    const double back_mutate_chance = 1;
    ref = ref <= 0 ? 0 : ref;
    ll += ref * log10(f * e * back_mutate_chance  /* p(alt mut to ref) */ + (1-f)  * (1-e) /* p(ref not mutated) */ ); // ref
    ll += alt * log10(f * (1 - e) /* p(alt not mut) */ + (1-f) * back_mutate_chance * e /* p (ref mut to alt) */); // alt

    return ll;

  }

  void SampleInfo::modelSelection(double er) {

    if (alt >= cov) {
      cov = alt;
    }

    af = cov > 0 ? (double)alt / (double)cov : 1;
    af = af > 1 ? 1 : af;

    // mutect log liklihood against error
    double ll_alt = __log_likelihood(cov - alt, alt, af    , er);
    double ll_err = __log_likelihood(cov - alt, alt, 0, er);
    LO = ll_alt - ll_err;

    //mutetct log likelihood normal
    er = 0.0005; // make this low, so that ALT in REF is rare and NORM in TUM gives low somatic prob
    double ll_alt_norm = __log_likelihood(cov - alt, alt, 0.5, er);
    double ll_ref_norm = __log_likelihood(cov - alt, alt, 0  , er);
    LO_n = ll_ref_norm - ll_alt_norm; // higher number means more likely to be REF than ALT

    // genotype NOT SUPPORTED CURRENTLY
    if (af < 0.2) 
      genotype = "0/1";
    else if (af < 0.8)
      genotype = "0/1";
    else
      genotype = "0/1";

  }
  

  void BreakPoint::addCovs(const std::unordered_map<std::string, STCoverage*>& covs, const std::unordered_map<std::string, STCoverage*>& clip_covs) {
    for (auto& i : covs) 
      allele[i.first].cov = std::max(i.second->getCoverageAtPosition(b1.gr.chr, b1.gr.pos1), i.second->getCoverageAtPosition(b2.gr.chr, b2.gr.pos1));
    for (auto& i : clip_covs)
      allele[i.first].clip_cov = std::max(i.second->getCoverageAtPosition(b1.gr.chr, b1.gr.pos1), i.second->getCoverageAtPosition(b2.gr.chr, b2.gr.pos1));
  }

  std::ostream& operator<<(std::ostream& out, const SampleInfo& a) {
    out << " split: " << a.split << " cigar " << a.cigar << " alt " << a.alt << " clip_cov " << a.clip_cov << " cov " << a.cov << " disc " << a.disc;
    return out;
  }

  SampleInfo operator+(const SampleInfo& a1, const SampleInfo& a2) {

    SampleInfo a;

    a.disc = a1.disc + a2.disc;
    a.split = a1.split + a2.split;
    a.cigar = a1.cigar + a2.cigar;
    a.clip_cov = a1.clip_cov + a2.clip_cov;
    a.cov = a1.cov + a2.cov;

    // add the reads
    for (auto& i : a1.supporting_reads)
      a.supporting_reads.insert(i);
    for (auto& i : a2.supporting_reads)
      a.supporting_reads.insert(i);
    
    a.alt = std::max((int)a.supporting_reads.size(), a.cigar);

    return a;
  }

  void BreakPoint::__combine_alleles() {

    for (auto& s : allele) {
      if (s.first.at(0) == 't') {
	t = t + s.second;
      } else {
	n = n + s.second;
      }
    }

    a = t + n;
    
    double er = repeat_seq.length() > 10 ? 0.04 : ERROR_RATES[repeat_seq.length()] * 5;
    if (evidence == "INDEL") {
      a.modelSelection(er);
      t.modelSelection(er);
      n.modelSelection(er);
    }

  }


  std::string SampleInfo::toFileString() const {

    std::stringstream ss;
    
    if (indel)
      ss << std::setprecision(4) << genotype << ":" << std::max(alt, cigar) << ":" << cov << ":" << GQ << ":" << PL << ":" << split << ":" << cigar 
	 << ":" << LO_n << ":" << LO << ":" << SLO;
    else
      ss << std::setprecision(4) << genotype << ":" << alt << ":" << cov << ":" << GQ << ":" << PL << ":" << split
	 << ":" << disc << ":" << LO_n << ":" << LO << ":" << SLO;
          
    return ss.str();

  }

  void SampleInfo::fromString(const std::string& s) { 

    std::string val;
    int count = 0;
    std::istringstream input(s);
    
    while (std::getline(input, val, ':')) {
      switch(++count) {
      case 6: split = std::stoi(val);
      case 7: cigar = std::stoi(val);
      case 9: LO_n =   std::stod(val);
      case 10: LO = std::stod(val);
      case 11: SLO = std::stod(val);
      case 8: 
	if (indel)
	  cigar = std::stoi(val);
	else
	  disc = std::stoi(val);
      case 2: alt = std::stoi(val);
      case 3: cov = std::stoi(val);
      case 4: GQ = std::stod(val);
      case 5: PL = std::stod(val);
      case 1: genotype = val;
      }
    }
    
  }

  void BreakPoint::__rep(int rep_num, std::string& rseq) {
    
    // move left and right from breakend 1
    int replen = 0;
    std::string mc = "N";
    int curr_replen = 1;
    std::string c = "N";
    rseq = "";

    const int REP_BUFF = 15;
    int i = std::max(0, b1.cpos - REP_BUFF);
    int stop = std::min((int)seq.length() - 1, b1.cpos + REP_BUFF);

    assert(stop - i < 80);

    for (;i < stop; i += rep_num) {

      if (i + rep_num > seq.length()) // shouldn't happen, but ensures no substr errors
	break;

      // terminating
      if (c != seq.substr(i,rep_num) || i >= (stop - rep_num)) {
	if (curr_replen >= replen) {
	  mc = c;
	  replen = curr_replen;
	}
	
	c = seq.substr(i,rep_num );
	curr_replen = rep_num;
      } else {
	curr_replen += rep_num;
      }

    }
    
    for (i = 0; i < replen / rep_num; ++i)
      rseq = rseq + mc;


  }

  std::string BreakEnd::hash() const {
    
    return (std::to_string(gr.chr) + ":" + std::to_string(gr.pos1));
    
  }

}
