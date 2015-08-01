#include "SnowTools/BreakPoint.h"

#include <getopt.h>

#include "SnowTools/gzstream.h"
//#include "SnowTools/HTSTools.h"
//#include "SnowTools/SnowToolsCommon.h"

#define SPLIT_BUFF 8

namespace opt {

  static std::string input_file = "";
  static std::string output_file = "";
  static std::string pon = "";
  static std::string analysis_id = "noid";
  static bool noreads = false;
  static std::string indel_mask = "";

  static std::string ref_index = "/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta";

  static int verbose = 1;
}

static const char* shortopts = "hxi:a:v:q:g:m:";
static const struct option longopts[] = {
  { "help",                    no_argument, NULL, 'h' },
  { "input-bps",               required_argument, NULL, 'i'},
  { "panel-of-normals",        required_argument, NULL, 'q'},
  { "indel-mask",              required_argument, NULL, 'm'},
  { "reference-genome",        required_argument, NULL, 'g'},
  { "analysis-id",             required_argument, NULL, 'a'},
  { "no-reads",                no_argument, NULL, 'x'},
  { "verbose",                 required_argument, NULL, 'v' },
  { NULL, 0, NULL, 0 }
};

// define repeats
static std::vector<std::string> repr = {"AAAAAAAA", "TTTTTTTT", "CCCCCCCC", "GGGGGGGG", 
				 "TATATATATATATATA", "ATATATATATATATAT", 
				 "GCGCGCGCGCGCGCGC", "CGCGCGCGCGCGCGCG", 
				 "TGTGTGTGTGTGTGTG", "GTGTGTGTGTGTGTGT", 
				 "TCTCTCTCTCTCTCTC", "CTCTCTCTCTCTCTCT", 
				 "CACACACACACACACA", "ACACACACACACACAC", 
				 "GAGAGAGAGAGAGAGA", "AGAGAGAGAGAGAGAG"};


static const char *BP_USAGE_MESSAGE =
"Usage: snowman refilter [OPTION] -i bps.txt.gz -o bps.new.txt.gz\n\n"
"  Description: \n"
"\n"
"  General options\n"
"  -v, --verbose                        Select verbosity level (0-4). Default: 1 \n"
"  -h, --help                           Display this help and exit\n"
"  -g, --reference-genome               Path to indexed reference genome to be used by BWA-MEM. Default is Broad hg19 (/seq/reference/...)\n"
"  Required input\n"
"  -i, --input-bps                      Original bps.txt.gz file\n"
"  Optional input\n"                       
"  -q, --panel-of-normals               Panel of normals files generated from snowman pon\n"                       
"  -a, --id-string                      String specifying the analysis ID to be used as part of ID common.\n"
"  -x, --no-reads                       Flag to turn off recording of supporting reads. Setting this flag greatly reduces file size.\n"
"  -m, --indel-mask                     BED-file with blacklisted regions for indel calling. Default none\n"
"\n";

namespace SnowTools {

  // send breakpoint to a string
  std::string BreakPoint::toString() const {
    std::stringstream out;
    
    out << gr1.chr+1 << ":" << gr1.pos1 << "(" << gr1.strand << ")" << "-" 
	<< gr2.chr+1 << ":" << gr2.pos1 << "(" << gr2.strand << ")" <<
      " SPAN: " << getSpan() << " MAPQ: " << 
      mapq1 << "/" << mapq2 << " HOM: " << 
      homology << " INS: " << insertion << " NS: " << 
      nsplit << " TS: " << tsplit << " TD: " << dc.tcount << " ND: " << dc.ncount 
	<< " NC " << ncigar << " TC " << tcigar << " TCOV " << tcov << " NCOV " << ncov << " -- " << cname; // << " isBest: " << isBest;
    return out.str();
  }
  
  // test whether are same break 
  /*bool BreakPoint::sameBreak(BreakPoint &bp) const {
    return bp.refID1 == refID1 && bp.pos1 == pos1 && bp.refID2 == refID2 && bp.pos2 == pos2;
    }*/
  
  // order them
  void BreakPoint::order() {
    
    if (gr1 < gr2)
      return;
    
    GenomicRegion tmp = gr1;
    gr1 = gr2;
    gr2 = tmp;
    unsigned tmptsplit1 = tsplit1;  tsplit1 = tsplit2;  tsplit2 = tmptsplit1;
    unsigned tmpnsplit1 = nsplit1;  nsplit1 = nsplit2;  nsplit2 = tmpnsplit1;
    
  }
  
  // make the file string
  std::string BreakPoint::toFileString(bool noreads) {
    
    std::string sep = "\t";
    std::stringstream ss;
    
    int discordant_tum = dc.tcount;
    int discordant_norm = dc.ncount;
    std::string contig_name = cname;

    /*if (discovar) {
      discordant_tum = disco_tum;
      discordant_norm = disco_norm;
      contig_name = "discovar_" + cname;
      } */
    
    // set the evidence string
    bool isdisc = (dc.tcount + dc.ncount) != 0;
    //bool issplit = (tsplit + nsplit) != 0;
    bool issplit = (tsplit + nsplit) != 0;//dc.m_contig != ""; 
    //if (discovar)
    //  evidence = "DSCVR";
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
    
    if (!evidence.length())
      std::cerr << "num_align " << num_align << " isdisc " << isdisc << " issplit " << issplit << std::endl;
    assert(evidence.length());
    
    confidence = "";
    // low split, high span, no discordant. Suspicous
    
    int disc_count = dc.tcount + dc.ncount;
    int split_count = tsplit + nsplit;
    bool germ = dc.ncount > 0 || nsplit > 0;
    int total_count = disc_count + split_count;
    
    // set the allelic fraction
    double af_n = -1;
    double af_t = -1;
  if (isindel && tcov > 0) 
    af_t = static_cast<double>(std::max(tsplit, tcigar)) / static_cast<double>(tcov);
  if (isindel && ncov > 0) 
    af_n = static_cast<double>(std::max(nsplit, ncigar)) / static_cast<double>(ncov);

  int span = getSpan();

  // check assembly -only ones
  if (num_align > 1 && !hasDiscordant()) {

    if (split_count < 6 && (span > 1500 || span == -1))  // large and inter chrom need 7+
      confidence = "NODISC";
    else if (std::max(mapq1, mapq2) != 60 || std::min(mapq1, mapq2) <= 50) 
      confidence = "LOWMAPQ";
    else if ( split_count <= 3 && (span <= 1500 && span != -1) ) // small with little split
      confidence = "WEAKASSEMBLY";
    else if (/*std::min_end_align_length <= 40 || */(germ && span == -1) || (germ && span > 1000000)) // super short alignemtns are not to be trusted. Also big germline events
      confidence = "WEAKASSEMBLY";
    else
      confidence = "PASS";

    // score ones with both assembly and discordant
  } else if (num_align > 1 && hasDiscordant()) {

    double min_disc_mapq = std::min(dc.getMeanMapq(true), dc.getMeanMapq(false));
    int min_assm_mapq = std::min(mapq1, mapq2);
    //double std::max_disc_mapq = std::max(dc.getMeanMapq(true), dc.getMeanMapq(false));
    int max_assm_mapq = std::max(mapq1, mapq2);
    
    if ( (min_disc_mapq < 10 && min_assm_mapq < 30) || (max_assm_mapq < 40))
      confidence = "LOWMAPQ";
    else if ( std::max(tsplit, nsplit) < 2 || total_count < 4 || (germ && (total_count <= 6) )) // stricter about germline
      confidence = "WEAKASSEMBLY";
    else if ( total_count < 15 && germ && span == -1) // be super strict about germline interchrom
      confidence = "WEAKASSEMBLY";
    else
      confidence = "PASS";

    // score ones with discordant only
  } else if (num_align == 0) {

    if (std::min(mapq1, mapq2) < 30 || std::max(mapq1, mapq2) <= 36) // mapq here is READ mapq (37 std::max)
      confidence = "LOWMAPQ";
    else if (std::min(mapq1, mapq2) == 37 && disc_count >= 6)
      confidence = "PASS";
    else if ( disc_count < 8 || (dc.ncount > 0 && disc_count < 12) )  // be stricter about germline disc only
      confidence = "WEAKDISC";
    else 
      confidence = "PASS";
    
    // indels
  } else if (num_align == 1) {

    double max_af = std::max(af_t, af_n);
    int cigar_count = ncigar+tcigar; 
    int max_count = std::max(split_count, cigar_count);
    bool blacklist_and_low_count = blacklist && (tsplit + nsplit) < 5 && (tcigar + ncigar) < 5;
    bool blacklist_and_low_AF = max_af < 0.3 && blacklist;

    if (rs.length())
      confidence="DBSNP";
    if (blacklist && pon > 0)
      confidence="BLACKLISTANDPON";
    else if (blacklist_and_low_count || blacklist_and_low_AF)
      confidence="BLACKLIST";
    else if ( (max_count < 4 && max_af < 0.3) || (max_count < 3))
      confidence="WEAKASSEMBLY";
    //else if ( cigar_count < 3 && getSpan() < 6)
    //  confidence="WEAKCIGARMATCH";
    else if (std::min(mapq1, mapq2) < 40)
      confidence="LOWMAPQ";
    //else if (seq.find("AAAAAAAAAAA") != string::npos || seq.find("TTTTTTTTTTT") != string::npos || seq.find("TGTGTGTGTGTGTGTGTGTGTGTGTG") != string::npos)
    else if (repeat_seq.length() > 1) // && pon > 0)
      confidence="REPEAT";
    else if ( (max_af < 0.2 && nsplit > 0) || (max_af < 0.15 && nsplit == 0)) // more strict for germline bc purity is not issue
      confidence = "LOWAF";
    else if (ncov <= 5)
      confidence = "LOWNORMCOV";
    else
      confidence="PASS";
  } 

  assert(confidence.length());

  assert(getSpan() > -2);

  if (!noreads) {
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
    
    // print reads to a string
    // delimit with a ,
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
  } else {
    read_names = "";
  }


  // TODO convert chr to string with treader
  ss << gr1.chr+1 << sep << gr1.pos1 << sep << gr1.strand << sep 
     << gr2.chr+1 << sep << gr2.pos1 << sep << gr2.strand << sep 
     << getSpan() << sep
     << mapq1 <<  sep << mapq2 << sep 
     << nsplit << sep << tsplit << sep
     << discordant_norm << sep << discordant_tum << sep
     << ncigar << sep << tcigar << sep
    //<< nall << sep << tall << sep 
     << (homology.length() ? homology : "x") << sep 
     << (insertion.length() ? insertion : "x") << sep 
     << contig_name << sep
     << num_align << sep 
     << confidence << sep << evidence << sep
     << pon << sep << (repeat_seq.length() ? repeat_seq : "x") << sep 
     << ncov << sep << tcov << sep << af_n << sep << af_t << sep
     << blacklist << sep << (rs.length() ? rs : "x") << sep 
     << (read_names.length() ? read_names : "x");

  return ss.str();
  
  }

  // make a breakpoint from a discordant cluster 
  BreakPoint::BreakPoint(const DiscordantCluster& tdc) {
    
    num_align = 0;
    dc = tdc;

    gr1.pos1 = (tdc.m_reg1.strand == '+') ? tdc.m_reg1.pos2 : tdc.m_reg1.pos1;
    gr1.pos2 = gr1.pos1;
    gr2.pos1 = (tdc.m_reg2.strand == '+') ? tdc.m_reg2.pos2 : tdc.m_reg2.pos1;
    gr2.pos2 = gr2.pos1;
    gr1.chr = tdc.m_reg1.chr;
    gr2.chr = tdc.m_reg2.chr;
    gr1.strand = tdc.m_reg1.strand;
    gr2.strand = tdc.m_reg2.strand;
    
    mapq1 = tdc.getMeanMapq(false); 
    mapq2 = tdc.getMeanMapq(true); // mate

    cname = tdc.toRegionString();
    
  }

  bool BreakPoint::hasDiscordant() const {
    return !dc.m_reg1.isEmpty();
  }
  
  
  /** 
   * Has at least two supporting reads
   */
  bool BreakPoint::hasMinimal() const {
    
    int total = tsplit + nsplit + dc.tcount + dc.ncount;
    
    if (total >= 2)
      return true;
    else
      return false;
    
  }
  
  bool BreakPoint::operator==(const BreakPoint &bp) const {
    
    return (gr1 == bp.gr1 && gr2 == bp.gr2); //gr1.chr == bp.gr1.chr && gr1.pos1 == bp.gr1.pos1 && gr2.pos1 == bp.gr2.pos1);
    
  }
  
  
  void runRefilterBreakpoints(int argc, char** argv) {
    
    parseBreakOptions(argc, argv);
    
    opt::output_file = opt::analysis_id + ".filtered.bps.txt.gz";
    if (opt::verbose > 0) {
      std::cout << "Input bps file:  " << opt::input_file << std::endl;
      std::cout << "Output bps file: " << opt::output_file << std::endl;
      std::cout << "Panel of normals file: " << opt::pon << std::endl;
      std::cout << "Indel mask BED:      " << opt::indel_mask << std::endl;
      std::cout << "Analysis id: " << opt::analysis_id << std::endl;
    }
    
    if (!read_access_test(opt::input_file)) {
      std::cerr << "ERROR: Cannot read file " << opt::input_file  << std::endl;
      exit(EXIT_FAILURE);
    }
    
    // load the reference index
    if (opt::verbose > 0)
      std::cout << "attempting to load: " << opt::ref_index << std::endl;
    faidx_t * findex = fai_load(opt::ref_index.c_str());  // load the reference
    
    // open the output file
    igzstream iz(opt::input_file.c_str());
    if (!iz) {
      std::cerr << "Can't read file " << opt::input_file << std::endl;
      exit(EXIT_FAILURE);
    }
    
    // read in the PON
    std::unique_ptr<PON> pmap;
    BreakPoint::readPON(opt::pon, pmap);
    
    // read the indel mask
    //GenomicIntervalTreeMap grm_mask;
    if (!read_access_test(opt::indel_mask) && opt::indel_mask.length()) {
      std::cerr << "indel mask " << opt::indel_mask << " does not exist / is not readable. Skipping indel masking."  << std::endl;
      opt::indel_mask = "";
    }
    
    GenomicRegionCollection<GenomicRegion> grv_mask;
    if (opt::indel_mask.length()) 
      grv_mask.regionFileToGRV(opt::indel_mask, 0, NULL);
    
    ogzstream oz(opt::output_file.c_str(), std::ios::out);
    if (!oz) {
      std::cerr << "Can't write to output file " << opt::output_file << std::endl;
      exit(EXIT_FAILURE);
    }
    
    if (opt::verbose)
      std::cout << "...refiltering variants" << std::endl;
    
    // set the header
    oz << BreakPoint::header() << std::endl;
    
    std::string line;
    //skip the header
    std::getline(iz, line, '\n');
    
    size_t count = 0;
    while (std::getline(iz, line, '\n')) {
      
      if (++count % 10000 == 1 && opt::verbose > 0)
	std::cout << "filtering breakpoint " << count << std::endl;
      
      BreakPoint bp(line);
      
      // check if in panel of normals
      bp.checkPon(pmap);
      
      // check for blacklist
      bp.checkBlacklist(grv_mask);
      
      // check for repeat sequence
      if (bp.gr1.chr < 24) {
	GenomicRegion gr = bp.gr1;
	gr.pad(20);

	// get the reference seqeuence for this piece
	int len;
	std::string chrstring = GenomicRegion::chrToString(gr.chr);
	char * seq = faidx_fetch_seq(findex, const_cast<char*>(chrstring.c_str()), gr.pos1-1, gr.pos2-1, &len);
	if (!seq) {
	  std::cerr <<  "faidx_fetch_seq fail" << std::endl;
	}
	std::string seqr = std::string(seq);
	
	// define repeats
	//std::vector<std::string> repr = {"AAAAAAAA", "TTTTTTTT", "CCCCCCCC", "GGGGGGGG", 
	//				 "TATATATATATATATA", "ATATATATATATATAT", 
	//				 "GCGCGCGCGCGCGCGC", "CGCGCGCGCGCGCGCG", 
	//				 "TGTGTGTGTGTGTGTG", "GTGTGTGTGTGTGTGT", 
	//				 "TCTCTCTCTCTCTCTC", "CTCTCTCTCTCTCTCT", 
	//				 "CACACACACACACACA", "ACACACACACACACAC", 
	//				 "GAGAGAGAGAGAGAGA", "AGAGAGAGAGAGAGAG"};
	//std::cout << "seq " << seqr << std::endl;
	for (auto& i : repr)
	  if (seqr.find(i) != std::string::npos)
	    bp.repeat_seq = i;
	
      }
      
      if (bp.hasMinimal() && bp.gr1.chr < 24 && bp.gr2.chr < 24)
	oz << bp.toFileString(opt::noreads) << std::endl;
    }
    
    oz.close();
    
    // make the VCF files
    /*
      if (opt::verbose)
      std::cout << "...converting filtered.bps.txt.gz to vcf files" << std::endl;
      VCFFile filtered_vcf(opt::output_file, opt::ref_index.c_str(), '\t', opt::analysis_id);
      
      // output the vcfs
      string basename = opt::analysis_id + ".broad-snowman.DATECODE.filtered.";
      filtered_vcf.writeIndels(basename, true); // zip them -> true
      filtered_vcf.writeSVs(basename, true);
    */
  }
  
  void BreakPoint::repeatFilter(faidx_t * f) {

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
    
  }

  BreakPoint::BreakPoint(std::string &line) {
    
    std::istringstream iss(line);
    std::string val;
    size_t count = 0;
    while (std::getline(iss, val, '\t')) {
      try{
	switch(++count) {
	case 1: gr1.chr = stoi(val) - 1; break;
	case 2: gr1.pos1 = stoi(val); gr1.pos2 = gr1.pos1; break;
	case 3: gr1.strand = val.at(0); break;
	case 4: gr2.chr = stoi(val) - 1; break;
	case 5: gr2.pos1 = stoi(val); gr2.pos2 = gr2.pos1; break;
	case 6: gr2.strand = val.at(0); break;
	  //case 7: span = stoi(val); break;
	case 8: mapq1 = stoi(val); break;
	case 9: mapq2 = stoi(val); break;
	case 10: nsplit = stoi(val); break;
	case 11: tsplit = stoi(val); break;
	case 12: dc.ncount = stoi(val); break;
	case 13: dc.tcount = stoi(val); break;
	case 14: ncigar = stoi(val); break;
	case 15: tcigar = stoi(val); break;
	case 16: homology = val; break;
	case 17: insertion = val; break;
	case 18: cname = val; break;
	case 19: num_align = stoi(val); break;
	case 20: confidence = val; break;
	case 21: evidence = val; break;
	case 22: pon = stoi(val); break;
	case 23: repeat_seq = val; break;
	case 24: ncov = stoi(val); break;
	case 25: tcov = stoi(val); break;
	  //case 26: n_af = stod(val); break;
	  //case 27: t_af = stod(val); break;
	case 28: blacklist = (val=="1"); break;
	case 29: read_names = val; break;
	  
	}
      } catch(...) {
	std::cerr << "caught stoi error on: " << val << std::endl;
	std::cerr << line << std::endl;
	exit(1);
      }
    }
    
    if (evidence == "INDEL")
      isindel = true;
    
    //debug
    if (num_align == 1) {
      dc.tcount = 0;
      dc.ncount = 0;
    }
    
    
  }
  
  
  // parse the command line options
  void parseBreakOptions(int argc, char** argv) {
    bool die = false;
    
    if (argc <= 2) 
      die = true;
    
    std::string tmp;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
      std::istringstream arg(optarg != NULL ? optarg : "");
      switch (c) {
      case 'h': die = true; break;
      case 'g': arg >> opt::ref_index; break;
      case 'm': arg >> opt::indel_mask; break;
      case 'i': arg >> opt::input_file; break;
      case 'v': arg >> opt::verbose; break;
      case 'q': arg >> opt::pon; break;
      case 'a': arg >> opt::analysis_id; break;
      case 'x': opt::noreads = true; break;
      }
    }

    if (opt::input_file.length() == 0)
      die = true;
    
    if (die) {
      std::cout << "\n" << BP_USAGE_MESSAGE;
      exit(1);
    }
  }
  
  void BreakPoint::splitCoverage(BamReadVector &bav) {
    
    assert(cname.length());
    
    // make sure we're not double counting
    assert(tsplit == 0 && nsplit == 0 && nsplit1 == 0 && nsplit2 == 0 && tsplit1 == 0 && tsplit2 == 0);
    
    // total number of valid split reads
    tsplit = 0;
    nsplit = 0;
    
    int rightbreak1= cpos1 + SPLIT_BUFF; // read must extend this far right of break1
    int leftbreak1 = cpos1 - SPLIT_BUFF; // read must extend this far left break2
    int rightbreak2= cpos2 + SPLIT_BUFF;
    int leftbreak2 = cpos2 - SPLIT_BUFF;

    for (auto& j : bav) {
      
      std::string sr = j.GetZTag("SR");
      assert(sr.at(0) == 't' || sr.at(0) == 'n');
      
      bool tumor_read = (sr.at(0) == 't');

      int startdum = 0;
      std::string seq = j.QualityTrimmedSequence(4, startdum);
      //r_get_trimmed_seq(j, seq);
      assert(seq.length() > 0);
      
      std::string qname = j.GetSmartStringTag("CN").back(); //GetStringTag(j, "CN").back();
      int pos = j.GetSmartIntTag("SL").back();
      int te = 0;
      try {
	te = j.GetSmartIntTag("SE").back();
      } catch (...) {
	std::cerr << "error grabbing SE tag for tag " << j.GetZTag("SE") << std::endl;
      }
      
      if (qname != cname)
	std::cerr << "qname " << qname << "cname " << cname << std::endl;

      if (qname == cname) { // make sure we're comparing the right alignment
	int rightend = te; //seq.length();
	int leftend  = pos;
	bool issplit1 = (leftend <= leftbreak1) && (rightend >= rightbreak1);
	bool issplit2 = (leftend <= leftbreak2) && (rightend >= rightbreak2);
	
	// add the split reads for each end of the break
	if (issplit1 || issplit2)
	  reads.push_back(j);

	// update the counters for each end
	if (issplit1 && tumor_read)
	  tsplit1++;
	if (issplit1 && !tumor_read)
	  nsplit1++;
	if (issplit2 && tumor_read)
	  tsplit2++;
	if (issplit2 && !tumor_read)
	  nsplit2++;
	
	// read spans both ends
	if ((issplit1 || issplit2) && tumor_read)
	  tsplit++;
	if ((issplit1 || issplit2) && !tumor_read)
	  nsplit++;
	
	/*if (issplit1 || issplit2) {
	  if (tumor_read)
	  tunique++;
	  else
	  nunique++;
	  }*/
	
      }
    }
    //std::cout << "tsplit1 " << tsplit1 << " nsplit1 " << nsplit1 << " tsplit2 " << tsplit2 << " nsplit2 " << nsplit2 << " mapq1 " << mapq1 << " mapq2 " << mapq2 << " cname " << cname << " bp " << *this << std::endl;
  }
  
  std::string BreakPoint::getHashString() const {
    
    int pos1 = gr1.pos1;
    bool isdel = insertion.length() == 0;
    if (isdel) // del breaks are stored as last non-deleted base. CigarMap stores as THE deleted base
      pos1++;
    std::string st = std::to_string(gr1.chr) + "_" + std::to_string(pos1) + "_" + std::to_string(this->getSpan()) + (isdel ? "D" : "I");
    return st;
  }
  
  int BreakPoint::checkPon(std::unique_ptr<PON> &p) {
    
    // only built for indels for now
    if (!isindel || !p)
      return 0;
    
    /*string chr = std::to_string(gr1.chr+1);
      if (chr == "23")
      chr = "X";
      if (chr == "24")
      chr = "Y";
      if (chr == "25")
      chr = "M";
    */
    std::string chr = std::to_string(gr1.chr);
    
    std::string key1, key2, key3, key4, key5;
    int span = getSpan();
    std::string type = (insertion != "") ? "I" : "D";
    key1 = chr + "_" + std::to_string(gr1.pos1-1) + "_" + std::to_string(span) + type;
    key2 = chr + "_" + std::to_string(gr1.pos1  ) + "_" + std::to_string(span) + type;
    key3 = chr + "_" + std::to_string(gr1.pos1+1) + "_" + std::to_string(span) + type;
    key4 = chr + "_" + std::to_string(gr1.pos1-2) + "_" + std::to_string(span) + type;
    key5 = chr + "_" + std::to_string(gr1.pos1+2) + "_" + std::to_string(span) + type;
    
    if (p->count(key1))
      pon = std::max((*p)[key1], pon);
    if (p->count(key2))
      pon = std::max((*p)[key2], pon);
    if (p->count(key3))
      pon = std::max((*p)[key3], pon);
    if (p->count(key4))
      pon = std::max((*p)[key4], pon);
    if (p->count(key5))
      pon = std::max((*p)[key5], pon);
    
    //if (pon > 0)
    //  std::cout << "pon key " << key1 << " val " << pon << std::endl;
    
    return pon;
    
  }
  
  
  void BreakPoint::readPON(std::string &file, std::unique_ptr<PON> &pmap) {
    
    // import the pon
    igzstream izp(file.c_str());
    if (!izp) {
      std::cerr << "Can't read file " << file << std::endl;
      exit(EXIT_FAILURE);
    }
    
    if (opt::verbose)
      std::cout << "...importing PON data" << std::endl;
    pmap = std::unique_ptr<PON>(new PON());
    std::string pval;
    while (std::getline(izp, pval, '\n')) {
      std::istringstream gg(pval);
      std::string tval;
      std::string key;
      //vector<int> scount;
      int sample_count_total = 0;
      //int read_count_total = 0;
      size_t c = 0;
      while (std::getline(gg, tval, '\t')) {
	c++;
	if (c == 1 && tval.length()) {
	  key = tval.substr(1,tval.length() - 1);
	  if (key.at(0) == 'T') // only accumulate normal
	    break;
	}
	else if (tval.length())
	  try { sample_count_total += (stoi(tval) > 0 ? 1 : 0); } catch(...) { std::cerr << "stoi error in PON read with val " << tval << " on line " << pval << std::endl; }
	//else if (tval.length() && c==2)
	//	try { read_count_total += stoi(tval); } catch(...) { std::cerr << "stoi error in PON read with val " << tval << " on line " << pval << std::endl; }
	
      }
      
      if (sample_count_total > 1)
	(*pmap)[key] = sample_count_total;
    }
    
    return;
    
  }
  
  void BreakPoint::addAllelicFraction(STCoverage * t_cov, STCoverage * n_cov) {
    
    if (t_cov)
      tcov = t_cov->getCoverageAtPosition(gr1.chr, gr1.pos1); 
    if (n_cov)
      ncov = n_cov->getCoverageAtPosition(gr1.chr, gr1.pos1); 
  }
 
  std::string BreakPoint::toPrintString() const {
    
    std::stringstream ss;
    
    // set the allelic fraction
    double af_n = -1;
    double af_t = -1;
    std::cerr << cname << " " << tcov_support << " " << ncov_support << std::endl;
    if (isindel && tcov > 0) 
      af_t = static_cast<double>(/*std::max(tsplit, tcigar)*/tcov_support) / static_cast<double>(tcov);
    if (isindel && ncov > 0) 
      af_n = static_cast<double>(/*std::max(nsplit, ncigar)*/ncov_support) / static_cast<double>(ncov);
    
    
    if (isindel)
      ss << ">>>> " << (insertion.size() ? "INSERTION" : "DELETION") <<  " of length " << getSpan() << " at " << gr1 << " contig " << cname 
	 << " T/N split: " << tsplit << "/" << nsplit << " T/N cigar: " 
	 << tcigar << "/" << ncigar << " T/N AF " << af_t << "/" << af_n << " T/N Cov " << tcov << "/" << ncov << " DBSNP " << rs;
    else
      ss << ">>>> STRUCTURAL VAR  at " << gr1.pointString() << " to " << gr2.pointString() << " SPAN " << getSpan() << " contig " << cname 
	 << " T/N split: " << tsplit << "/" << nsplit << " T/N discordant: " 
	 << dc.tcount << "/" << dc.ncount << " evidence " << evidence;
    
    
    return ss.str();
    
  }
  
  void BreakPoint::checkBlacklist(GRC &grv) {
    
    if (!isindel)
      return;
    
    if (grv.findOverlapping(gr1))
      blacklist = true;
    
  }

  void BreakPoint::__combine_with_discordant_cluster(DiscordantClusterMap& dmap)
  {
    int PAD = 400;
    GenomicRegion bp1 = gr1;
    GenomicRegion bp2 = gr2;
    bp1.pad(PAD);
    bp2.pad(PAD);

    for (auto& d : dmap)
      {
	bool bp1reg1 = bp1.getOverlap(d.second.m_reg1) > 0;
	bool bp2reg2 = bp2.getOverlap(d.second.m_reg2) > 0;
	
	//debug
	bool bp1reg2 = bp1.getOverlap(d.second.m_reg2) > 0;
	bool bp2reg1 = bp2.getOverlap(d.second.m_reg1) > 0;
	
	bool pass = (bp1reg1 && bp2reg2) || (bp1reg2 && bp2reg1);

	/*	std::cerr << " gr1 " << gr1 << " gr2 " << gr2 << std::endl << 
	  " m_reg1 " << d.second.m_reg1 << " m_reg2 " << 
	  d.second.m_reg2 << std::endl << " bp1 " << bp1 << 
	  " bp2 " << bp2 << " pass " << pass << 
	  " bp1reg1 " << bp1reg1 << " bp2reg2 " << bp2reg2 << std::endl << 
	  " bp1reg2 " << bp1reg2 << " bp2reg1 " << bp2reg1 << std::endl;
	*/

	if (pass)
	  {
	    // check that we haven't already added a cluster to this breakpoint
	    // if so, chose the one with more tumor support
	    if (dc.isEmpty() || dc.tcount < d.second.tcount) {
	      dc = d.second;
	      d.second.m_contig = cname;
	    } 
	    
	  }
      }
    
  }

}
