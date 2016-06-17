#include "SnowTools/MiniRules2.h"

#include <cassert>
#include "htslib/khash.h"
#include <unordered_set>

//#define QNAME "tumor-10385-DUP-866"
//#define QFLAG 161

//#define DEBUG_MINI 1

namespace SnowTools {

// define what is a valid condition
static const std::unordered_set<std::string> valid = 
  { 
  "duplicate", "supplementary", "qcfail", "hardclip", "fwd_strand",
  "rev_strand", "mate_fwd_strand", "mate_rev_strand", "mapped",
  "mate_mapped", "isize","clip", "phred", "length","nm",
  "mapq", "all", "ff", "xp","fr","rr","rf",
  "ic", "discordant","motif","nbases",
  "ins","del",  "sub",  "subsample", "RG"
};

  static const std::unordered_set<std::string> allowed_region_annots = 
    { "region","pad", "matelink", "exclude"};

  static const std::unordered_set<std::string> allowed_flag_annots = 
    {"duplicate", "supplementary", "qcfail", "hardclip", 
     "fwd_strand", "rev_strand", "mate_fwd_strand", "mate_rev_strand",
     "mapped", "mate_mapped", "ff", "fr", "rr", "rf", "ic"};

bool MiniRules::isValid(BamRead &r) {

  for (auto& it : m_abstract_rules)
    if (it.isValid(r)) 
       return true; // it is includable in at least one. 
      
  return false;

}

  bool __convert_to_bool(const Json::Value& value, const std::string& name) {

    Json::Value null(Json::nullValue);
    Json::Value v = value.get(name, null);
    if (v != null) {
      try {
	if (v.asBool())
	  return true;
	else if (!v.asBool())
	  return false;
      } catch (...) {
	std::cerr << " trouble converting " << name << " to bool on " << value << std::endl;
      }
    }

    return false;
    
  }


  bool __validate_json(const Json::Value value, const std::unordered_set<std::string>& valid_vals) {

    Json::Value null(Json::nullValue);
    for (auto& r : value) {
      bool ok = false;
      for (auto& v : valid_vals) {
	if (r.get(v, null) != null) {
	  ok = true;
	  break;
	}
      }
      if (!ok) {
	std::cerr << " unexpected JSON element of " << value << std::endl;
	std::cerr << " For this scope, must be one of: " << std::endl;
	for (auto& k : valid_vals)
	  std::cerr << "   " << k << std::endl;
	return false;
      }
    }
    
    return true;

  }

// check whether a BamAlignment (or optionally it's mate) is overlapping the regions
// contained in these rules
bool MiniRules::isReadOverlappingRegion(BamRead &r) {

  // if this is a whole genome rule, it overlaps
  if (m_whole_genome)
    return true;

  assert(!m_grv.empty());

  if (m_grv.findOverlapping(GenomicRegion(r.ChrID(), r.Position(), r.PositionEnd())))
    return true;
  
  if (!m_applies_to_mate)
    return false;
  
  if (m_grv.findOverlapping(GenomicRegion(r.MateChrID(), r.MatePosition(), r.MatePosition() + r.Length())))
    return true;

  return false;
}

// checks which rule a read applies to (using the hiearchy stored in m_regions).
// if a read does not satisfy a rule it is excluded.
  bool MiniRulesCollection::isValid(BamRead &r) {

  ++m_count_seen;

  if (m_regions.size() == 0) {
    return "all";
    //std::cerr << "Empty MiniRules" << std::endl;
    //exit(EXIT_FAILURE);
  }

  // need to run all rules if there is an excluder
  if (!m_fall_through)
    for (auto& i : m_regions)
      m_fall_through = m_fall_through || i.excluder;

  size_t which_region = 0;
  size_t which_rule = 0;
  
  // find out which rule it is a part of
  // lower number rules dominate

  bool is_valid = false;
  bool exclude_hit = false; // did we hit excluder rule

  for (auto& it : m_regions) {
    which_rule = 0;
    bool rule_hit = false;
    if (it.isReadOverlappingRegion(r)) { // read overlaps a region

      // empty rule. It's a pass
      if (!it.m_abstract_rules.size()) { 
	is_valid = true;
	if (!rule_hit)
	  ++it.m_count;
	rule_hit = true;
      }
      
      // non-empty, need to check
      for (auto& jt : it.m_abstract_rules) { 
	if (jt.isValid(r)) {

	  // this whole read is valid or not valid
	  // depending on if this is an excluder region
	  if (it.excluder)
	    exclude_hit = true;
	  // if it excluded already, can never be included
	  // if its a normal rule, then exclude_hit is F, it.excluder is F and is_valid get T
	  is_valid = !exclude_hit && !it.excluder; 

	  // update the region counter
	  if (!rule_hit) // first hit for this region?
	    ++it.m_count;

	  rule_hit = true;
	  
	  // update the rule counter within this region
	  ++jt.m_count;
	  
	  if (!m_fall_through)
	    break;
	} 
	++which_rule;
      } // end rules loop
    }
    
    // found a hit in a rule
    if ( (rule_hit && !m_fall_through) || exclude_hit)
      break;
    
    // didnt find hit (or fall through checking), move it up one
    ++which_region;
  }
  
  // isn't in a rule or it never satisfied one. Remove
  if (!is_valid)
    return false; 
  
  ++m_count;
  return true; 
  
}

// convert a region BED file into an interval tree map
void MiniRules::setRegionFromFile(const std::string& file) {
  
  m_region_file = file;

  // parse if it's a file
  if (SnowTools::read_access_test(file))
    m_grv.regionFileToGRV(file, pad, mrc->h);
  else {
    if (mrc->h) {
      GenomicRegion gr(file, mrc->h);
      gr.pad(pad);
      m_grv.add(gr);
    } else {
      std::cerr << "!!!!!!!!MiniRules region parsing: Header from BAM not set!!!!!!!!!" << std::endl;
    }
  }

 
  if (m_grv.empty()) {
    std::cerr << "!!!!!!!!!!!!!!!!!!" << std::endl;
    std::cerr << "!!!!!!!!!!!!!!!!!!" << std::endl;
    std::cerr << "!!!!!!!!!!!!!!!!!!" << std::endl;
    std::cerr << "Warning: No regions detected in region/file: " << file << std::endl;
    std::cerr << "!!!!!!!!!!!!!!!!!!" << std::endl;
    std::cerr << "!!!!!!!!!!!!!!!!!!" << std::endl;
    std::cerr << "!!!!!!!!!!!!!!!!!!" << std::endl;
    return;
  }
  
  // create the interval tree 
  m_grv.createTreeMap();

  return;
}

  MiniRulesCollection::MiniRulesCollection(const std::string& file, bam_hdr_t *b)
  {
    // set the header
    h = b;
    __construct_MRC(file);
  }

  MiniRulesCollection::MiniRulesCollection(const std::string& file)
  {
    __construct_MRC(file);
  }

// constructor to make a MiniRulesCollection from a rules file.
// This will reduce each individual BED file and make the 
// GenomicIntervalTreeMap
void MiniRulesCollection::__construct_MRC(const std::string& script) {

  // parse the rules file
  //std::cerr << " FILE " << file << std::endl;
  //const std::string document = GetScriptContents(file);

  // set up JsonCPP reader and attempt to parse script
  Json::Value root;
  Json::Reader reader;
  if ( !reader.parse(script, root)) {
    
    if (script.empty()) {
      std::cerr << "JSON script is empty. Setting default to filter all reads" << std::endl;
      return;
    }
      
    // use built-in error reporting mechanism to alert user what was wrong with the script
    std::cerr  << "ERROR: failed to parse JSON script" << std::endl;
    std::cerr << script << std::endl;
    exit(EXIT_FAILURE);
  }

  Json::Value null(Json::nullValue);

  int level = 1;

  Json::Value glob = root.removeMember("global");
  if (!glob.isNull()) {
    rule_all.parseJson(glob);
  }

  // iterator over regions
  for (auto& regions : root) {
      
    MiniRules mr;
    mr.mrc = this;
    
    // add global rules (if there are any)
    //for (auto& a : all_rules)
    // mr.m_abstract_rules.push_back(a);

    // check if mate applies
    mr.m_applies_to_mate = __convert_to_bool(regions, "matelink");

    // check for region padding
    mr.pad = regions.get("pad", 0).asInt();

    // set the region
    std::string reg;
    Json::Value v  = regions.get("region", null);
    if (v != null) 
      reg = v.asString();

    // actually parse the region
    if (reg == "WG" || reg.empty())
      mr.m_whole_genome = true;
    else
      mr.setRegionFromFile(reg);

    // check if its excluder region
    mr.excluder = false; // default is no exclude
    v = regions.get("exclude", null);
    if (v != null) 
      mr.excluder = v.asBool();

    // set the rules
    v = regions.get("rules", null);
    // loop through the rules
    for (auto& vv : v) {
      if (vv != null) {
    	AbstractRule ar = rule_all; // always start with the global rule
    	ar.parseJson(vv);
	// add the rule to the region
	mr.m_abstract_rules.push_back(ar);
      }
    }

    // check that the regions have at least one rule
    // if it it doesn't, give it the global WG all
    if (!mr.m_abstract_rules.size())
      mr.m_abstract_rules.push_back(rule_all);

    mr.m_level = level++;
    mr.id = std::to_string(level);

    m_regions.push_back(mr);

  }
  
  // check that there is at least one non-excluder region. 
  // if not, give global includer
  bool has_includer = false;
  for (auto& kk : m_regions)
    if (!kk.excluder)
      has_includer = true;
  if (!has_includer) {
    MiniRules mr;
    mr.m_whole_genome = true;
    mr.m_abstract_rules.push_back(rule_all);
    mr.mrc = this; // set the pointer to the collection
    m_regions.push_back(mr);
  }

}

// print the MiniRulesCollection
std::ostream& operator<<(std::ostream &out, const MiniRulesCollection &mr) {

  out << "----------MiniRulesCollection-------------" << std::endl;
  out << "--- counting all rules (fall through): " << (mr.m_fall_through ? "ON" : "OFF") << std::endl;

  for (auto& it : mr.m_regions)
    out << it;
  out << "------------------------------------------";
  /*  std::cerr << "--- Rule counts " << std::endl;
  for (auto& g : mr.m_regions)
    for (auto& r : g.m_abstract_rules)
      std::cerr << g.id << "\t" << g.m_count << "\t" << r.id << "\t" << r.m_count << std::endl;
  */
  return out;

}

// print a MiniRules information
std::ostream& operator<<(std::ostream &out, const MiniRules &mr) {
  
  std::string file_print = mr.m_whole_genome ? "WHOLE GENOME" : mr.m_region_file;
  out << (mr.excluder ? "--Exclude Region: " : "--RegionInput:") << file_print;
  if (!mr.m_whole_genome) {
    //out << " --Size: " << AddCommas<int>(mr.m_width); 
    out << " --Pad: " << mr.pad;
    out << " --Include Mate: " << (mr.m_applies_to_mate ? "ON" : "OFF");
    if (mr.m_grv.size() == 1)
      out << " --Region : " << mr.m_grv[0] << std::endl;
    else
      out << " --Region: " << mr.m_grv.size() << " regions" << std::endl;      
  } else {
    out << std::endl;
  }

  for (auto& it : mr.m_abstract_rules) 
    out << it << std::endl;
  
  return out;
}

// merge all of the intervals into one and send to a bed file
void MiniRulesCollection::sendToBed(std::string file) {

  std::ofstream out(file);
  if (!out) {
    std::cerr << "Cannot write BED file: " << file << std::endl;
    return;
  }

  // make a composite from all the rules
  GenomicRegionCollection<GenomicRegion> comp;
  for (auto& it : m_regions)
    comp.concat(it.m_grv);
  
  // merge it down
  comp.mergeOverlappingIntervals();

  // send to BED file
  out << comp.sendToBED();
  out.close();

  return;
}

  
  bool Flag::parseJson(const Json::Value& value, const std::string& name) {

    if (value.isMember(name.c_str())) {
      __convert_to_bool(value, name) ? setOn() : setOff();
      return true;
    }
    
    return false;

  }


  void FlagRule::parseJson(const Json::Value& value) {

    // have to set the na if find flag so that rule knows it cant skip checking
    if (dup.parseJson(value, "duplicate")) na = false;
    if (supp.parseJson(value, "supplementary")) na = false;
    if (qcfail.parseJson(value, "qcfail")) na = false;
    if (hardclip.parseJson(value, "hardclip")) na = false;
    if (fwd_strand.parseJson(value, "fwd_strand")) na = false;
    if (mate_rev_strand.parseJson(value, "mate_rev")) na = false;
    if (mate_fwd_strand.parseJson(value, "mate_fwd")) na = false;
    if (mate_mapped.parseJson(value, "mate_mapped")) na = false;
    if (mapped.parseJson(value, "mapped")) na = false;
    if (ff.parseJson(value, "ff")) na = false;
    if (fr.parseJson(value, "fr")) na = false;
    if (rf.parseJson(value, "rf")) na = false;
    if (rr.parseJson(value, "rr")) na = false;
    if (ic.parseJson(value, "ic")) na = false;

  }
  
  void Range::parseJson(const Json::Value& value, const std::string& name) {
    Json::Value null(Json::nullValue);
    Json::Value v = value.get(name, null);

    if (v != null) {
      if (v.size() > 2) {
	std::cerr << " ERROR. Not expecting array size " << v.size() << " for Range " << name << std::endl;
      } else {
	every = false;
	none = false;
	inverted = false;

	
	if (v.isArray()) {
	  min = v[0].asInt();
	  max = v[1].asInt();
	} else if (v.isInt()) {
	  min = v.asInt();
	  max = INT_MAX;
	} else if (v.isBool()) {
	  min = v.asBool() ? 1 : INT_MAX; // if true, [1,MAX], if false [MAX,1] (not 1-MAX)
	  max = v.asBool() ? INT_MAX : 1;
	} else {
	  std::cerr << "Unexpected type for range flag: " << name << std::endl;
	  exit(EXIT_FAILURE);
	}

	if (min > max) {
	  inverted = true;
	  std::swap(min, max); // make min always lower
	}
      }
	
    }
  }

  void AbstractRule::parseJson(const Json::Value& value) {

    // verify that it has appropriate values
    for (auto& i : value.getMemberNames()) {
      if (!valid.count(i)) {
	std::cerr << "Invalid key value in JSON: " << i << std::endl;
	exit(EXIT_FAILURE);
      }
    }
	
    // parse read group
    const std::string rg = "RG";
    if (value.isMember(rg.c_str())) {
      Json::Value null(Json::nullValue);
      Json::Value v = value.get(rg, null);
      assert(v != null);
      read_group = v.asString();
    }
      
    // parse the flags
    fr.parseJson(value);
    
    isize.parseJson(value, "isize");
    mapq.parseJson(value, "mapq");
    len.parseJson(value, "length");
    clip.parseJson(value, "clip");
    phred.parseJson(value, "phred");
    nbases.parseJson(value, "nbases");
    ins.parseJson(value, "ins");
    del.parseJson(value, "del");
    nm.parseJson(value, "nm");
    xp.parseJson(value, "xp");
    
    // parse the subsample data
    parseSubLine(value);
    
#ifndef __APPLE__
    // parse aho corasick file, if not already inheretid
    if (!atm) {
	parseSeqLine(value);
	if (atm) {
	  ac_automata_finalize(atm);
	  std::cerr << "Done generating Aho-Corasick tree" << std::endl;  
	}
      }
#endif

  }


// main function for determining if a read is valid
  bool AbstractRule::isValid(BamRead &r) {
    

#ifdef QNAME
    if (r.Qname() == QNAME && (r.AlignmentFlag() == QFLAG || QFLAG == -1)) {
      std::cerr << "MINIRULES Read seen " << " ID " << id << " " << r << std::endl;
      std::cerr << " PAIR ORIENTATION " << r.PairOrientation() << " PAIR MAPPED FLAG " << r.ReverseFlag() << " MR " << r.MateReverseFlag() << " POS < MAT " << (r.Position() < r.MatePosition()) << std::endl;
    }
#endif

    // check if its keep all or none
    if (isEvery())
      return true;
    
    // check if it is a subsample
    if (subsam_frac < 1) {
      uint32_t k = __ac_Wang_hash(__ac_X31_hash_string(r.QnameChar()) ^ subsam_seed);
      if ((double)(k&0xffffff) / 0x1000000 >= subsam_frac) 
	return false;
    }
    /*if (subsample < 100) {
      int randn = (rand() % 100); // random number between 1 and 100
      if (subsample < randn)
      return false;
      }*/
    
    // check if is discordant
    bool isize_pass = isize.isValid(r.FullInsertSize());

#ifdef QNAME
    if (r.Qname() == QNAME  && (r.AlignmentFlag() == QFLAG || QFLAG == -1))
      std::cerr << "MINIRULES isize_pass " << isize_pass << " " << " ID " << id << " " << r << std::endl;
#endif
    
    if (!isize_pass) {
      return false;
    }
    
    // check for valid read name 
    if (!read_group.empty()) {
      std::string RG = r.ParseReadGroup();
      if (!RG.empty() && RG != read_group)
	return false;
    }

    // check for valid mapping quality
    if (!mapq.isEvery())
      if (!mapq.isValid(r.MapQuality())) 
	return false;

#ifdef QNAME
    if (r.Qname() == QNAME && (r.AlignmentFlag() == QFLAG || QFLAG == -1) )
      std::cerr << "MINIRULES mapq pass " << " ID " << id << " " << r << std::endl;
#endif
    
    // check for valid flags
    if (!fr.isValid(r))
      return false;

#ifdef QNAME
    if (r.Qname() == QNAME && (r.AlignmentFlag() == QFLAG || QFLAG == -1))
      std::cerr << "MINIRULES flag pass " << " " << id << " " << r << std::endl;
#endif
    
    // check the CIGAR
    if (!ins.isEvery() || !del.isEvery()) {
      if (!ins.isValid(r.MaxInsertionBases()))
	return false;
      if (!del.isValid(r.MaxDeletionBases()))
	return false;
    }


#ifdef QNAME
    if (r.Qname() == QNAME && (r.AlignmentFlag() == QFLAG || QFLAG == -1))
      std::cerr << "MINIRULES cigar pass " << " ID " << id << " "  << r << std::endl;
#endif

    
    // if we dont need to because everything is pass, just just pass it
    bool need_to_continue = !nm.isEvery() || !clip.isEvery() || !len.isEvery() || !nbases.isEvery() || atm_file.length() || !xp.isEvery();
#ifdef QNAME
    if (r.Qname() == QNAME && (r.AlignmentFlag() == QFLAG || QFLAG == -1))
      std::cerr << "MINIRULES is need to continue " << need_to_continue << " ID " << id << " " << r << std::endl;
    if (r.Qname() == QNAME && (r.AlignmentFlag() == QFLAG || QFLAG == -1) && !need_to_continue)
      std::cerr << "****** READ ACCEPTED ****** " << std::endl;
      
#endif

    if (!need_to_continue)
      return true;

#ifdef QNAME
    if (r.Qname() == QNAME && (r.AlignmentFlag() == QFLAG || QFLAG == -1))
      std::cerr << "MINIRULES moving on. ID " << id << " " << r << std::endl;
#endif
    
    // now check if we need to build char if all we want is clip
    unsigned clipnum = 0;
    if (!clip.isEvery()) {
      clipnum = r.NumClip();

#ifdef QNAME
    if (r.Qname() == QNAME && (r.AlignmentFlag() == QFLAG || QFLAG == -1))
      std::cerr << "MINIRULES CLIPNUM " << clipnum << " " << r << " NM&&LEN&&CLIP " << (nm.isEvery() && len.isEvery() && !clip.isValid(clipnum)) << " NM " << nm.isEvery() << " CLIPVALID " << clip.isValid(clipnum) << " LEN " << len.isEvery() << "  CLUIP " << clip << std::endl;
#endif

      if (nm.isEvery() && len.isEvery() && !clip.isValid(clipnum)) // if clip fails, its not going to get better by trimming. kill it now before building teh char data
	return false;
    }

#ifdef QNAME
    if (r.Qname() == QNAME && (r.AlignmentFlag() == QFLAG || QFLAG == -1))
      std::cerr << "MINIRULES pre-phred filter clip pass. ID " << id  << " " << r << std::endl;
#endif
    
    // check for valid NM
    if (!nm.isEvery()) {
      int32_t nm_val = r.GetIntTag("NM");
      if (!nm.isValid(nm_val))
	return false;
    }
    
    // trim the read, then check length
    int32_t new_len, new_clipnum; 
    if (phred.isEvery()) {
      new_len = r.Length(); //a.QueryBases.length();
      new_clipnum = clipnum;
    }
    
    if (!phred.isEvery()) {
      
      int32_t startpoint = 0, endpoint = 0;
      r.QualityTrimmedSequence(phred.min, startpoint, endpoint);
      new_len = endpoint - startpoint;
      
      if (endpoint != -1 && new_len < r.Length() && new_len > 0 && new_len - startpoint >= 0 && startpoint + new_len <= r.Length()) { 
	try { 
	  r.AddZTag("GV", r.Sequence().substr(startpoint, new_len));
	  assert(r.GetZTag("GV").length());
	} catch (...) {
	  std::cerr << "Subsequence failure with sequence of length "  
		    << r.Sequence().length() << " and startpoint "
		    << startpoint << " endpoint " << endpoint 
		    << " newlen " << new_len << std::endl;
	}
	// read is fine
      } else {
	r.AddZTag("GV", r.Sequence());
      }

#ifdef QNAME
      if (r.Qname() == QNAME && (r.AlignmentFlag() == QFLAG || QFLAG == -1))
	std::cerr << "MINIRULES pre-phred filter. ID " << id << " start " << startpoint << " endpoint " << endpoint << " " << r << std::endl;
#endif

      // all the bases are trimmed away 
      if (endpoint == -1 || new_len == 0)
	return false;
      
      new_clipnum = std::max(0, static_cast<int>(clipnum - (r.Length() - new_len)));
      
      // check the N
      if (!nbases.isEvery()) {
	
	size_t n = 0;
	//assert((new_len + start - 1) < r.Length()); //debug
	//r_count_sub_nbases(r, n, start, new_len + start); // TODO factor in trimming
	n = r.CountNBases();
	if (!nbases.isValid(n))
	  return false;
      }
      
    }
    
    // check the N if we didn't do phred trimming
    if (!nbases.isEvery() && phred.isEvery()) {
      size_t n = r.CountNBases();
      if (!nbases.isValid(n))
	return false;
    }

#ifdef QNAME
    if (r.Qname() == QNAME && (r.AlignmentFlag() == QFLAG || QFLAG == -1))
      std::cerr << "MINIRULES NBASES PASS. id: " << id << r << " new_len " << new_len << " new_clipnum" << new_clipnum << std::endl;
#endif

    
    // check for valid length
    if (!len.isValid(new_len))
      return false;

#ifdef QNAME
    if (r.Qname() == QNAME && (r.AlignmentFlag() == QFLAG || QFLAG == -1) )
      std::cerr << "MINIRULES LENGTH PASS. id: " << id << r << " newlen " << new_len << std::endl;
#endif
    
    // check for valid clip
    if (!clip.isValid(new_clipnum))
      return false;

#ifdef QNAME
    if (r.Qname() == QNAME && (r.AlignmentFlag() == QFLAG || QFLAG == -1))
      std::cerr << "MINIRULES CLIP AND LEN PASS. ID " << id << " "  << r << std::endl;
#endif

    
    // check for secondary alignments
    if (!xp.isEvery()) 
      if (!xp.isValid(r.CountSecondaryAlignments()))
	return false;
    
    //#ifdef HAVE_AHOCORASICK_AHOCORASICK_H
#ifndef __APPLE__
    if (atm_file.length()) {
      bool m = ahomatch(r);
      if ( (!m && !atm_inv) || (m && atm_inv) )
	return false;
    }
#endif
    
#ifdef QNAME
    if (r.Qname() == QNAME && (r.AlignmentFlag() == QFLAG || QFLAG == -1))
      std::cerr << "****** READ ACCEPTED ****** " << std::endl;
#endif


    return true;
  }
  
  bool FlagRule::isValid(BamRead &r) {
    
    if (isEvery())
      return true;
    
    if (!dup.isNA()) 
      if ((dup.isOff() && r.DuplicateFlag()) || (dup.isOn() && !r.DuplicateFlag()))
      return false;
  if (!supp.isNA()) 
    if ((supp.isOff() && r.SecondaryFlag()) || (supp.isOn() && !r.SecondaryFlag()))
      return false;
  if (!qcfail.isNA())
    if ((qcfail.isOff() && r.QCFailFlag()) || (qcfail.isOn() && !r.QCFailFlag()))
      return false;
  if (!mapped.isNA())
    if ( (mapped.isOff() && r.MappedFlag()) || (mapped.isOn() && !r.MappedFlag()))
      return false;
  if (!mate_mapped.isNA())
    if ( (mate_mapped.isOff() && r.MateMappedFlag()) || (mate_mapped.isOn() && !r.MateMappedFlag()) )
      return false;

  // check for hard clips
  if (!hardclip.isNA())  {// check that we want to chuck hard clip
    if (r.CigarSize() > 1) {
      bool ishclipped = r.NumHardClip() > 0;
      if ( (ishclipped && hardclip.isOff()) || (!ishclipped && hardclip.isOn()) )
	return false;
    }
  }

  // check for orientation
  // check first if we need to even look for orientation
  bool ocheck = !ff.isNA() || !fr.isNA() || !rf.isNA() || !rr.isNA() || !ic.isNA();

  // now its an orientation pair check. If not both mapped, chuck
  if (!r.PairMappedFlag() && ocheck)
    return false;

  if ( ocheck ) {

    //    bool first = r.Position() < r.MatePosition();
    //bool bfr = (first && (!r.ReverseFlag() && r.MateReverseFlag())) || (!first &&  r.ReverseFlag() && !r.MateReverseFlag());
    //bool brr = r.ReverseFlag() && r.MateReverseFlag();
    //bool brf = (first &&  (r.ReverseFlag() && !r.MateReverseFlag())) || (!first && !r.ReverseFlag() &&  r.MateReverseFlag());
    //bool bff = !r.ReverseFlag() && !r.MateReverseFlag();
      
    bool bic = r.Interchromosomal();

    int PO = r.PairOrientation();

    // its FR and it CANT be FR (off) or its !FR and it MUST be FR (ON)
    // orienation not defined for inter-chrom, so exclude these with !ic
    if (!bic) { // PROCEED IF INTRA-CHROMOSOMAL
      //if ( (bfr && fr.isOff()) || (!bfr && fr.isOn())) 
      if ( (PO == FRORIENTATION && fr.isOff()) || (PO != FRORIENTATION && fr.isOn())) 
	return false;
      //if ( (brr && rr.isOff()) || (!brr && rr.isOn())) 
      if ( (PO == RRORIENTATION && rr.isOff()) || (PO != RRORIENTATION && rr.isOn())) 
	return false;
      //if ( (brf && rf.isOff()) || (!brf && rf.isOn())) 
      if ( (PO == RFORIENTATION && rf.isOff()) || (PO != RFORIENTATION&& rf.isOn())) 
	return false;
      //if ( (bff && ff.isOff()) || (!bff && ff.isOn())) 
      if ( (PO == FFORIENTATION && ff.isOff()) || (PO != FFORIENTATION && ff.isOn())) 
	return false;
    }
    if ( (bic && ic.isOff()) || (!bic && ic.isOn()))
      return false;
      
  }

  return true;
  
}

// define how to print
std::ostream& operator<<(std::ostream &out, const AbstractRule &ar) {

  out << "  Rule: " << ar.name << " -- ";;
  if (ar.isEvery()) {
    out << "  ALL";
  } else if (ar.isNone()) {
    out << "  KEEPING NONE";  
  } else {
    if (!ar.read_group.empty())
      out << "Read Group: " << ar.read_group << " -- ";
    if (!ar.isize.isEvery())
      out << "isize:" << ar.isize << " -- " ;
    if (!ar.mapq.isEvery())
      out << "mapq:" << ar.mapq << " -- " ;
    if (!ar.len.isEvery())
      out << "length:" << ar.len << " -- ";
    if (!ar.clip.isEvery())
      out << "clip:" << ar.clip << " -- ";
    if (!ar.phred.isEvery())
      out << "phred:" << ar.phred << " -- ";
    if (!ar.nm.isEvery())
      out << "nm:" << ar.nm << " -- ";
    if (!ar.xp.isEvery())
      out << "xp:" << ar.xp << " -- ";
    if (!ar.nbases.isEvery())
      out << "nbases:" << ar.nbases << " -- ";
    if (!ar.ins.isEvery())
      out << "ins:" << ar.ins << " -- ";
    if (!ar.del.isEvery())
      out << "del:" << ar.del << " -- ";
    if (ar.subsam_frac < 1)
      out << "sub:" << ar.subsam_frac << " -- ";
#ifndef __APPLE__
    //#ifdef HAVE_AHOCORASICK_AHOCORASICK_H
    if (ar.atm_file != "")
      out << (ar.atm_inv ? "NOT " : "") << "matching on " << ar.atm_count << " motifs from " << ar.atm_file << " -- ";
#endif
    out << ar.fr;
  }
  return out;
}

// define how to print
std::ostream& operator<<(std::ostream &out, const FlagRule &fr) {

  if (fr.isEvery()) {
    out << "Flag: ALL";
    return out;
  } 

  std::string keep = "Flag ON: ";
  std::string remo = "Flag OFF: ";

  if (fr.dup.isOff())
    remo += "duplicate,";
  if (fr.dup.isOn())
    keep += "duplicate,";

  if (fr.supp.isOff())
    remo += "supplementary,";
  if (fr.supp.isOn())
    keep += "supplementary,";

  if (fr.qcfail.isOff())
    remo += "qcfail,";
  if (fr.qcfail.isOn())
    keep += "qcfail,";

  if (fr.hardclip.isOff())
    remo += "hardclip,";
  if (fr.hardclip.isOn())
    keep += "hardclip,";

  if (fr.paired.isOff())
    remo += "paired,";
  if (fr.paired.isOn())
    keep += "paired,";



  if (fr.ic.isOff())
    remo += "ic,";
  if (fr.ic.isOn())
    keep += "ic,";

  if (fr.ff.isOff())
    remo += "ff,";
  if (fr.ff.isOn())
    keep += "ff,";

  if (fr.fr.isOff())
    remo += "fr,";
  if (fr.fr.isOn())
    keep += "fr,";

  if (fr.rr.isOff())
    remo += "rr,";
  if (fr.rr.isOn())
    keep += "rr,";

  if (fr.rf.isOff())
    remo += "rf,";
  if (fr.rf.isOn())
    keep += "rf,";

  if (fr.mapped.isOff())
    remo += "mapped,";
  if (fr.mapped.isOn())
    keep += "mapped,";

  if (fr.mate_mapped.isOff())
    remo += "mate_mapped,";
  if (fr.mate_mapped.isOn())
    keep += "mate_mapped,";

  keep = keep.length() > 10 ? keep.substr(0, keep.length() - 1) : ""; // remove trailing comment
  remo = remo.length() > 10 ? remo.substr(0, remo.length() - 1) : ""; // remove trailing comment
  
  if (!keep.empty() && !remo.empty())
    out << keep << " -- " << remo;
  else if (!keep.empty())
    out << keep;
  else
    out << remo;

  return out;
}

// define how to print
std::ostream& operator<<(std::ostream &out, const Range &r) {
  if (r.isEvery())
    out << "all";
  else if (r.min == 1 && r.max == INT_MAX && !r.inverted)
    out << "ALL";
  else if (r.min == 1 && r.max == INT_MAX && r.inverted)
    out << "NONE";
  else
    out << (r.inverted ? "NOT " : "") << "[" << r.min << "," << (r.max == INT_MAX ? "MAX" : std::to_string(r.max))  << "]";
  return out;
}

  //#ifdef HAVE_AHOCORASICK_AHOCORASICK_H
#ifndef __APPLE__
// check if a string contains a substring using Aho Corasick algorithm
//bool AbstractRule::ahomatch(const string& seq) {
bool AbstractRule::ahomatch(BamRead &r) {

  // make into Ac strcut
  std::string seq = r.Sequence();
  AC_TEXT_t tmp_text = {seq.c_str(), static_cast<unsigned>(seq.length())};
  ac_automata_settext (atm, &tmp_text, 0);

  // do the check
  AC_MATCH_t * matchp;  
  matchp = ac_automata_findnext(atm);

  if (matchp) 
    return true;
  else 
    return false;
  
}
#endif

  //#ifdef HAVE_AHOCORASICK_AHOCORASICK_H
#ifndef __APPLE__
// check if a string contains a substring using Aho Corasick algorithm
bool AbstractRule::ahomatch(const char * seq, unsigned len) {

  // make into Ac strcut
  AC_TEXT_t tmp_text = {seq, len}; //, static_cast<unsigned>(seq.length())};
  ac_automata_settext (atm, &tmp_text, 0);

  // do the check
  AC_MATCH_t * matchp;  
  matchp = ac_automata_findnext(atm);

  if (matchp) 
    return true;
  else 
    return false;
}
#endif

  void AbstractRule::parseSeqLine(const Json::Value& value) {
    Json::Value null(Json::nullValue);
    if (value.get("motif", null) != null) 
      atm_file = value.get("motif", null).asString();
    else
      return;

#ifdef __APPLE__
    std::cerr << "NOT AVAILBLE ON APPLE -- You are attempting to perform motif matching without Aho-Corasick library. Need to link to lahocorasick to do this." << std::endl;
    exit(EXIT_FAILURE);
#endif


  // open the sequence file
  igzstream iss(atm_file.c_str());
  if (!iss || !read_access_test(atm_file)) {
    std::cerr << "ERROR: Cannot read the sequence file: " << atm_file << std::endl;
    exit(EXIT_FAILURE);
  }

  // should it be inverted?
  std::string inv;
  //if (line.at(0) == '!') {
  //  atm_inv = true;
  //  inv = " -- Inverted -- ";
  //}

//#ifdef HAVE_AHOCORASICK_AHOCORASICK_H
#ifndef __APPLE__
  // initialize it
  if (!atm)
    atm = ac_automata_init(); //atm_ptr(ac_automata_init(), atm_free_delete);
  // make the Aho-Corasick key
  std::cerr << "...generating Aho-Corasick key"  << inv << " from file " << atm_file << std::endl;
  std::string pat;
  size_t count = 0;
  while (getline(iss, pat, '\n')) {
    ++count;
    std::istringstream iss2(pat);
    std::string val;
    size_t count2 = 0;
    /// only look at second element
    while (getline(iss2, val, '\t')) {
	++count2;
	if (count2 > 1)
	  continue;
	AC_PATTERN_t tmp_pattern;
	tmp_pattern.astring = val.c_str();
	tmp_pattern.length = static_cast<unsigned>(val.length());
	ac_automata_add(atm, &tmp_pattern);
      }
  }
  //ac_automata_finalize(atm);
  //std::cerr << "Done generating Aho-Corasick key of size " << count << std::endl;  

  atm_count = count;
#endif

  return;

  }


  void AbstractRule::parseSubLine(const Json::Value& value) {
    Json::Value null(Json::nullValue);
    if (value.get("sub", null) != null) 
      subsam_frac = value.get("sub", null).asDouble();
  }
  
GRC MiniRulesCollection::getAllRegions() const
{
  GRC out;

  for (auto& i : m_regions)
    out.concat(i.m_grv);

  out.mergeOverlappingIntervals();
  return out;
}

  void MiniRulesCollection::countsToFile(const std::string& file) const
  {
    std::ofstream of;
    of.open(file);
    char sep = '\t';
    of << "total_seen_count" << sep << "total_passed_count" << sep << "region" << sep << "region_passed_count" << sep << "rule" << sep << "rule_passed_count" << std::endl;
    for (auto& g : m_regions)
      for (auto& r : g.m_abstract_rules)
	of << m_count_seen << sep << m_count << sep << g.id << sep << g.m_count << sep << r.id << sep << r.m_count << std::endl;
    of.close();
    
  }

const std::string MiniRulesCollection::GetScriptContents(const std::string& script) {
  
  std::ifstream iss_file(script.c_str());

  std::string output;
  std::string line;
  while(getline(iss_file, line, '\n')) {
    output += line;
  }

  return(output);

}

}

