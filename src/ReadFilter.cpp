#include "SeqLib/ReadFilter.h"

#include <cassert>
#include "htslib/khash.h"
#include <unordered_set>

//#define QNAME "D0EN0ACXX111207:7:2306:6903:136511"
//#define QFLAG -1

//#define DEBUG_MINI 1

namespace SeqLib {

  // return if this rule accepts all reads
  bool AbstractRule::isEvery() const {
    return read_group.empty() && ins.isEvery() && del.isEvery() && isize.isEvery() && mapq.isEvery() && len.isEvery() && clip.isEvery() && phred.isEvery() && nm.isEvery() && nbases.isEvery() && fr.isEvery() && 
      (subsam_frac >= 1) && xp.isEvery();
  }


// define what is a valid condition
static const std::unordered_set<std::string> valid = 
  { 
  "duplicate", "supplementary", "qcfail", "hardclip", "fwd_strand",
  "rev_strand", "mate_fwd_strand", "mate_rev_strand", "mapped",
  "mate_mapped", "isize","clip", "phred", "length","nm",
  "mapq", "all", "ff", "xp","fr","rr","rf",
  "ic", "discordant","motif","nbases","!motif","allflag", "!allflag", "anyflag", "!anyflag",
  "ins","del",  "sub",  "subsample", "rg"
};

  static const std::unordered_set<std::string> allowed_region_annots = 
    { "region","pad", "matelink", "exclude", "rules"};

  static const std::unordered_set<std::string> allowed_flag_annots = 
    {"duplicate", "supplementary", "qcfail", "hardclip", 
     "fwd_strand", "rev_strand", "mate_fwd_strand", "mate_rev_strand",
     "mapped", "mate_mapped", "ff", "fr", "rr", "rf", "ic"};

bool ReadFilter::isValid(BamRecord &r) {

  // empty default is pass
  if (!m_abstract_rules.size())
    return true;

  for (auto& it : m_abstract_rules)
    if (it.isValid(r)) {
      ++it.m_count; //update this rule counter
      ++m_count;
       return true; // it is includable in at least one. 
    }
      
  return false;

}

  int FlagRule::__parse_json_int(const Json::Value& v) {

      try {
	if (v.asInt())
	 return v.asInt();
	else if (v.isString())
	  return std::stoi(v.asString());
      } catch (...) {
	std::cerr << " trouble converting flag to int on " << v << std::endl;
      }

      return 0;
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


  bool ReadFilterCollection::__validate_json_value(const Json::Value value, const std::unordered_set<std::string>& valid_vals) {

    // verify that it has appropriate values
    for (auto& i : value.getMemberNames()) {
      if (!valid_vals.count(i)) {
	std::cerr << "Invalid key value in JSON: " << i << std::endl;
	return false;
      }
    }

    return true;
  }

// check whether a BamAlignment (or optionally it's mate) is overlapping the regions
// contained in these rules
bool ReadFilter::isReadOverlappingRegion(BamRecord &r) {

  // if this is a whole genome rule, it overlaps
  if (!m_grv.size()) 
    return true;

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
  bool ReadFilterCollection::isValid(BamRecord &r) {

  ++m_count_seen;

#ifdef QNAME 
  if (r.Qname() == QNAME)
    std::cerr << "...starting valid checking in RFC. Regions size : "  << m_regions.size() << std::endl;
#endif

  if (m_regions.size() == 0)
    return true;


#ifdef QNAME 
  if (r.Qname() == QNAME)
    std::cerr << "...starting valid checking in array of filters in RFC" << std::endl;
#endif
  
  // need to run all rules if there is an excluder
  if (!m_fall_through)
    for (auto& i : m_regions)
      m_fall_through = m_fall_through || i.excluder;

  bool is_valid = false;
  bool exclude_hit = false; // did we hit excluder rule

  for (auto& it : m_regions) {

    // only check read validity if it overlaps region
    if (!it.isReadOverlappingRegion(r)) 
      continue;
    
    // check the region with all its rules
    if (it.isValid(r)) {
     
      // if this is excluder region, exclude read
      if (it.excluder)
	exclude_hit = true;
      
      // if we don't need to check all regions, quit here
      if (!m_fall_through) {
	++m_count;
	return true;
      }

      // in case we do fall through, track that we passed here
      is_valid = true;

    }
  }

    // found a hit in a rule
    if (!exclude_hit && is_valid) {
      ++m_count;
      return true;
    }
    
    return false;
}

  void ReadFilter::AddRule(const AbstractRule& ar) {
    m_abstract_rules.push_back(ar);
  }

  // constructor to make a ReadFilterCollection from a rules file.
  // This will reduce each individual BED file and make the 
  // GenomicIntervalTreeMap
  ReadFilterCollection::ReadFilterCollection(const std::string& script, const BamHeader& hdr) {

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

    // meake sure it at least has a rule
    /*std::stringstream buffer;
    if (read_access_test(script)) {
      std::ifstream t(script);
      buffer << t.rdbuf();
    } else {
      buffer << script;
    }
    if (buffer.str().find("\"rules\"") == std::string::npos) {
      std::cerr << " !!! JSON must be formated as {\"region\" : \{\"rules\" : [{...}]}}" << std::endl 
		<< " where \"rules\" is a keyword " << std::endl;
      exit(EXIT_FAILURE);
      }*/
    
    Json::Value null(Json::nullValue);
    
    int level = 1;

    // assign the global rule if there is one
    // remove from the rest of the rules
    Json::Value glob = root.removeMember("global");
    if (!glob.isNull()) 
      rule_all.parseJson(glob);
    
    // iterator over regions
    for (auto& regions : root) {
     
      if (!__validate_json_value(regions, allowed_region_annots))
	exit(EXIT_FAILURE);

      ReadFilter mr;
      
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
      if (v != null) {
	reg = v.asString();
	mr.id = mr.id + reg;
      }
      
      // actually parse the region
      if (reg == "WG" || reg.empty())
	mr.m_grv.clear(); // ensure it is whole-genome
      else
	mr.setRegions(GRC(reg, hdr));
	// debug mr.setRegionFromFile(reg, hdr);
      
      // check if its excluder region
      mr.excluder = false; // default is no exclude
      v = regions.get("exclude", null);
      if (v != null) {
	mr.excluder = v.asBool();
	if (mr.excluder)
	  mr.id = mr.id + "_exclude";
      }
      
      // set the rules
      v = regions.get("rules", null);
      if (!v.size()) {
	//std::cerr << " !!!! RULES size is zero. !!!! " << std::endl;
	//exit(EXIT_FAILURE);
      }

      // loop through the rules
      for (auto& vv : v) {
	if (vv != null) {
	  if (!__validate_json_value(vv, valid))
	    exit(EXIT_FAILURE);
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
      ReadFilter mr;
      mr.m_abstract_rules.push_back(rule_all);
      mr.id = "WG_includer";
      m_regions.push_back(mr);
    }
    
  }
  
  void ReadFilter::setRegions(const GRC& g) {
    m_grv = g;
    m_grv.createTreeMap();
  }

  void ReadFilter::addRegions(const GRC& g) {
    m_grv.concat(g);
    m_grv.mergeOverlappingIntervals();
    m_grv.createTreeMap();
  }


  // print the ReadFilterCollection
  std::ostream& operator<<(std::ostream &out, const ReadFilterCollection &mr) {
    
    out << "----------ReadFilterCollection-------------" << std::endl;
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

// print a ReadFilter information
std::ostream& operator<<(std::ostream &out, const ReadFilter &mr) {
  
  std::string file_print = !mr.m_grv.size() ? "WHOLE GENOME" : mr.m_region_file;
  out << (mr.excluder ? "--Exclude Region: " : "--Include Region: ") << file_print;
  if (mr.m_grv.size()) {
    //out << " --Size: " << AddCommas<int>(mr.m_width); 
    out << " Pad: " << mr.pad;
    out << " Matelink: " << (mr.m_applies_to_mate ? "ON" : "OFF");
    if (mr.m_grv.size() == 1)
      out << " Region : " << mr.m_grv[0] << std::endl;
    else
      out << " " << mr.m_grv.size() << " regions" << std::endl;      
  } else {
    out << std::endl;
  }

  for (auto& it : mr.m_abstract_rules) 
    out << it << std::endl;
  return out;
}

  void ReadFilterCollection::AddReadFilter(const ReadFilter& rf) {
    m_regions.push_back(rf);
  }

  ReadFilter::~ReadFilter() {}

  bool Flag::parseJson(const Json::Value& value, const std::string& name) {

    if (value.isMember(name.c_str())) {
      __convert_to_bool(value, name) ? setOn() : setOff();
      return true;
    }
    
    return false;

  }

  void FlagRule::parseJson(const Json::Value& value) {

    Json::Value null(Json::nullValue);
    if (value.isMember("allflag"))
      setAllOnFlag(__parse_json_int(value.get("allflag", null)));
    if (value.isMember("!allflag"))
      setAllOffFlag(__parse_json_int(value.get("!allflag", null)));
    if (value.isMember("anyflag"))
      setAnyOnFlag(__parse_json_int(value.get("anyflag", null)));
    if (value.isMember("!anyflag"))
      setAnyOffFlag(m_any_off_flag = __parse_json_int(value.get("!anyflag", null)));
    
    // have to set the every if find flag so that rule knows it cant skip checking
    if (dup.parseJson(value, "duplicate")) every = false;
    if (supp.parseJson(value, "supplementary")) every = false;
    if (qcfail.parseJson(value, "qcfail")) every = false;
    if (hardclip.parseJson(value, "hardclip")) every = false;
    if (fwd_strand.parseJson(value, "fwd_strand")) every = false;
    if (mate_rev_strand.parseJson(value, "mate_rev")) every = false;
    if (mate_fwd_strand.parseJson(value, "mate_fwd")) every = false;
    if (mate_mapped.parseJson(value, "mate_mapped")) every = false;
    if (mapped.parseJson(value, "mapped")) every = false;
    if (ff.parseJson(value, "ff")) every = false;
    if (fr.parseJson(value, "fr")) every = false;
    if (rf.parseJson(value, "rf")) every = false;
    if (rr.parseJson(value, "rr")) every = false;
    if (ic.parseJson(value, "ic")) every = false;

  }
  
  void Range::parseJson(const Json::Value& value, const std::string& name) {
    Json::Value null(Json::nullValue);
    Json::Value v = value.get(name, null);

    if (v != null) {
      if (v.size() > 2) {
	std::cerr << " ERROR. Not expecting array size " << v.size() << " for Range " << name << std::endl;
      } else {
	m_every = false;
	m_inverted = false;

	if (v.isArray()) {
	  m_min = v[0].asInt();
	  m_max = v[1].asInt();
	} else if (v.isInt()) {
	  m_min = v.asInt();
	  m_max = INT_MAX;
	} else if (v.isBool()) {
	  m_min = v.asBool() ? 1 : INT_MAX; // if true, [1,MAX], if false [MAX,1] (not 1-MAX)
	  m_max = v.asBool() ? INT_MAX : 1;
	} else {
	  std::cerr << "Unexpected type for range flag: " << name << std::endl;
	  exit(EXIT_FAILURE);
	}

	if (m_min > m_max) {
	  m_inverted = true;
	  std::swap(m_min, m_max); // make min always lower
	}
      }
	
    }
  }

  void AbstractRule::parseJson(const Json::Value& value) {

    // parse read group
    const std::string rg = "rg";
    if (value.isMember(rg.c_str())) {
      Json::Value null(Json::nullValue);
      Json::Value v = value.get(rg, null);
      assert(v != null);
      read_group = v.asString();
    }
      
    // set the ID
    for (auto& i : value.getMemberNames()) {
      id += i + ";";
    }
    if (id.length())
      id.pop_back();

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
    
  }


// main function for determining if a read is valid
  bool AbstractRule::isValid(BamRecord &r) {
    

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
    bool need_to_continue = !nm.isEvery() || !clip.isEvery() || !len.isEvery() || !nbases.isEvery() || 
      !xp.isEvery();

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
      r.QualityTrimmedSequence(phred.lowerBound(), startpoint, endpoint);
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
      if (!xp.isValid(r.CountBWASecondaryAlignments()))
	return false;
    
#ifdef QNAME
    if (r.Qname() == QNAME && (r.AlignmentFlag() == QFLAG || QFLAG == -1))
      std::cerr << "****** READ ACCEPTED ****** " << std::endl;
#endif


    return true;
  }
  
  bool FlagRule::isValid(BamRecord &r) {
    
#ifdef QNAME
    if (r.Qname() == QNAME && r.AlignmentFlag() == QFLAG)
      std::cerr << " IN FLAG RULE: EVERY? " << isEvery() << " ON NUMBER FLAG " << m_all_on_flag << " OFF NUMEBR FLAG " << m_all_off_flag <<  " ON NUM RESULT "  << std::endl;
#endif

    if (isEvery())
      return true;

    // if have on or off flag, use that
    // 0001100 - all flag
    // 0101000 - flag
    // -------
    // 0001000 - should fail all flag. Should pass any flag
    if (m_all_on_flag && !( (r.AlignmentFlag() & m_all_on_flag)  == m_all_on_flag) ) // if all on, pass
      return false;
    if (m_all_off_flag && ( (r.AlignmentFlag() & m_all_off_flag) == m_all_off_flag) ) // if all on, fail
      return false;

    // if have on or off flag, use that
    if (m_any_on_flag && !(r.AlignmentFlag() & m_any_on_flag) ) // if ANY on, pass
      return false;
    if (m_any_off_flag && (r.AlignmentFlag() & m_any_off_flag)) // if ANY on, fail
      return false;
       
#ifdef QNAME
    if (r.Qname() == QNAME) // && r.AlignmentFlag() == QFLAG)
      std::cerr << " CHECKING NAMED FLAGS " << std::endl;
#endif

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

  out << "  Rule: ";
  if (ar.isEvery()) {
    out << "  ALL";
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

  if (fr.m_all_on_flag)
    keep += "[(all)" + std::to_string(fr.m_all_on_flag) + "],";
  if (fr.m_all_off_flag)
    remo += "[(all)" + std::to_string(fr.m_all_off_flag) + "],";

  if (fr.m_any_on_flag)
    keep += "[(any)" + std::to_string(fr.m_any_on_flag) + "],";
  if (fr.m_any_off_flag)
    remo += "[(any)" + std::to_string(fr.m_any_off_flag) + "],";

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
    out << "ALL";
  else
    out << (r.m_inverted ? "NOT " : "") << "[" << r.m_min << "," << (r.m_max == INT_MAX ? "MAX" : std::to_string(r.m_max))  << "]";
  return out;
}

#ifdef HAVE_AHO_CORASICK
#ifndef __APPLE__
// check if a string contains a substring using Aho Corasick algorithm
//bool AbstractRule::ahomatch(const string& seq) {
bool AbstractRule::ahomatch(BamRecord &r) {

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
#endif

#ifdef HAVE_AHO_CORASICK
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
#endif

#ifdef HAVE_AHO_CORASICK
  void AbstractRule::parseSeqLine(const Json::Value& value) {
    
    bool i = false; // invert motif?
    std::string motif_file;
    Json::Value null(Json::nullValue);
    if (value.get("motif", null) != null) 
      motif_file = value.get("motif", null).asString();
    else if (value.get("!motif", null) != null) {
      motif_file = value.get("!motif", null).asString();
      i = true;
    }
    else
      return;

#ifdef __APPLE__
    std::cerr << "NOT AVAILBLE ON APPLE -- You are attempting to perform motif matching without Aho-Corasick library. Need to link to lahocorasick to do this." << std::endl;
    exit(EXIT_FAILURE);
#endif

    addMotifRule(motif_file, i);

  return;

  }
#endif

#ifdef AHO_CORASICK
  void AbstractRule::addMotifRule(const std::string& f, bool inverted) {
    std::cerr << "...making the AhoCorasick trie from " << f << std::endl;
    aho.TrieFromFile(f);
    std::cerr << "...finished making AhoCorasick trie with " << AddCommas(aho.count) << " motifs" << std::endl;
    aho.inv = inverted;
  }
#endif  

  void AhoCorasick::TrieFromFile(const std::string& f) {

    file = f;

    // open the sequence file
    std::ifstream iss(f.c_str());
    if (!iss || !read_access_test(f)) 
      throw std::runtime_error("AhoCorasick::TrieFromFile - Cannot read file: " + f);
    
    // make the Aho-Corasick trie
    std::string pat;
    while (getline(iss, pat, '\n')) {
      ++count;
      AddMotif(pat);
    }
  }
  
  void AbstractRule::parseSubLine(const Json::Value& value) {
    Json::Value null(Json::nullValue);
    if (value.get("sub", null) != null) 
      subsam_frac = value.get("sub", null).asDouble();
  }
  
GRC ReadFilterCollection::getAllRegions() const
{
  GRC out;

  for (auto& i : m_regions)
    out.concat(i.m_grv);

  return out;
}

  std::string ReadFilterCollection::EmitCounts() const
  {
    std::stringstream of;
    char sep = '\t';
    of << "total_seen_count" << sep << "total_passed_count" << sep << "region" << sep << "region_passed_count" << sep << "rule" << sep << "rule_passed_count" << std::endl;
    for (auto& g : m_regions)
      for (auto& r : g.m_abstract_rules)
	of << m_count_seen << sep << m_count << sep << g.id << sep << g.m_count << sep << r.id << sep << r.m_count << std::endl;

    return of.str();
    
  }

  ReadFilter::ReadFilter() {}

  /*  ReadFilter::ReadFilter(const CommandLineRegion& c, const SeqLib::BamHeader& hdr) {

    m_region_file = c.f;

    // set a whole genome ALL rule
    if (c.type < 0) {
      //m_grv.clear();
      id = "WG";
    } else {
      // set the genomic region this rule applies to
      setRegions(GRC(c.f, hdr));
      //debug setRegionFromFile(c.f, hdr);
    }

    // add the abstract rule
    AbstractRule ar;

    // set the flag
    if (c.i_flag || c.e_flag) {
      ar.fr.setAllOnFlag(c.i_flag);
      ar.fr.setAllOffFlag(c.e_flag);
    }

    // set the other fields
    if (c.len)
      ar.len = Range(c.len, INT_MAX, false);
    if (c.nbases != INT_MAX)
      ar.nbases = Range(0, c.nbases, false);
    if (c.phred)
      ar.phred = Range(c.phred, INT_MAX, false);
    if (c.mapq)
      ar.mapq = Range(c.mapq, INT_MAX, false);
    if (c.clip)
      ar.clip = Range(c.clip, INT_MAX, false);
    if (c.del)
      ar.del = Range(c.del, INT_MAX, false);
    if (c.ins)
      ar.ins = Range(c.ins, INT_MAX, false);

    // set the id
    ar.id = id + "_CMD_RULE";

#ifdef HAVE_AHO_CORASICK
    // add a motif rule
    if (!c.motif.empty())
      ar.addMotifRule(c.motif, false);
#endif

    // add read group rule
    ar.read_group = c.rg;

    m_abstract_rules.push_back(ar);

    // set the properties of the region
    if (c.type >= 0) {
      switch(c.type) {
      case MINIRULES_MATE_LINKED:
	m_applies_to_mate = true;
	excluder = false;
	break;
      case MINIRULES_MATE_LINKED_EXCLUDE:
	m_applies_to_mate = true;
	excluder = true;
	break;
      case MINIRULES_REGION:
	m_applies_to_mate = false;
	excluder = false;
	break;
      case MINIRULES_REGION_EXCLUDE:
	m_applies_to_mate = false;
	excluder = true;
	break;
      default:
	std::cerr << "Unexpected type in ReadFilter. Exiting" << std::endl;
	exit(EXIT_FAILURE);
      }
    }
    
  }
  */

}

