#include "SnowTools/MiniRules.h"

#include <cassert>
#include "htslib/khash.h"

//#define DEBUG_GIVEN_READ
//#define QNAME "HFCNNADXX150330:2:1203:21337:34241"
//#define QFLAG 99

//#define DEBUG_MINI 1

namespace SnowTools {

// define what is a valid condition
static const std::unordered_map<std::string,bool> valid = 
  { 
  {"duplicate",     true},
  {"supplementary", true},
  {"qcfail",        true},
  {"hardclip",      true},
  {"fwd_strand",    true},
  {"rev_strand",    true},
  {"mate_fwd_strand",  true},
  {"mate_rev_strand",  true},
  {"mapped",           true},
  {"mate_mapped",      true},
  {"isize", true},
  {"clip",  true},
  {"phred", true},
  {"length",   true},
  {"nm",    true},
  {"mapq",  true},
  {"all",   true},
  {"ff", true},
  {"xp", true},
  {"fr", true},
  {"rr", true},
  {"rf", true},
  {"ic", true},
  {"discordant", true},
  {"motif", true},
  {"nbases", true},
  {"ins", true},
  {"del", true},
  {"sub", true},
  {"subsample", true}
};


bool MiniRules::isValid(BamRead &r) {

  for (auto& it : m_abstract_rules)
    if (it.isValid(r)) 
       return true; // it is includable in at least one. 
      
  return false;

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

  size_t which_region = 0;
  size_t which_rule = 0;
  
  // find out which rule it is a part of
  // lower number rules dominate

  bool is_valid = false;

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
	  // this whole read is valid
	  is_valid = true;

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
    if (rule_hit && !m_fall_through)
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
void MiniRulesCollection::__construct_MRC(const std::string& file) {

  // parse the rules file
  std::vector<std::string> region_files;
  std::ifstream iss_file(file.c_str());
  char delim = '%';

  std::string ffile = file;;

  if (iss_file) {
    std::string temp;
    ffile = "";
    while(getline(iss_file, temp)) {
      ffile += temp + "\n";
    }
    iss_file.close();
    delim = '\n';
  }
  std::istringstream iss_rules(ffile.c_str());
  
  // loop through the rules file and grab the rules
  std::string line;
  int level = 1;

  // define a default rule set
  std::vector<AbstractRule> all_rules;

  while(getline(iss_rules, line, delim)) {

#ifdef DEBUG_MINI
    std::cerr << line << std::endl;
#endif
    //exclude comments and empty lines
    bool line_empty = line.find_first_not_of("\t\n ") == std::string::npos;
    bool line_comment = false;
    if (!line_empty)
      line_comment = line.at(0) == '#';
    
    if (line_comment || line_empty)
      continue;
    
    // check that it doesn't have too many rules on it
    if (count(line.begin(), line.end(), '@') > 1) {
      std::cerr << "ERROR: Every line must start with region@, global@ or specify a rule, and only one region/rule per line" << std::endl;
      std::cerr << "  If separating lines in -r flag, separate with %. If in file, use \\n" << std::endl;
      std::cerr << "  Offending line: " << line << std::endl;
      exit(EXIT_FAILURE);
    }
    
    //////////////////////////////////
    // its a rule line, get the region
    //////////////////////////////////
    if (line.find("region@") != std::string::npos) {
      
      // check that the last one isnt empty. 
      // if it is, add the global to it
      //if (m_regions.size() > 0)
      //	if (m_regions.back().m_abstract_rules.size() == 0) {
      //  m_regions.back().m_abstract_rules.push_back(rule_all);
      //	}
      
      // start a new MiniRule set
      MiniRules mr;
      mr.mrc = this; // set the pointer to the collection
      
      // add the defaults
      //mr->m_abstract_rules = all_rules;
      
      // check if the mate aplies
      if (line.find(";mate") != std::string::npos || line.find("mlregion") != std::string::npos) {
	mr.m_applies_to_mate = true;
      }

      // check if we should pad 
      std::regex reg_pad(".*?pad\\[([0-9]+)\\].*"); 
      
      std::smatch pmatch;
      if (std::regex_search(line, pmatch, reg_pad))
      	try { mr.pad = std::stoi(pmatch[1].str()); } catch (...) { std::cerr << "Cant read pad value for line " << line << ", setting to 0" << std::endl; }

      if (line.find("@WG") != std::string::npos) {
	  mr.m_whole_genome = true;
      } else {
	std::regex file_reg(".*?region@(.*?)(;|$)");
	std::smatch match;
	if (std::regex_search(line, match, file_reg))
	  mr.setRegionFromFile(match[1].str());
	else {
	  std::cerr << "Could not parse line: " << line << " to grab region " << std::endl;
	  exit(EXIT_FAILURE);
	}
      }

      mr.m_level = level++;
      mr.id = line;
      m_regions.push_back(mr);

    }
    ////////////////////////////////////
    // its a global rule
    ///////////////////////////////////
    else if (line.find("global@") != std::string::npos) {
      rule_all.parseRuleLine(line);
    }
    ////////////////////////////////////
    // its a rule
    ////////////////////////////////////
    else {

      AbstractRule ar = rule_all; // always start with the global rule
      // parse the line
      ar.parseRuleLine(line);

      m_regions.back().m_abstract_rules.push_back(ar);      
      m_regions.back().parseDiscordantShortcut(line, ar);

#ifdef DEBUG_MINI
      std::cerr << "trying to make rule for " << line << std::endl;
      std::cerr << " rule is " << ar << std::endl;
#endif
    }
    
  } // end \n parse
  
  // check that the last region has at least one reuls
  // if it it doesn't, give it the global
  if (m_regions.size() > 0)
    if (m_regions.back().m_abstract_rules.size() == 0) {
      m_regions.back().m_abstract_rules.push_back(rule_all);
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
  out << "--Region: " << file_print;
  if (!mr.m_whole_genome) {
    //out << " --Size: " << AddCommas<int>(mr.m_width); 
    out << " --Pad: " << mr.pad;
    out << " --Include Mate: " << (mr.m_applies_to_mate ? "ON" : "OFF") << std::endl;
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

// parse a rule line looking for flag values
void FlagRule::parseRuleLine(std::string line) {

  std::istringstream iss(line);
  std::string val;
  while (getline(iss, val, ';')) {
    std::regex reg_dup("^!?dup.*");
    std::regex reg_sup("^!?supp.*");
    std::regex reg_qc("^!?qcfail$");
    std::regex reg_fs("^!?fwd_strand$");
    std::regex reg_hc("^!?hardclip.*");
    std::regex reg_rs("^!?rev_strand$");
    std::regex reg_mf("^!?mate_fwd_strand$");
    std::regex reg_mr("^!?mate_rev_strand$");
    std::regex reg_mp("^!?mapped$");
    std::regex reg_mm("^!?mate_mapped$");
    std::regex reg_ff("^!?ff$");
    std::regex reg_fr("^!?fr$");
    std::regex reg_rf("^!?rf$");
    std::regex reg_rr("^!?rr$");
    std::regex reg_ic("^!?ic$");

    if (dup.parseRuleLine(val, reg_dup))   na = false;
    if (supp.parseRuleLine(val, reg_sup))  na = false;
    if (qcfail.parseRuleLine(val, reg_qc)) na = false;
    if (hardclip.parseRuleLine(val, reg_hc))        na = false;
    if (fwd_strand.parseRuleLine(val, reg_fs))      na = false;
    if (mate_rev_strand.parseRuleLine(val, reg_mr)) na = false;
    if (mate_fwd_strand.parseRuleLine(val, reg_mf)) na = false;
    if (mate_mapped.parseRuleLine(val, reg_mm))     na = false;
    if (mapped.parseRuleLine(val, reg_mp)) na = false;
    if (ff.parseRuleLine(val, reg_ff))     na = false;
    if (fr.parseRuleLine(val, reg_fr))     na = false;
    if (rf.parseRuleLine(val, reg_rf))     na = false;
    if (rr.parseRuleLine(val, reg_rr))     na = false;
    if (ic.parseRuleLine(val, reg_ic))     na = false;
  }

}

// modify the rules based on the informaiton provided in the line
void AbstractRule::parseRuleLine(std::string line) {

  id += line + ";";

  // get everything but the global keyword, if there is one
  std::regex reg_noname("global@(.*)");
  std::smatch nnmatch;
  std::string noname;
  if (std::regex_search(line, nnmatch, reg_noname)) {
    noname = nnmatch[1].str();
  } else {
    noname = line;
  }

  // check that the conditoins are valid
  std::istringstream iss_c(noname);
  std::string tmp;

  while (getline(iss_c, tmp, ';')) {
    std::regex reg(".*?!?([a-z_]+).*");
    std::smatch cmatch;

    if (std::regex_search(tmp, cmatch, reg)) {
      if (valid.count(cmatch[1].str()) == 0) {
	std::cerr << "Invalid condition of: " << tmp << std::endl;
	exit(EXIT_FAILURE);
      }
    } else {
      std::cerr << "Invalid condition of: " << tmp << std::endl;
      exit(EXIT_FAILURE);
    }
  }
  
  // check for every/none flags
  if (noname.find("all") != std::string::npos) 
    setEvery();
  if (noname.find("!all") != std::string::npos)
    setNone();

  // modify the ranges if need to
  if (noname.find("isize") != std::string::npos)
    isize.parseRuleLine(noname);
  if (noname.find("mapq") != std::string::npos)
    mapq.parseRuleLine(noname);
  if (noname.find("len") != std::string::npos)
    len.parseRuleLine(noname);
  if (noname.find("clip") != std::string::npos)
    clip.parseRuleLine(noname);
  if (noname.find("phred") != std::string::npos)
    phred.parseRuleLine(noname);
  if (noname.find("nbases") != std::string::npos)
    nbases.parseRuleLine(noname);
  if (noname.find("ins") != std::string::npos)
    ins.parseRuleLine(noname);
  if (noname.find("del") != std::string::npos)
    del.parseRuleLine(noname);
  if (noname.find("nm") != std::string::npos)
    nm.parseRuleLine(noname);
  if (noname.find("xp") != std::string::npos)
    xp.parseRuleLine(noname);


  // parse the subsample data
  parseSubLine(noname);

  // parse the line for flag rules (also checks syntax)
  fr.parseRuleLine(noname);

  // parse aho corasick file, if not already inheretid
  if (!atm) {
    std::istringstream iss_m(noname);
    while (getline(iss_m, tmp, ';')) 
      if (tmp.find("motif") != std::string::npos)
	parseSeqLine(tmp);
    if (atm) {
      ac_automata_finalize(atm);
      std::cerr << "Done generating Aho-Corasick tree" << std::endl;  
    }
  }
    
}

// parse for range
void Range::parseRuleLine(std::string line) {
  
  std::istringstream iss(line);
  std::string val;
  
  while (getline(iss, val, ';')) {
    
    std::string i_reg_str = "!" + pattern + ":?\\[(.*?),(.*?)\\]";
    std::string   reg_str = pattern + ":?\\[(.*?),(.*?)\\]";
    
    std::string n_reg_str = pattern + ":?!all";
    std::string a_reg_str = pattern + ":?all";
    
    std::regex ireg(i_reg_str);
    std::regex  reg(reg_str);
    std::regex nreg(n_reg_str);
    std::regex areg(a_reg_str);
    
    std::smatch match;
    if (std::regex_search(val, match, areg)) {
      setEvery();
    } else if (std::regex_search(val, match, ireg)) {
      try {
	min = std::stoi(match[1].str());
	max = std::stoi(match[2].str());
	inverted = true;
	every = false; none = false;
	return;
      } catch (...) {
	std::cerr << "Caught error trying to parse inverted for " << pattern << " on line " << line << " match[1] " << match[1].str() << " match[2] " << match[2].str() << std::endl;     
	exit(EXIT_FAILURE);
      }
    } else if (std::regex_search(val, match, reg)) {
      try {
	min = std::stoi(match[1].str());
	max = std::stoi(match[2].str());
	inverted = false;
	every = false; none = false;
	return;
      } catch (...) {
	std::cerr << "Caught error trying to parse for " << pattern << " on line " << line << " match[1] " << match[1].str() << " match[2] " << match[2].str() << std::endl;     
	exit(EXIT_FAILURE);
      }

    }
    //else {
    //  std::cerr << "Could not parse rule: " << val << " for range " << pattern << std::endl;
    //  exit(EXIT_FAILURE);
    //}
    
  } // end getline
}

// main function for determining if a read is valid
  bool AbstractRule::isValid(BamRead &r) {
    

#ifdef DEBUG_GIVEN_READ
    if (r.Qname() == QNAME && r.AlignmentFlag() == QFLAG)
      std::cerr << "MINIRULES Read seen " << " ID " << id << " " << r << std::endl;
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
    bool isize_pass = isize.isValid(abs(r.InsertSize()));

#ifdef DEBUG_GIVEN_READ
    if (r.Qname() == QNAME && r.AlignmentFlag() == QFLAG)
      std::cerr << "MINIRULES isize_pass " << isize_pass << " " << " ID " << id << " " << r << std::endl;
#endif
    
    if (!isize_pass) {
      return false;
    }
    
    // check for valid mapping quality
    if (!mapq.isEvery())
      if (!mapq.isValid(r.MapQuality())) 
	return false;

#ifdef DEBUG_GIVEN_READ
    if (r.Qname() == QNAME && r.AlignmentFlag() == QFLAG)
      std::cerr << "MINIRULES mapq pass " << " ID " << id << " " << r << std::endl;
#endif
    
    // check for valid flags
    if (!fr.isValid(r))
      return false;

#ifdef DEBUG_GIVEN_READ
    if (r.Qname() == QNAME && r.AlignmentFlag() == QFLAG)
      std::cerr << "MINIRULES flag pass " << " " << id << " " << r << std::endl;
#endif
    
    // check the CIGAR
    if (!ins.isEvery() || !del.isEvery()) {
      if (!ins.isValid(r.MaxInsertionBases()))
	return false;
      if (!del.isValid(r.MaxDeletionBases()))
	return false;
    }


#ifdef DEBUG_GIVEN_READ
    if (r.Qname() == QNAME && r.AlignmentFlag() == QFLAG)
      std::cerr << "MINIRULES cigar pass " << " ID " << id << " "  << r << std::endl;
#endif

    
    // if we dont need to because everything is pass, just just pass it
    bool need_to_continue = !nm.isEvery() || !clip.isEvery() || !len.isEvery() || !nbases.isEvery() || atm_file.length() || !xp.isEvery();
    if (!need_to_continue)
      return true;

#ifdef DEBUG_GIVEN_READ
    if (r.Qname() == QNAME  && r.AlignmentFlag() == QFLAG)
      std::cerr << "MINIRULES moving on. ID " << id << " " << r << std::endl;
#endif
    
    // now check if we need to build char if all we want is clip
    unsigned clipnum = 0;
    if (!clip.isEvery()) {
      clipnum = r.NumClip();
      if (nm.isEvery() && len.isEvery() && !clip.isValid(clipnum)) // if clip fails, its not going to get better by trimming. kill it now before building teh char data
	return false;
    }

#ifdef DEBUG_GIVEN_READ
    if (r.Qname() == QNAME && r.AlignmentFlag() == QFLAG)
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
      
      if (endpoint != -1 && new_len < r.Length() && new_len > 0 && new_len - startpoint >= 0 && startpoint + new_len < r.Length()) { 
	try { 
	  r.AddZTag("GV", r.Sequence().substr(startpoint, new_len));
	  assert(r.GetZTag("GV").length());
	} catch (...) {
	  std::cerr << "Subsequence failure with sequence of length "  
		    << r.Sequence().length() << " and startpoint "
		    << startpoint << " endpoint " << endpoint 
		    << " newlen " << new_len << std::endl;
	}
      }

#ifdef DEBUG_GIVEN_READ
    if (r.Qname() == QNAME && r.AlignmentFlag() == QFLAG)
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

#ifdef DEBUG_GIVEN_READ
    if (r.Qname() == QNAME && r.AlignmentFlag() == QFLAG)
      std::cerr << "MINIRULES NBASES PASS. id: " << id << r << " new_len " << new_len << " new_clipnum" << new_clipnum << std::endl;
#endif

    
    // check for valid length
    if (!len.isValid(new_len))
      return false;

#ifdef DEBUG_GIVEN_READ
    if (r.Qname() == QNAME && r.AlignmentFlag() == QFLAG)
      std::cerr << "MINIRULES LENGTH PASS. id: " << id << r << " newlen " << new_len << std::endl;
#endif
    
    // check for valid clip
    if (!clip.isValid(new_clipnum))
      return false;

#ifdef DEBUG_GIVEN_READ
    if (r.Qname() == QNAME && r.AlignmentFlag() == QFLAG)
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
    
#ifdef DEBUG_GIVEN_READ
    if (r.Qname() == QNAME && r.AlignmentFlag() == QFLAG)
      std::cerr << "MINIRULES PASS EVERYTHING. ID " << id << " " << r << std::endl;
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

    bool first = r.Position() < r.MatePosition();
    bool bfr = (first && (!r.ReverseFlag() && r.MateReverseFlag())) || (!first &&  r.ReverseFlag() && !r.MateReverseFlag());
    bool brr = r.ReverseFlag() && r.MateReverseFlag();
    bool brf = (first &&  (r.ReverseFlag() && !r.MateReverseFlag())) || (!first && !r.ReverseFlag() &&  r.MateReverseFlag());
    bool bff = !r.ReverseFlag() && !r.MateReverseFlag();
      
    bool bic = r.Interchromosomal();

    // its FR and it CANT be FR (off) or its !FR and it MUST be FR (ON)
    // orienation not defined for inter-chrom, so exclude these with !ic
    if (!bic) { // PROCEED IF INTRA-CHROMOSOMAL
      if ( (bfr && fr.isOff()) || (!bfr && fr.isOn())) 
	return false;
      // etc....
      if ( (brr && rr.isOff()) || (!brr && rr.isOn())) 
	return false;
      if ( (brf && rf.isOff()) || (!brf && rf.isOn())) 
	return false;
      if ( (bff && ff.isOff()) || (!bff && ff.isOn())) 
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
    out << "  KEEPING ALL";
  } else if (ar.isNone()) {
    out << "  KEEPING NONE";  
  } else {
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


  /*

  // get the strings
  for (auto it : fr.flags) {
    if (it.second.isNA())
      na += it.first + ",";
    else if (it.second.isOn())
      keep += it.first + ",";
    else if (it.second.isOff())
      remo += it.first + ",";
    else // shouldn't get here
      exit(1); 
  }
  */
  if (!fr.isEvery())
    out << keep << " -- " << remo;

  return out;
}

// define how to print
std::ostream& operator<<(std::ostream &out, const Range &r) {
  if (r.isEvery())
    out << "all";
  else
    out << (r.inverted ? "NOT " : "") << "[" << r.min << "," << r.max << "]";
  return out;
}

bool Flag::parseRuleLine(std::string &val, std::regex &reg) {

  std::smatch match;
  if (std::regex_search(val, match, reg)) {
    //auto ff = flags.find(match[1].str()); 
    if (val.at(0) == '!') { // it is a val in flags and is off
      setOff();
      return true;
    } else  { // is in a val in flags and is on
      setOn();
      return true;
    } //else if (ff == flags.end() && valid.count(match[1].str()) == 0) { // its not anything and its bad
      //std::cerr << "Not a valid condition: " << match[1].str() << " on val " << val << std::endl;
      //exit(EXIT_FAILURE);
    //}
  }

  return false;
  
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

// add the aho subsequences by reading in the sequence file
void AbstractRule::parseSeqLine(std::string line) {

  // get the sequence file out
  std::regex reg("^!?motif\\[(.*)\\].*");
  std::smatch match;
  if (std::regex_search(line, match, reg)) {
    line = match[1].str();
    atm_file = line;

    //#ifndef HAVE_AHOCORASICK_AHOCORASICK_H
#ifdef __APPLE__
    std::cerr << "NOT AVAILBLE ON APPLE -- Attempting to perform motif matching without Aho-Corasick library. Need to link to lahocorasick to do this." << std::endl;
    exit(EXIT_FAILURE);
#endif

  } else {
    return;
  }

  // open the sequence file
  igzstream iss(atm_file.c_str());
  if (!iss || !read_access_test(atm_file)) {
    std::cerr << "ERROR: Cannot read the sequence file: " << atm_file << std::endl;
    exit(EXIT_FAILURE);
  }

  // should it be inverted?
  std::string inv;
  if (line.at(0) == '!') {
    atm_inv = true;
    inv = " -- Inverted -- ";
  }

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

// parse the subsample line
void AbstractRule::parseSubLine(std::string line) {

  std::regex reg("^!?sub\\[(.*)\\].*");
  std::smatch match;
  if (std::regex_search(line, match, reg)) {
    try {
      subsam_frac = stod(match[1].str());
    } catch (...) {
      std::cerr << "ERROR parsing string for subsample. Line is: " << line << std::endl;
      exit(EXIT_FAILURE);
    }
  } else {
    return;
  }

  // check that it is valid
  assert(subsam_frac >= 0 && subsam_frac <= 1);

}

GRC MiniRulesCollection::getAllRegions() const
{
  GRC out;

  for (auto& i : m_regions)
    out.concat(i.m_grv);

  out.mergeOverlappingIntervals();
  return out;
}

  void MiniRules::parseDiscordantShortcut(const std::string& line, const AbstractRule& ar)
{

  // check for "discordant" shortcut
  std::regex  regex_disc( ".*?discordant\\[([0-9]+),([0-9]+)\\]($|;)");
  std::smatch omatch;
  if (std::regex_search(line, omatch, regex_disc)) {
    bool isneg = line.find("!discordant[") != std::string::npos;
    try {
      // fill in the isize condition 
      m_abstract_rules.back().isize.min = std::stoi(omatch[1].str());
      m_abstract_rules.back().isize.max = std::stoi(omatch[2].str());
      m_abstract_rules.back().isize.inverted = !isneg;
      m_abstract_rules.back().isize.every = false;
      m_abstract_rules.back().id += "INSERT_SIZE_CONDITION";
      // use the template to set a sequene of orientation rules
      AbstractRule aro = ar;
      if (isneg)
	aro.fr.ff.setOff();
      else
	aro.fr.ff.setOn();
      aro.fr.na = false;
      aro.id += "FF_CONDITION";
      m_abstract_rules.push_back(aro);
      // set another rr rule
      aro = ar;
      if (isneg)
	aro.fr.rr.setOff();
      else
	aro.fr.rr.setOn();
      aro.fr.na = false;
      aro.id += "RR_CONDITION";
      m_abstract_rules.push_back(aro);
      // set another rf rule
      aro = ar;
      if (isneg)
	aro.fr.rf.setOff();
      else
	aro.fr.rf.setOn();
      aro.fr.na = false;
      aro.id += "RF_CONDITION";
      m_abstract_rules.push_back(aro);
      // set another ic rule
      aro = ar;
      if (isneg)
	aro.fr.ic.setOff();
      else
	aro.fr.ic.setOn();
      aro.fr.na = false;
      aro.id += "INTERCHROM_CONDITION";
      m_abstract_rules.push_back(aro);
    } catch (...) {
      std::cerr << "Caught error trying to parse for discordant " << " on line " << line << " match[1] " << omatch[1].str() << " match[2] " << omatch[2].str() << std::endl;     
      exit(EXIT_FAILURE);
    }
  } // end discorant regex
  
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
}

