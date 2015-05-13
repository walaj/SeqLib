#ifndef MINI_RULES_H
#define MINI_RULES_H

#include <string>
#include <vector>
#include <unordered_map>
#include <regex>

#include "SnowTools/GenomicRegion.h"
#include "SnowTools/GenomicRegionCollection.h"
#include "SnowTools/SnowUtils.h"
#include "SnowTools/HTSTools.h"

#ifdef HAVE_AHOCORASICK_AHOCORASICK_H
#include "ahocorasick/ahocorasick.h"
// custom deleter for aho-corasick
struct atm_free_delete {
  void operator()(void* x) { free((AC_AUTOMATA_t*)x); }
};
typedef unique_ptr<AC_AUTOMATA_t> atm_ptr;
#endif

namespace SnowTools {

/** Stores a rule for a single alignment flag.
 *
 * Rules for alignment flags can be one of three states:
 * - NA - All flag values are valid
 * - Off - Flag is valid if OFF
 * - On - Flag is valid if ON
 */
class Flag {
  
 public:

  /** Construct a new Flag with NA rule 
   */
  Flag() : on(false), off(false), na(true) {}
    
    /** Set the flag to NA
     */
  void setNA() { on = false; off = false; na = true; } 
  
  /** Set the flag to ON
   */
  void setOn() { on = true; off = false; na = false; } 

  /** Set the flag to OFF
   */
  void setOff() { on = false; off = true; na = false; } 

  /** Query if the flag is NA
   */
  bool isNA()  const { return na; } 

  /** Query if the flag is ON
   */
  bool isOn()  const { return on; } 

  /** Query if the flag is OFF
   */
  bool isOff() const { return off; } 

  /** Set the values of this flag from an input line.
   * 
   * This takes as input a rule string (e.g. !hardclip)
   * and parses it to set the rule. The type of rule it 
   * will look for is specified by the regex (e.g. !?hardclip)
   * @return Returns true if the regex was successfully parsed
   */
  bool parseRuleLine(std::string &val, std::regex &reg);

 private: 

  bool on;
  bool off; 
  bool na;

};

/** Hold a range of valid numeric values (e.g. isize). 
 *
 * Can optionally invert the range to make rule the complement of the range
 * (eg isize NOT in [300,600]
 */
struct Range {

  /** Construct a Range from 0 to 0, non-inverted
   */
  Range() : min(0), max(0), inverted(false), pattern("") {}

  /** Construct a Range from mn to mx, inclusive
   */
  Range(int mn, int mx, int in, std::string p) : min(mn), max(mx), inverted(in), pattern(p) {}

  int min;
  int max;
  bool inverted;
  std::string pattern;
  bool every = true;
  bool none = false;
  
  bool isValid(int val) {
    if (every)
      return true;
    if (none)
      return true;
    if (!inverted)
      return (val >= min && val <= max);
    else
      return (val < min || val > max);
  }

  void parseRuleLine(std::string line);

  friend std::ostream& operator<<(std::ostream &out, const Range &r);

  // set that this ranges accepts everything
  void setEvery() {
    every = true;
    none = false;
  }

  // set that this range accepts nothing
    void setNone() {
    every = false;
    none = true;
   }
  
  // return if this range accepts all values
  bool isEvery() const { return every; }
  
  // return if this range accepts no values
  bool isNone() const { return none; }


};

// a container to hold boolean rules based mostly on alignment flag
struct FlagRule {
  
  FlagRule() {
    dup  = Flag();
    supp       = Flag();
    qcfail     = Flag();
    hardclip   = Flag();
    fwd_strand = Flag();
    rev_strand = Flag();
    mate_fwd_strand = Flag();
    mate_rev_strand = Flag();
    mapped          = Flag();
    mate_mapped     = Flag();
    ff = Flag();
    fr = Flag();
    rf = Flag();
    rr = Flag();
    ic = Flag();
  }
  

  /**
   * if inv is true, then if flag rule is ON and read is ON, return FALSE
   */ 
  /*bool inline flagCheck(Flag &f, bam1_t *b, int bamflag, bool inv) {
    
    if (!f.isNA()) {
      bool val = (b->core.flag & bamflag);
      if ( (f.isOff() && val) || (f.isOn() && !val))
	return inv ? false : true;
    }
    return true; 
    }*/

  Flag dup, supp, qcfail, hardclip, fwd_strand, rev_strand,
    mate_fwd_strand, mate_rev_strand, mapped, mate_mapped, ff, fr, rf, rr, ic;

  bool na = true;
  void parseRuleLine(std::string line);
  
  // ask whether a read passes the rule
  bool isValid(Read &r);

  friend std::ostream& operator<<(std::ostream &out, const FlagRule &fr);

  // set every flag to NA (most permissive)
  void setEvery() {
    dup.setOn();
    supp.setOn();
    qcfail.setOn();
    hardclip.setOn();
    fwd_strand.setOn();
    rev_strand.setOn();
    mate_fwd_strand.setOn();
    mate_rev_strand.setOn();
    mapped.setOn();
    mate_mapped.setOn();
    ff.setOn();
    fr.setOn();
    rf.setOn();
    rr.setOn();
    ic.setOn();
    na = true;
  }

  // set every flag to OFF everythign off)
  void setNone() {
    dup.setOff();
    supp.setOff();
    qcfail.setOff();
    hardclip.setOff();
    fwd_strand.setOff();
    rev_strand.setOff();
    mate_fwd_strand.setOff();
    mate_rev_strand.setOff();
    mapped.setOff();
    mate_mapped.setOff();
    ff.setOff();
    fr.setOff();
    rf.setOff();
    rr.setOff();
    ic.setOff();
  }


  // ask if every flag is set to NA (most permissive)
  bool isEvery() const { return na; }

};

//
class AbstractRule {

 public:

  AbstractRule() {}
  ~AbstractRule() {}

  std::string name = "";
  Range isize = {-1, -1, true, "isize"}; // include all
  Range mapq =  {-1, -1, true, "mapq"}; 
  Range len =   {-1, -1, true, "length"};
  Range clip =  {-1, -1, true, "clip"};
  Range phred = {-1, -1, true, "phred"};
  Range nm = {-1, -1, true, "nm"};
  Range nbases = {-1,-1,true, "nbases"};
  Range ins = {-1,-1,true, "ins"};
  Range del = {-1,-1,true, "del"};
  std::unordered_map<std::string,bool> orientation;

  std::string atm_file = "";
  bool atm_inv = false;
  size_t atm_count = 0;

#ifdef HAVE_AHOCORASICK_AHOCORASICK_H
  //atm_ptr atm;
  AC_AUTOMATA_t * atm = 0;
#endif

  uint32_t subsam_seed = 999;
  double subsam_frac = 1;

  bool none = false;
  // set to true if you want a read to belong to the region if its mate does
  //bool mate = false; 

  FlagRule fr;

  bool isValid(Read &r);

  void parseRuleLine(std::string line);

  void parseSubLine(std::string line);

#ifdef HAVE_AHOCORASICK_AHOCORASICK_H
  bool ahomatch(Read &r);

  bool ahomatch(const char * seq, unsigned len);
#endif

  void parseSeqLine(std::string line);

  friend std::ostream& operator<<(std::ostream &out, const AbstractRule &fr);

  void setEvery() {
    isize.setEvery();
    mapq.setEvery();
    len.setEvery();
    clip.setEvery();
    phred.setEvery();
    nm.setEvery();
    nbases.setEvery();
    fr.setEvery();
    ins.setEvery();
    del.setEvery();
    atm_file = "";
    subsam_frac = 1;
  }
  
  void setNone() { 
    isize.setNone();
    mapq.setNone();
    len.setNone();
    clip.setNone();
    phred.setNone();
    nm.setNone();
    nbases.setNone();
    fr.setNone();
    del.setNone();
    ins.setNone();
    none = true;
  }

  // return if this rule accepts all reads
  bool isEvery() const {
    return ins.isEvery() && del.isEvery() && isize.isEvery() && mapq.isEvery() && len.isEvery() && clip.isEvery() && phred.isEvery() && nm.isEvery() && nbases.isEvery() && fr.isEvery() && (atm_file.length() == 0) && (subsam_frac >= 1);
  }

  // return if this rule accepts no reads
  bool isNone() const {
    return none;
    //return isize.isNone() && mapq.isNone() && len.isNone() && clip.isNone() && phred.isNone() && nm.isNone() && fr.isNone();
  }


};

class MiniRulesCollection;

/** Define a set of rules for creating a variant bam. The syntax is:
   all@!isize:[0,800],mapq:[0,60]
   region@REGION_FILE
   rule1@isize:[0,800],mapq:[0,60]
   rule2@!isize[0,800]:mapq[0,60],:ardclip:supplementary:duplicate:qcfail
   rule3@
   
   A file of NA indicates that the rule should be applied genome-wide.
   The ordering of the lines sets the hierarchical rule. For instance, a rule on line 2 will be applied 
   before a rule on line 3 for all regions that are the union of regions in level 3 and below.
   
   e.g. Level 3 region file has region chr1   100   1000
        Level 2 region file has region chr1   150   1200
	The union of these will produce a new region chr1   100   1200, with level 2
*/
class MiniRules {
  
  friend class MiniRulesCollection;

  public:
  MiniRules() {}
  ~MiniRules() {}
    
  bool isValid(Read &r);
   
  void setRegionFromFile(const std::string& file);

  bool isReadOverlappingRegion(Read &r);

  friend std::ostream& operator<<(std::ostream& out, const MiniRules &mr);
 
  size_t size() const {
    return m_abstract_rules.size();
  }

  void parseDiscordantShortcut(const std::string& line, const AbstractRule& ar);

  bool m_whole_genome = false;

  int m_width = 0;  

  std::string m_region_file;
  //private:

  GRC m_grv;

  int m_level = -1;

  int pad = 0; // how much should we pad the region?

  std::vector<AbstractRule> m_abstract_rules;

  // rule applies to mate too
  bool m_applies_to_mate = false;

  // count the total number of valid reads
  int m_count = 0;

  // pointer to its containing MiniRulesCollection
  MiniRulesCollection * mrc; 
};

// a hierarchy of mini rules to operate on
class MiniRulesCollection {

 public: 
  
  /** Construct an empty MiniRulesCollection 
   * that will pass all reads
   */
  MiniRulesCollection() {}

  MiniRulesCollection(std::string file, bam_hdr_t *b);

  std::string isValid(Read &r);
  
  friend std::ostream& operator<<(std::ostream& out, const MiniRulesCollection &mr);
  
  void sendToBed(std::string file);

  // check if we should do the whole genome
  bool hasWholeGenome() const {
    for (auto it : m_regions)
      if (it.m_whole_genome)
	return true;
    return false;
  }

  std::vector<int> rule_counts;

  GRC getAllRegions() const;

  size_t size() const { return m_regions.size(); } 

  size_t numRules() const {
    size_t num = 0;
    for (auto& it : m_regions)
      num += it.size();
    return num;
  }

  std::vector<MiniRules> m_regions;

  bam_hdr_t * h;// in case we need to convert from text chr to id chr

 private:


  
  
};

}

#endif
