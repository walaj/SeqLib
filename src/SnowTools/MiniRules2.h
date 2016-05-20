#ifndef MINI_RULES_H
#define MINI_RULES_H

#include <string>
#include <vector>
#include <unordered_map>

#include "json/json.h"

#include "SnowTools/GenomicRegionCollection.h"
#include "SnowTools/BamRead.h"

//#define HAVE_AHOCORASICK_AHOCORASICK_H 1
//#ifdef HAVE_AHOCORASICK_AHOCORASICK_H
#ifndef __APPLE__

#include "ahocorasick/ahocorasick.h"
#include <memory>

// custom deleter for aho-corasick
struct atm_free_delete {
  void operator()(void* x) { free((AC_AUTOMATA_t*)x); }
};
typedef std::unique_ptr<AC_AUTOMATA_t> atm_ptr;
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

  bool parseJson(const Json::Value& value, const std::string& name);

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

  void parseJson(const Json::Value& value, const std::string& name);

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

/** Stores a set of Flag objects for filtering alignment flags
 *
 * An alignment can be queried against a FlagRule to check if it 
 * satisfies the requirements for its alignment flag.
 */
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
    paired = Flag();
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
    mate_fwd_strand, mate_rev_strand, mapped, mate_mapped, ff, fr, rf, rr, ic, paired;

  bool na = true;

  void parseJson(const Json::Value& value);
  
  // ask whether a read passes the rule
  bool isValid(BamRead &r);

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
    paired.setOn();
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
    paired.setOff();
  }


  // ask if every flag is set to NA (most permissive)
  bool isEvery() const { return na; }

};

/** Stores a full rule (Flag + Range + motif etc)
 *
 * An alignment can be queried with an AbstractRule object
 * to check if it passes that rule.
 */
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
  Range xp = {-1,-1,true, "xp"};
  std::unordered_map<std::string,bool> orientation;

  std::string atm_file = "";
  bool atm_inv = false;
  size_t atm_count = 0;

  std::string id;

  // how many reads pass this rule?
  size_t m_count = 0;

  //#ifdef HAVE_AHOCORASICK_AHOCORASICK_H
#ifndef __APPLE__
  //atm_ptr atm;
  AC_AUTOMATA_t * atm = 0;
#endif

  uint32_t subsam_seed = 999;
  double subsam_frac = 1;

  bool none = false;
  // set to true if you want a read to belong to the region if its mate does
  //bool mate = false; 

  FlagRule fr;

  bool isValid(BamRead &r);



  //#ifdef HAVE_AHOCORASICK_AHOCORASICK_H
#ifndef __APPLE__
  bool ahomatch(BamRead &r);

  bool ahomatch(const char * seq, unsigned len);
#endif

  void parseJson(const Json::Value& value);
  void parseSubLine(const Json::Value& value);
  void parseSeqLine(const Json::Value& value);

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
    xp.setEvery();
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
    xp.setNone();
    ins.setNone();
    none = true;
  }

  // return if this rule accepts all reads
  bool isEvery() const {
    return ins.isEvery() && del.isEvery() && isize.isEvery() && mapq.isEvery() && len.isEvery() && clip.isEvery() && phred.isEvery() && nm.isEvery() && nbases.isEvery() && fr.isEvery() && (atm_file.length() == 0) && (subsam_frac >= 1) && xp.isEvery();
  }

  // return if this rule accepts no reads
  bool isNone() const {
    return none;
    //return isize.isNone() && mapq.isNone() && len.isNone() && clip.isNone() && phred.isNone() && nm.isNone() && fr.isNone();
  }


};

class MiniRulesCollection;

/** 
 * A set of AbstractRules on a region united by logi rules.
 *
 * MiniRules stores an arbitrarily complex collection of AbstractRules
 * (e.g. (Mapped && Clipped) || (Unmapped)).
 */
class MiniRules {
  
  friend class MiniRulesCollection;

  public:
  MiniRules() {}
  ~MiniRules() {}
  
  std::string id;
  
  bool excluder = false; // this region is for excluding

  bool isValid(BamRead &r);
   
  void setRegionFromFile(const std::string& file);

  bool isReadOverlappingRegion(BamRead &r);

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

  // pointer to its containing MiniRulesCollection
  MiniRulesCollection * mrc; 

  // how many reads pass this MiniRule
  size_t m_count = 0;

};

/** A full set of rules across any number of regions
 *
 * Stores the entire set of MiniRules, each defined on a unique interval.
 * A single MiniRulesCollection object is sufficient to store any combination of rules,
 * and is the highest in the rule hierarchy. (MiniRulesCollection stores MiniRules 
 * stores AbstractRules stores FlagRule/Ranges).
 */
class MiniRulesCollection {

 public: 
  
  /** Construct an empty MiniRulesCollection 
   * that will pass all reads
   */
  MiniRulesCollection() {}

  AbstractRule rule_all;
  
  MiniRulesCollection(const std::string& file, bam_hdr_t *b);

  MiniRulesCollection(const std::string& file);

  void addGlobalRule(const std::string& rule);

  bool isValid(BamRead &r);
  
  friend std::ostream& operator<<(std::ostream& out, const MiniRulesCollection &mr);
  
  void sendToBed(std::string file);

  size_t m_count = 0; // passed
  size_t m_count_seen = 0; // tested

  GRC getAllRegions() const;

  size_t size() const { return m_regions.size(); } 

  size_t numRules() const {
    size_t num = 0;
    for (auto& it : m_regions)
      num += it.size();
    return num;
  }

  std::vector<MiniRules> m_regions;

  bam_hdr_t * h = nullptr;// in case we need to convert from text chr to id chr

  void countsToFile(const std::string& file) const;

  // should we keep checking rules, even it passed? (useful for counting)
  bool m_fall_through = false;

 private:  

  void __construct_MRC(const std::string& file);  

  const std::string GetScriptContents(const std::string& script);

  bool ParseFilterObject(const std::string& filterName, const Json::Value& filterObject);
};

}

#endif
