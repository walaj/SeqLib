#ifndef MINI_RULES_H
#define MINI_RULES_H

#include <string>
#include <vector>
#include <unordered_map>
#include <climits>

#include "json/json.h"

#include "SnowTools/GenomicRegionCollection.h"
#include "SnowTools/BamRead.h"

// motif matching with ahocorasick not available on OSX
#ifndef __APPLE__
#include "ahocorasick/ahocorasick.h"
#include <memory>

// custom deleter for aho-corasick
struct atm_free_delete {
  void operator()(void* x) { free((AC_AUTOMATA_t*)x); }
};
typedef std::unique_ptr<AC_AUTOMATA_t> atm_ptr;
#endif

#define MINIRULES_MATE_LINKED 1
#define MINIRULES_MATE_LINKED_EXCLUDE 2
#define MINIRULES_REGION 3
#define MINIRULES_REGION_EXCLUDE 4

namespace SnowTools {

struct CommandLineRegion {
  
CommandLineRegion(const std::string& mf, int t) : f(mf), type(t), pad(0), i_flag(0), e_flag(0) {}

  std::string f; // file
  int type; // mate linked, excluder, etc
  int pad;
  uint32_t i_flag; // inclusive flags
  uint32_t e_flag; // exclusive flags

  int len = 0;
  int mapq = 0;
  int nbases = INT_MAX;
  int phred = 0;
  int clip = 0;
  int ins = 0;
  int del = 0;
  std::string rg, motif;

  bool all() const { 
    return !len && !mapq && !nbases && !phred && rg.empty() && !i_flag && !e_flag; 
  }

};


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

  /** Construct a default range with everything accepted
   */
  Range() : m_min(0), m_max(0), m_inverted(false), m_every(true) {}

  /** Construct a Range from min to max, inclusive
   * @param min Minimum for range 
   * @param max Maximum for range 
   * @param inverted Declare if this should be an inverted range (do NOT accept vals in range)
   */
  Range(int min, int max, bool inverted) : m_min(min), m_max(max), m_inverted(inverted), m_every(false) {}
  
  /** Given a query value, determine if the value passes this Range
   * @param val Query value (e.g. mapping quality)
   * @return true if the value passes this Range rule
   */
  bool isValid(int val) {
    if (m_every)
      return true;
    if (!m_inverted)
      return (val >= m_min && val <= m_max);
    else
      return (val < m_min || val > m_max);
  }

  /** Parse a JSON value 
   * @param value 
   * @param name
   */
  void parseJson(const Json::Value& value, const std::string& name);

  /** Print the contents of this Range */
  friend std::ostream& operator<<(std::ostream &out, const Range &r);

  /** Return if this range accepts all values */
  bool isEvery() const { return m_every; }
  
  /** Return the lower bound of the range */
  int lowerBound() const { return m_min; }

  /** Return the upper bound of the range */
  int upperBound() const { return m_max; }

  /** Return true if the range is inverted (e.g. do NOT accept i in [min,max] */
  bool isInverted() const { return m_inverted; }
  
private:
  
  int m_min;
  int m_max;
  bool m_inverted;
  bool m_every;

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
    m_on_flag = 0;
    m_off_flag = 0;
  }
  
  Flag dup, supp, qcfail, hardclip, fwd_strand, rev_strand,
    mate_fwd_strand, mate_rev_strand, mapped, mate_mapped, ff, fr, rf, rr, ic, paired;

  bool na = true;

  void parseJson(const Json::Value& value);

  void setOnFlag(uint32_t f) { m_on_flag = f; na = na && f == 0; } 

  void setOffFlag(uint32_t f) { m_off_flag = f; na = na && f == 0; } 

  // ask whether a read passes the rule
  bool isValid(BamRead &r);

  /** Print the flag rule */
  friend std::ostream& operator<<(std::ostream &out, const FlagRule &fr);

  // ask if every flag is set to NA (most permissive)
  bool isEvery() const { return na; }

private:

  uint32_t m_on_flag;
  uint32_t m_off_flag;

  int __parse_json_int(const Json::Value& v);

};

/** Stores a full rule (Flag + Range + motif etc)
 *
 * An alignment can be queried with an AbstractRule object
 * to check if it passes that rule.
 */
class AbstractRule {

 public:

  /** Create empty rule with default to accept all */
  AbstractRule() {}

  /** Destroy */
  ~AbstractRule() {}

  Range isize, mapq, len, clip, phred, nm, nbases, ins, del, xp;

  void addMotifRule(const std::string& f, bool inverted);


  std::string id;

  // read group 
  std::string read_group;

  // how many reads pass this rule?
  size_t m_count = 0;

#ifndef __APPLE__
  //atm_ptr atm;
  AC_AUTOMATA_t * atm = 0;
#endif


  double subsam_frac = 1;

  bool isValid(BamRead &r);

  void parseJson(const Json::Value& value);

  friend std::ostream& operator<<(std::ostream &out, const AbstractRule &fr);

  // return if this rule accepts all reads
  bool isEvery() const {
    return read_group.empty() && ins.isEvery() && del.isEvery() && isize.isEvery() && mapq.isEvery() && len.isEvery() && clip.isEvery() && phred.isEvery() && nm.isEvery() && nbases.isEvery() && fr.isEvery() && (atm_file.length() == 0) && (subsam_frac >= 1) && xp.isEvery();
  }

  FlagRule fr;

 private:

  // data
  uint32_t subsam_seed = 999; // random seed for subsampling

  // motif data
  std::string atm_file; // sequence file
  bool atm_inv = false; // is this inverted
  size_t atm_count = 0; // number of motifs

#ifndef __APPLE__
  bool ahomatch(BamRead &r);
  bool ahomatch(const char * seq, unsigned len);
#endif

  void parseSubLine(const Json::Value& value);

  void parseSeqLine(const Json::Value& value);

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

  /** Make a MiniRules with a an all exclude or include rule
   * @param Samtools style string, BED file or VCF 
   * @param reg_type The type of rule this will be
   * @param h BAM header that defines available chromosomes
   */
  MiniRules(const CommandLineRegion& c, bam_hdr_t * h);

  std::string id;
  
  bool excluder = false; // this region is for excluding

  bool isValid(BamRead &r);
   
  void setRegionFromFile(const std::string& file, bam_hdr_t * h);

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

  GRC m_grv; // the interval tree with the regions this rule applies to

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
  
  MiniRulesCollection(const std::string& script, bam_hdr_t *h);

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

  const std::string GetScriptContents(const std::string& script);

  bool ParseFilterObject(const std::string& filterName, const Json::Value& filterObject);

  bool __validate_json_value(const Json::Value value, const std::unordered_set<std::string>& valid_vals);
};

}

#endif
