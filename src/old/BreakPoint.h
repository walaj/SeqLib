#ifndef SNOWTOOLS_BREAKPOINT_H__
#define SNOWTOOLS_BREAKPOINT_H__

#include <cstdlib>
#include <string>
#include <unordered_map>
#include <map>
#include "htslib/faidx.h"

#include "SnowTools/GenomicRegion.h"
#include "SnowTools/GenomicRegionCollection.h"
#include "SnowTools/BamRead.h"

#include <unordered_set>
#include "SnowTools/STCoverage.h"

namespace SnowTools {

  // forward declares
  struct BreakPoint;
  
  typedef std::vector<BreakPoint> BPVec;
 
  /** Minimalist BreakEnd used to reduce memory
   *
   * ReducedBreakEnd contains a subset of BreakEnd. 
   * Its principle use is for reading in large files
   * of BreakPoint lines, and de-duping and sending to VCF
   */
 struct ReducedBreakEnd {
   
   ReducedBreakEnd() {}

   ReducedBreakEnd(const GenomicRegion& g, int mq, const std::string & chr_n);
   
   std::string chr_name;
   //char * chr_name;
   GenomicRegion gr;
   uint32_t mapq:8, sub_n:8, nm:16;

 };

 /** Represents a single break-end of a structural variant
  *
  * Created from a contig alignment or directly, stores the 
  * location of a single break-end, and the alignment details
  * of the contig used to create it. A BreakPoint typically consistents
  * of a pair of BreakEnd objects.
  */
 struct BreakEnd {
   
   BreakEnd() {}

   BreakEnd(const GenomicRegion& g, int mq, const std::string & chr_n);
   
   BreakEnd(const BamRead& b);
   
   void checkLocal(const GenomicRegion& window);

   std::string hash() const;

   std::string id;
   std::string chr_name;
   GenomicRegion gr;

   int mapq = -1;
   int cpos = -1;
   int nm = -1;
   int matchlen = -1;

 };

 /** Minimalist BreakPoint used for reading large breakpoint files
  */
 struct ReducedBreakPoint {

   // some helper functions
   char* __string_alloc2char(const std::string& str, char * p) {
     if (!str.empty() && str != "x") {
       p = (char*)malloc(str.length() + 1);
       strcpy(p, str.c_str());
       return p;
     } else {
       return nullptr;
     }
   }
   
   void __smart_check_free(char * p) {
     if (p)
       free(p);
   }

   int getSpan() const {
     if (indel && !insertion) // deletion
       return (abs((int)b1.gr.pos1 - (int)b2.gr.pos1) - 1);
     if (indel) // insertion
       return (strlen(insertion)); // insertion
     if (b1.gr.chr == b2.gr.chr)
       return abs((int)b1.gr.pos1-(int)b2.gr.pos1);
     else
       return -1;

   }

   ReducedBreakPoint() {}
   ~ReducedBreakPoint() {
     __smart_check_free(ref);
     __smart_check_free(alt);
     __smart_check_free(cname);
     __smart_check_free(homology);
     __smart_check_free(insertion);
     __smart_check_free(evidence);
     __smart_check_free(confidence);
     __smart_check_free(repeat);
     //__smart_check_free(read_names);
   }
   ReducedBreakPoint(const std::string &line, bam_hdr_t* h);

   char * ref;
   char * alt;
   char * cname;
   char * evidence;
   char * confidence;
   char * insertion;
   char * homology;
   char * repeat;
   //char * read_names;
   std::string read_names;

   std::vector<std::string> format_s;

   //std::string ref;
   //std::string alt;
   //std::string cname;
   //std::string evidence;
   //std::string confidence;
   //std::string insertion;
   //std::string homology;

   ReducedBreakEnd b1, b2;
   double somatic_score = 0;
   double somatic_lod = 0; // LogOdds that variant not in normal
   double true_lod = 0;

   uint32_t nsplit:8, tsplit:8, af_n:7, num_align:5, secondary:1, dbsnp:1, pass:1, blacklist:1, indel:1, imprecise:1;
   uint32_t tcov_support:8, ncov_support:8, tcov:8, ncov:8;
   uint32_t tcigar:8, ncigar:8, quality:8, af_t:8; 
   uint8_t pon;

 };

 
 /** A structural variation (connnecting A->B) on the genome
  *
  * BreakPoint stores information about an alignment of a contig
  * that indicates a structural variation. It is composed of 
  * two BreakEnd objects and support for the variant
  */
 struct BreakPoint {

   std::string seq, cname, insertion, homology; 

   // the evidence per break-end
   BreakEnd b1, b2;
   
   // reads spanning this breakpoint
   BamReadVector reads;

   bool secondary = false;

   int num_align = 0;
   
   bool isindel = false;

   BreakPoint() {}
   
  /*! @function get the span of the breakpoints (in bp). -1 for interchrom
   * @return int distance between breakpoints
   */
   int getSpan() const;

   std::string getHashString() const;
   
   void order();
   
   bool isEmpty() const { return (b1.gr.pos1 == 0 && b2.gr.pos1 == 0); }
   
   bool operator==(const BreakPoint& bp) const;
   
   bool operator < (const BreakPoint& bp) const { 

     if (b1.gr < bp.b1.gr)
       return true;
     else if (bp.b1.gr < b1.gr)
       return false;
     
     if (b2.gr < bp.b2.gr)
       return true;
     else if (bp.b2.gr < b2.gr)
       return false;
     
     if (cname > bp.cname)
       return true;
     else if (cname < bp.cname)
       return false;
     
     return false;
  }

   friend std::ostream& operator<<(std::ostream& out, const BreakPoint& bp);
   
   bool valid() const;

};

}

#endif
