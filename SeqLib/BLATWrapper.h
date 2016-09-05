#ifndef SEQLIB_BLAT_WRAPPER_H__
#define SEQLIB_BLAT_WRAPPER_H__

/****************************************************************
 ******************* LICENSE AND COPYRIGHT **********************
A large portion of this file and BLATWrapper.cpp are copied
from the BLAT source code, which is copyright of Jim Kent. See 
the BLAT license below:

CONTENTS AND COPYRIGHT

This archive contains the entire source tree for BLAT and
associated utilities.  All files are copyrighted, but license 
is hereby granted for personal, academic, and non-profit use.  
A license is also granted for the contents of the top level 
lib and inc directories for commercial users.  Commercial 
users should contact jim_kent@pacbell.net for access to other modules.
**/

#include <string>
#include <iostream>
#include <algorithm>
#include "SeqLib/BamRecord.h"
#include "SeqLib/BamHeader.h"

extern "C" {

#include "common.h"
#include "memalloc.h"
#include "linefile.h"
#include "bits.h"
#include "hash.h"
#include "dnautil.h"
#include "dnaseq.h"
#include "fa.h"
#include "nib.h"
#include "twoBit.h"
#include "psl.h"
#include "sig.h"
#include "options.h"
#include "obscure.h"
#include "genoFind.h"
#include "trans3.h"
#include "gfClientLib.h"

#include "net.h"
#include "sqlNum.h"
#include "supStitch.h"
#include "genoFind.h"
#include "gfInternal.h"
#include "errabort.h"

}

namespace SeqLib {

  /** Call Blat on a single sequence
   * 
   * BLATWrapper is an interface to the BLAT program (Jim Kent, UCSC, 2002). 
   * The purpose of this class is to allow rapid querying of a sequence
   * against a reference entirely in memory, with no file intermediates
   * or the need to call BLAT externally. 
   */
class BLATWrapper {

 public:

  BLATWrapper() {}

  /** Provide BLAT with information about the reference genome 
   * @param h Header object containing referene genome information
   * @note This is required before running QuerySequence();
   */
  void SetHeaderInfo(const BamHeader& h);

  /** Load the BLAT indexed reference genome and over-represented k-mer file 
   * @param file Path to reference genome
   * @param oocfile Path to the BLAT ooc file (over-represented k-kmer file)
   */
  void LoadIndex(const std::string& file, const std::string& oocfile);
  
  /** Read and query a FASTA file, and write output to disk
   * @param file Path to a FASTA file to BLAT
   * @param ofile Output file to write. Writes in the BLAT "psl" format
   */
  void QueryFile(const std::string& file, const std::string& ofile);

  /** Query a single sequence and return the BLAT alignments 
   * @param name Name of the sequence
   * @param sequence Sequence to query (composed of ACTGN)
   * @param brv Vector of aligned reads that will hold all of the output alignments
   * @note Mapping quality is automatically set to 0, but the alignment score provided
   * by BLAT is provided in the AS tag.
   */
  void AlignSequence(const std::string& name, const std::string& sequence, BamRecordVector& brv);
  
  //template<class Archive>
  // void serialize(Archive & ar, const unsigned int version) {
  //  for (int i = 0; i < slCount(dbSeqList); ++i)
  //    ar & gf->lists[i];
  //}

 private:

  std::unordered_map<std::string, int> m_name2id;

  char** dbFiles;
  struct dnaSeq *dbSeqList;
  struct genoFind *gf;

  struct gfOutput *gvo;		/* Overall output controller */

  int tileSize = 11;
  int stepSize = 0;	/* Default (same as tileSize) */
  int minMatch = 2;
  int minScore = 30;
  int maxGap = 2;
  int repMatch = 1024; //*4;
  int dotEvery = 0;
  boolean oneOff = false;
  boolean noHead = false;
  boolean trimA = false;
  boolean trimHardA = false;
  boolean trimT = false;
  boolean fastMap = false;
  char *makeOoc = NULL;
  char *ooc = NULL;
  enum gfType qType = gftDna;
  enum gfType tType = gftDna;
  char *mask = NULL;
  char *repeats = NULL;
  char *qMask = NULL;
  double minRepDivergence = 15;
  double minIdentity = 90;
  char *outputFormat = "psl";
  
  void searchOneStrand(struct dnaSeq *seq, struct genoFind *gf, 
		       boolean isRc, Bits *qMaskBits, BamRecordVector& brv);

  Bits* maskQuerySeq(struct dnaSeq *seq, boolean isProt, 
				  boolean maskQuery, boolean lcMask);

  void searchOne(bioSeq *seq, struct genoFind *gf, struct hash *maskHash, Bits *qMaskBits, BamRecordVector& brv);

  void trimSeq(struct dnaSeq *seq, struct dnaSeq *trimmed);

  void __searchOneIndex(int fileCount, char *files[], struct genoFind *gf, char *outName, 
		        struct hash *maskHash, FILE *outFile, boolean showStatus);
    

  void searchOneMaskTrim(struct dnaSeq *seq, boolean isProt,
		       struct genoFind *gf,
		       struct hash *maskHash,
			 long long *retTotalSize, int *retCount, BamRecordVector& brv);

  void __gfLongDnaInMem(struct dnaSeq *query, struct genoFind *gf, 
			boolean isRc, int minScore, Bits *qMaskBits, 
			boolean fastMap, boolean band, BamRecordVector& brv);

  struct ssBundle* __fastMapClumpsToBundles(struct genoFind *gf, struct gfClump *clumpList, bioSeq *qSeq);

  void __clumpToHspRange(struct gfClump *clump, bioSeq *qSeq, int tileSize,
			 int frame, struct trans3 *t3, struct gfRange **pRangeList, 
			 boolean isProt, boolean fastMap);
  
  int __scoreAli(struct ffAli *ali, boolean isProt, 
		 enum ffStringency stringency, 
		 struct dnaSeq *tSeq, struct trans3 *t3List);
  
  struct ssBundle* __gfClumpsToBundles(struct gfClump *clumpList, 
				       boolean isRc, struct dnaSeq *seq, int minScore,  
				       struct gfRange **retRangeList);
  
  void __addToBigBundleList(struct ssBundle **pOneList, struct hash *bunHash, 
			    struct ssBundle **pBigList, struct dnaSeq *query);
  

  struct gfRange* __seqClumpToRangeList(struct gfClump *clumpList, int frame);

  boolean __alignComponents(struct gfRange *combined, struct ssBundle *bun, enum ffStringency stringency);

  char* __clumpTargetName(struct gfClump *clump);

  void __extendHitLeft(int qMax, int tMax,
				  char **pStartQ, char **pStartT, int (*scoreMatch)(char a, char b),
		       int maxDown);
  

  void __extendHitRight(int qMax, int tMax,
			char **pEndQ, char **pEndT, int (*scoreMatch)(char a, char b), 
			int maxDown);

  struct gfHit* __gfFindHitsWithQmask(struct genoFind *gf, bioSeq *seq,
				      Bits *qMaskBits, int qMaskOffset, struct lm *lm, int *retHitCount, 
				      struct gfSeqSource *target, int tMin, int tMax);
  

  struct gfHit* __gfFastFindDnaHits(struct genoFind *gf, struct dnaSeq *seq, 
				    Bits *qMaskBits,  int qMaskOffset, struct lm *lm, int *retHitCount,
				    struct gfSeqSource *target, int tMin, int tMax);
    

  struct gfSeqSource* __findSource(struct genoFind *gf, bits32 targetPos);
  
};
 
}

#endif
