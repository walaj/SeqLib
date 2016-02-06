#ifndef SNOWTOOLS_BLAT_WRAPPER_H__
#define SNOWTOOLS_BLAT_WRAPPER_H__

#include <string>
#include <iostream>
#include <algorithm>
#include "SnowTools/BamRead.h"
#include "SnowTools/BamWalker.h"

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

namespace SnowTools {

class BLATWrapper {

 public:

  BLATWrapper() {}

  void loadIndex(const std::string& file, const std::string& oocfile);
  
  void queryFile(const std::string& file);

  void querySequence(const std::string& name, const std::string& sequence, BamReadVector& brv);
  
 private:

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
  
  void searchOneStrand(struct dnaSeq *seq, struct genoFind *gf, FILE *psl, 
		       boolean isRc, struct hash *maskHash, Bits *qMaskBits, BamReadVector& brv);

  Bits* maskQuerySeq(struct dnaSeq *seq, boolean isProt, 
				  boolean maskQuery, boolean lcMask);

  void searchOne(bioSeq *seq, struct genoFind *gf, FILE *f, boolean isProt, struct hash *maskHash, Bits *qMaskBits, BamReadVector& brv);

  void trimSeq(struct dnaSeq *seq, struct dnaSeq *trimmed);

  void __searchOneIndex(int fileCount, char *files[], struct genoFind *gf, char *outName, 
			boolean isProt, struct hash *maskHash, FILE *outFile, boolean showStatus);
    

  void searchOneMaskTrim(struct dnaSeq *seq, boolean isProt,
		       struct genoFind *gf, FILE *outFile,
		       struct hash *maskHash,
			 long long *retTotalSize, int *retCount, BamReadVector& brv);

  void __gfLongDnaInMem(struct dnaSeq *query, struct genoFind *gf, 
			boolean isRc, int minScore, Bits *qMaskBits, 
			struct gfOutput *out, boolean fastMap, boolean band, BamReadVector& brv);

  struct ssBundle* __fastMapClumpsToBundles(struct genoFind *gf, struct gfClump *clumpList, bioSeq *qSeq);

  struct ffAli* __refineSmallExons(struct ffAli *ff, 
				   struct dnaSeq *nSeq, struct dnaSeq *hSeq);
  
  void __clumpToHspRange(struct gfClump *clump, bioSeq *qSeq, int tileSize,
			 int frame, struct trans3 *t3, struct gfRange **pRangeList, 
			 boolean isProt, boolean fastMap);
  
  void __refineSmallExonsInBundle(struct ssBundle *bun);
  
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
  
  
};
 
}

#endif
