#ifndef SNOWTOOLS_HTSTOOLS_H__
#define SNOWTOOLS_HTSTOOLS_H__

#include <memory>
#include <vector>
#include <unordered_map>
#include <string>
#include <algorithm>

#include "htslib/hts.h"
#include "htslib/sam.h"
#include "htslib/bgzf.h"
#include "htslib/kstring.h"


static const char BASES[16] = {' ', 'A', 'C', ' ',
                               'G', ' ', ' ', ' ', 
                               'T', ' ', ' ', ' ', 
                               ' ', ' ', ' ', 'N'};

struct free_delete {
  void operator()(void* x) { bam_destroy1((bam1_t*)x); }
};

typedef std::shared_ptr<bam1_t> Read;
typedef std::vector<Read> ReadVec;

#define r_is_rev(b) (((b)->core.flag&BAM_FREVERSE) != 0)
#define r_is_mrev(b) (((b)->core.flag&BAM_FMREVERSE) != 0)
#define r_pos(b) ((b)->core.pos)
#define r_mpos(b) ((b)->core.mpos)
#define r_id(b) ((b)->core.tid)
#define r_mid(b) ((b)->core.mtid)
#define r_is_mapped(b) (!((b)->core.flag&BAM_FUNMAP))
#define r_is_mmapped(b) (!((b)->core.flag&BAM_FMUNMAP))
#define r_mapq(b) ((b)->core.qual)
#define r_cig_size(b) ((b)->core.n_cigar)
#define r_is_first(b) ((b)->core.flag&BAM_FREAD1)
#define r_qname(b) (std::string(bam_get_qname((b).get())))
#define r_remove_tag(b, t) do { uint8_t *p  = bam_aux_get((b).get(), t); if (p) bam_aux_del((b).get(), p); } while(0)
#define r_flag(b) ((b)->core.flag)
#define r_is_pmapped(b) (!((b)->core.flag&BAM_FUNMAP) && !((b)->core.flag&BAM_FMUNMAP) && ((b)->core.flag&BAM_FPAIRED))
#define r_length(b) ((b)->core.l_qseq)
#define r_isize(b) ((b)->core.isize)
#define r_seq(b, s) do { uint8_t * p = bam_get_seq(b); char c[(b)->core.l_qseq+1]; for (int ww = 0; ww < (b)->core.l_qseq; ww++) { c[ww] = BASES[bam_seqi(p, ww)]; } c[(b)->core.l_qseq] = '\0'; s = std::string(c); } while(0)
#define r_get_clip(b, p) do { uint32_t* c = bam_get_cigar(b); for (int ww = 0; ww < (b)->core.n_cigar; ww++) if (c[ww] & BAM_CSOFT_CLIP || c[ww] & BAM_CHARD_CLIP) p += bam_cigar_oplen(c[ww]); } while(0)

#define r_get_Z_tag(b, t, v) do { uint8_t* p = bam_aux_get((b).get(), t); if (p) v = std::string(bam_aux2Z(p)); else v = std::string(); } while(0) 
#define r_get_SR(b, v) do { uint8_t* p = bam_aux_get((b).get(), "SR"); assert(p); v = bam_aux2Z(p); } while(0) 
#define r_get_int32_tag(b, t, v) do { uint8_t* p = bam_aux_get((b).get(), t); if (p) v = bam_aux2i(p); else v = 0; } while(0) 

#define r_cig_type(b, i) (bam_cigar_opchr(bam_get_cigar((b).get())[i]))
#define r_cig_len(b, i) (bam_cigar_oplen(bam_get_cigar((b).get())[i]))

#define r_add_Z_tag(b, t, v) bam_aux_append((b).get(), t, 'Z', (v).length()+1, (uint8_t*)(v).c_str())
#define r_add_int32_tag(b, t, v) bam_aux_append((b).get(), t, 'i', 4, (uint8_t*)&(v))
#define r_count_nbases(b, n) do { uint8_t* p = bam_get_seq(b); for (int ww = 0; ww < (b)->core.l_qseq; ww++) if (bam_seqi(p,ww) == 15) n++; } while(0)
#define r_count_sub_nbases(b, n, s, e) do { uint8_t* p = bam_get_seq(b); for (int ww = (s); ww < (e); ww++) if (bam_seqi(p,ww) == 15) n++; } while(0)

#define r_is_rstrand(b) (((b)->core.flag&BAM_FREVERSE) > 0) 
#define r_is_mrstrand(b) (((b)->core.flag&BAM_FMREVERSE) > 0)
#define r_is_dup(b) (((b)->core.flag&BAM_FDUP) > 0)
#define r_is_primary(b) (!((b)->core.flag&BAM_FSECONDARY))
#define r_is_qc_fail(b) (((b)->core.flag&BAM_FQCFAIL) > 0)
#define r_strand(b) ((b)->core.flag&BAM_FREVERSE ? '-' : '+')
#define r_endpos(b) (bam_endpos((b).get()))

#define r_get_trimmed_seq(b, s) do { int32_t ts, tl; r_get_int32_tag(b, "TS", ts);  r_get_int32_tag(b, "TL", tl); r_seq(b, s); if (tl != 0 && tl != r_length(b))  s = s.substr(ts, tl); } while (0) 

  namespace SnowTools {

    int32_t qualityTrimRead(int qualTrim, int32_t &startpoint, Read &r);

    std::vector<std::string> GetStringTag(const Read& a, const std::string tag);
    
    void SmartAddTag(Read &a, const std::string tag, const std::string val);
    
    std::vector<int> GetIntTag(const Read& a, const std::string tag);

    void rcomplement(std::string& a);

    void removeAllTags(Read& a);
  }

#endif
