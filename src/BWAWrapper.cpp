#include "SnowTools/BWAWrapper.h"

#include <sstream>
#include <iostream>

extern "C" {
  #include <string.h>
}

#define DEBUG_BWATOOLS 1

#define _set_pac(pac, l, c) ((pac)[(l)>>2] |= (c)<<((~(l)&3)<<1))
#define _get_pac(pac, l) ((pac)[(l)>>2]>>((~(l)&3)<<1)&3)

namespace SnowTools {

  void BWAWrapper::constructIndex(const USeqVector& v) {

    if (idx) {
      std::cout << "...clearing old index" << std::endl;
      bwa_idx_destroy(idx);
      idx = 0;
    }
    
    // allocate memory for idx
    idx = (bwaidx_t*)calloc(1, sizeof(bwaidx_t));;

    // construct the forward-only pac
    uint8_t* fwd_pac = __make_pac(v, true, false); //true->for_only, false->write_file

    // construct the forward-reverse pac ("packed" 2 bit sequence)
    uint8_t* pac = __make_pac(v, false, false); // don't write, becasue only used to make BWT

    size_t tlen = 0;
    for (auto& i : v)
      tlen += i.seq.length();

#ifdef DEBUG_BWATOOLS
    //std::cerr << "ref seq length: " << tlen << std::endl;
#endif

    // make the bwt
    bwt_t *bwt;
    bwt = __bwt_pac2bwt(pac, tlen*2); // *2 for fwd and rev
    bwt_bwtupdate_core(bwt);
    free(pac); // done with fwd-rev pac 
    
    // construct sa from bwt and occ. adds it to bwt struct
    bwt_cal_sa(bwt, 32);
    bwt_gen_cnt_table(bwt);
        
    // make the bns
    bntseq_t * bns = (bntseq_t*) calloc(1, sizeof(bntseq_t));
    bns->l_pac = tlen;
    bns->n_seqs = v.size();
    bns->seed = 11;
    bns->n_holes = 0;

    // make the anns
    bns->anns = (bntann1_t*)calloc(v.size(), sizeof(bntann1_t));
    size_t offset = 0;
    for (size_t k = 0; k < v.size(); ++k) {
      __add_to_anns(v[k].name, v[k].seq, &bns->anns[k], offset);
      offset += v[k].seq.length();
    }
    
    //ambs is "holes", like N bases
    bns->ambs = 0; //(bntamb1_t*)calloc(1, sizeof(bntamb1_t));
    
    // make the in-memory idx struct
    idx->bwt = bwt;
    idx->bns = bns;
    idx->pac = fwd_pac;

    return;
    
  }

  void BWAWrapper::alignSingleSequence(const std::string& seq) {

    mem_alnreg_v ar;
    ar = mem_align1(memopt, idx->bwt, idx->bns, idx->pac, seq.length(), seq.c_str()); // get all the hits

#ifdef DEBUG_BWATOOLS
    std::cout << "num hits: " << ar.n << std::endl;
    //std::cout << __print_bns() << std::endl;
#endif    

    // loop through the hits
    for (size_t i = 0; i < ar.n; ++i) {
      mem_aln_t a;
      if (ar.a[i].secondary >= 0) 
	continue; // skip secondary alignments
      a = mem_reg2aln(memopt, idx->bns, idx->pac, seq.length(), seq.c_str(), &ar.a[i]); // get forward-strand position and CIGAR

      BamRead b;
      b.SetQname("testing123");
      b.SetSequence(seq);
      b.SetMapQuality(a.mapq);
      //std::cout << b.toSam(h) << std::endl;
#ifdef DEBUG_BWATOOLS
      // print alignment
      printf("\t%c\t%s\t%ld\t%d\t", /*ks->name.s,*/ "+-"[a.is_rev], idx->bns->anns[a.rid].name, (long)a.pos, a.mapq);
      for (int k = 0; k < a.n_cigar; ++k) // print CIGAR
	printf("%d%c", a.cigar[k]>>4, "MIDSH"[a.cigar[k]&0xf]);
      printf("\t%d\n", a.NM); // print edit distance
#endif
      free(a.cigar); // don't forget to deallocate CIGAR
    }
    free (ar.a); // dealloc the hit list
}

  //void BWAWrapper::alignReads(const std::vector<BamRead>& reads){}

uint8_t* BWAWrapper::__add1(const kseq_t *seq, bntseq_t *bns, uint8_t *pac, int64_t *m_pac, int *m_seqs, int *m_holes, bntamb1_t **q)
{
  bntann1_t *p;
  int lasts;
  if (bns->n_seqs == *m_seqs) {
    *m_seqs <<= 1;
    bns->anns = (bntann1_t*)realloc(bns->anns, *m_seqs * sizeof(bntann1_t));
  }
  p = bns->anns + bns->n_seqs;
  p->name = strdup((char*)seq->name.s);
  p->anno = seq->comment.l > 0? strdup((char*)seq->comment.s) : strdup("(null)");
  p->gi = 0; p->len = seq->seq.l;
  p->offset = (bns->n_seqs == 0)? 0 : (p-1)->offset + (p-1)->len;
  p->n_ambs = 0;
  for (size_t i = lasts = 0; i < seq->seq.l; ++i) {
    int c = nst_nt4_table[(int)seq->seq.s[i]];
    if (c >= 4) { // N
      if (lasts == seq->seq.s[i]) { // contiguous N
	++(*q)->len;
      } else {
	if (bns->n_holes == *m_holes) {
	  (*m_holes) <<= 1;
	  bns->ambs = (bntamb1_t*)realloc(bns->ambs, (*m_holes) * sizeof(bntamb1_t));
	}
	*q = bns->ambs + bns->n_holes;
	(*q)->len = 1;
	(*q)->offset = p->offset + i;
	(*q)->amb = seq->seq.s[i];
	++p->n_ambs;
	++bns->n_holes;
      }
    }
    lasts = seq->seq.s[i];
    { // fill buffer
      if (c >= 4) c = lrand48()&3;
      if (bns->l_pac == *m_pac) { // double the pac size
	*m_pac <<= 1;
	pac = (uint8_t*)realloc(pac, *m_pac/4);
	memset(pac + bns->l_pac/4, 0, (*m_pac - bns->l_pac)/4);
      }
      _set_pac(pac, bns->l_pac, c);
      ++bns->l_pac;
    }
  }
  ++bns->n_seqs;

  return pac;
}

uint8_t* BWAWrapper::__make_pac(const USeqVector& v, bool for_only, bool write_file)
{

  bntseq_t * bns = (bntseq_t*)calloc(1, sizeof(bntseq_t));
  uint8_t *pac = 0; 
  int32_t m_seqs, m_holes;
  int64_t m_pac, l;
  bntamb1_t *q;

  bns->seed = 11; // fixed seed for random generator
  m_seqs = m_holes = 8; m_pac = 0x10000;
  bns->anns = (bntann1_t*)calloc(m_seqs, sizeof(bntann1_t));
  bns->ambs = (bntamb1_t*)calloc(m_holes, sizeof(bntamb1_t));
  pac = (uint8_t*) calloc(m_pac/4, 1);
  q = bns->ambs;

  // move through the unaligned sequences
  for (size_t k = 0; k < v.size(); ++k) {

    // make the ref name kstring
    kstring_t * name = (kstring_t*)malloc(1 * sizeof(kstring_t));
    name->l = v[k].name.length();
    name->m = v[k].name.length() + 2;
    name->s = (char*)calloc(v[k].name.length(), sizeof(char));
    strncpy(name->s, v[k].name.c_str(), v[k].name.length());
    
    // make the sequence kstring
    kstring_t * t = (kstring_t*)malloc(sizeof(kstring_t));
    t->l = v[k].seq.length();
    t->m = v[k].seq.length() + 2;
    t->s = (char*)calloc(v[k].seq.length(), sizeof(char));
    strncpy(t->s, v[k].seq.c_str(), v[k].seq.length());

    // put into a kstring
    kseq_t *ks = (kseq_t*)calloc(1, sizeof(kseq_t));  
    ks->seq = *t;
    ks->name = *name;
    
    // make the forward only pac
    pac = __add1(ks, bns, pac, &m_pac, &m_seqs, &m_holes, &q);

    // clear it out
    free(name->s);
    free(name);
    free(t->s);
    free(t);
    //free(ks->name.s); 
    //free(ks->seq.s);
    //free(ks->f->buf);
    //free(
    free(ks);
    // NOTE free kstring_t?
    //kseq_destroy(s);
  }

  if (!for_only) 
    {
      // add the reverse complemented sequence
      m_pac = (bns->l_pac * 2 + 3) / 4 * 4;
      pac = (uint8_t*)realloc(pac, m_pac/4);
      memset(pac + (bns->l_pac+3)/4, 0, (m_pac - (bns->l_pac+3)/4*4) / 4);
      for (l = bns->l_pac - 1; l >= 0; --l, ++bns->l_pac)
	_set_pac(pac, bns->l_pac, 3-_get_pac(pac, l));
    }

  if (write_file)
  { // finalize .pac file
    FILE *fp;
    fp = xopen("mem_test.pac", "wb");
    ubyte_t ct;
    err_fwrite(pac, 1, (bns->l_pac>>2) + ((bns->l_pac&3) == 0? 0 : 1), fp);
    // the following codes make the pac file size always (l_pac/4+1+1)
    if (bns->l_pac % 4 == 0) {
      ct = 0;
      err_fwrite(&ct, 1, 1, fp);
    }
    ct = bns->l_pac % 4;
    err_fwrite(&ct, 1, 1, fp);
    // close .pac file
    err_fflush(fp);
    err_fclose(fp);
  }

  bns_destroy(bns);
  
  return pac;
}

bwt_t *BWAWrapper::__bwt_pac2bwt(const uint8_t *pac, int bwt_seq_lenr)
{

  bwt_t *bwt;
  ubyte_t *buf;
  int i;
  //FILE *fp;

  // initialization
  bwt = (bwt_t*)calloc(1, sizeof(bwt_t));
  bwt->seq_len = bwt_seq_lenr; //bwa_seq_len(fn_pac); //dummy
  bwt->bwt_size = (bwt->seq_len + 15) >> 4;
  //fp = xopen(fn_pac, "rb");

  // prepare sequence
  //pac_size = (bwt->seq_len>>2) + ((bwt->seq_len&3) == 0? 0 : 1);
  //buf2 = (ubyte_t*)calloc(pac_size, 1);
  //err_fread_noeof(buf2, 1, pac_size, fp);
  //err_fclose(fp);
  memset(bwt->L2, 0, 5 * 4);
  buf = (ubyte_t*)calloc(bwt->seq_len + 1, 1);
  for (i = 0; i < (int)bwt->seq_len; ++i) {
    buf[i] = pac[i>>2] >> ((3 - (i&3)) << 1) & 3;
    ++bwt->L2[1+buf[i]];
  }
  for (i = 2; i <= 4; ++i) 
    bwt->L2[i] += bwt->L2[i-1];
  //free(buf2);

  // Burrows-Wheeler Transform
  bwt->primary = is_bwt(buf, bwt->seq_len);
  bwt->bwt = (u_int32_t*)calloc(bwt->bwt_size, 4);
  for (i = 0; i < (int)bwt->seq_len; ++i)
    bwt->bwt[i>>4] |= buf[i] << ((15 - (i&15)) << 1);
  free(buf);
  return bwt;
}

  bntann1_t* BWAWrapper::__add_to_anns(const std::string& name, const std::string& seq, bntann1_t* ann, size_t offset) 
  {

    ann->offset = offset;
    ann->name = (char*)malloc(name.length()+1); // +1 for \0
    strncpy(ann->name, name.c_str(), name.length()+1);
    ann->anno = (char*)malloc(7);
    strcpy(ann->anno, "(null)\0");
    ann->len = seq.length();
    ann->n_ambs = 0; // number of "holes"
    ann->gi = 0; // gi?
    ann->is_alt = 0;
    
    return ann;
  }

  
  void BWAWrapper::retrieveIndex(const std::string& file)
  {

    if (idx) {
      std::cout << "...clearing old index" << std::endl;
      bwa_idx_destroy(idx);
      idx = 0;
    }
    
    // read in the bwa index
    std::cout << "...loading in the index for BWA from location: " << file << std::endl;
    idx = bwa_idx_load(file.c_str(), BWA_IDX_ALL);
    if (idx == NULL) {
      std::cerr << "Could not load the reference: " << file << std::endl;
      exit(EXIT_FAILURE);
    }
  }


  void BWAWrapper::writeIndexToFiles(const std::string& index_name)
  {
    
    if (!idx) {
      std::cerr << "Error: No index initialized. Can't write to file " << std::endl;
      return;
    }

    std::string bwt_name = index_name + ".bwt";
    std::string sa_name = index_name + ".sa";
    bwt_dump_bwt(bwt_name.c_str(), idx->bwt); 
    bwt_dump_sa(sa_name.c_str(), idx->bwt);
    bns_dump(idx->bns, index_name.c_str());
    __write_pac_to_file(index_name);
  }


  void BWAWrapper::__write_pac_to_file(const std::string& file)
  {
    // finalize .pac file
    FILE *fp;
    std::string nm = file + ".pac";
    fp = xopen(nm.c_str(), "wb");
    ubyte_t ct;
    err_fwrite(idx->pac, 1, (idx->bns->l_pac>>2) + ((idx->bns->l_pac&3) == 0? 0 : 1), fp);

    // the following codes make the pac file size always (l_pac/4+1+1)
    if (idx->bns->l_pac % 4 == 0) {
      ct = 0;
      err_fwrite(&ct, 1, 1, fp);
    }
    ct = idx->bns->l_pac % 4;
    err_fwrite(&ct, 1, 1, fp);

    // close .pac file
    err_fflush(fp);
    err_fclose(fp);
  }

  std::string BWAWrapper::__print_bns()
  {
    std::stringstream  ss;
    
    ss << "BNS: l_pac: " << idx->bns->l_pac << " n_seqs: " << idx->bns->n_seqs <<
      " seed: " << idx->bns->seed << " n_holes " << idx->bns->n_holes;
    return ss.str();
  }
}

