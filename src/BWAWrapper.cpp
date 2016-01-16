#include "SnowTools/BWAWrapper.h"

#include <sstream>
#include <iostream>

extern "C" {
  #include <string.h>
}

//#define DEBUG_BWATOOLS 1

#define _set_pac(pac, l, c) ((pac)[(l)>>2] |= (c)<<((~(l)&3)<<1))
#define _get_pac(pac, l) ((pac)[(l)>>2]>>((~(l)&3)<<1)&3)

namespace SnowTools {

  int BWAWrapper::refCount() const {
    
    if (!idx)
      return 0;

    return idx->bns->n_seqs;
    
  }

  std::string BWAWrapper::ChrIDToName(int id) const {
    assert(idx);
    assert(id >= 0);
    assert(id < idx->bns->n_seqs);
    return std::string(idx->bns->anns[id].name);
  }

  bam_hdr_t * BWAWrapper::HeaderFromIndex() const 
  {

    std::string my_hdr = bwa_print_sam_hdr2(idx->bns, "");
    bam_hdr_t * hdr = bam_hdr_init();
    hdr = sam_hdr_read2(my_hdr); 
    //hdr->n_targets = idx->bns->n_seqs;
    //hdr->target_name = (char**)malloc(hdr->n_targets * sizeof(char*));
    //for (int i = 0; i < idx->bns->n_seqs; ++i) {
    //  hdr->target_name[i] = (char*)malloc( (strlen(idx->bns->anns[i].name) + 1) * sizeof(char));
    //  strcpy(hdr->target_name[i], idx->bns->anns[i].name);
    //}
    return hdr;
  }

  std::string BWAWrapper::bwa_print_sam_hdr2(const bntseq_t *bns, const char *hdr_line) const
  {
    std::string out;
    int i, n_SQ = 0;
    //extern char *bwa_pg;
    if (hdr_line) {
      const char *p = hdr_line;
      while ((p = strstr(p, "@SQ\t")) != 0) {
	if (p == hdr_line || *(p-1) == '\n') ++n_SQ;
	p += 4;
      }
    }

    //JEREMIAH
    // get makx size
    size_t max_s = 0;
    for (i = 0; i < bns->n_seqs; ++i)
      max_s = std::max(strlen(bns->anns[i].name), max_s);

    if (n_SQ == 0) {
      char buffer[max_s + 30];
      for (i = 0; i < bns->n_seqs; ++i) {
	//err_printf("@SQ\tSN:%s\tLN:%d\n", bns->anns[i].name, bns->anns[i].len);
	sprintf(buffer, "@SQ\tSN:%s\tLN:%d\n", bns->anns[i].name, bns->anns[i].len);
	out.append(buffer);
      }
    } else if (n_SQ != bns->n_seqs && bwa_verbose >= 2)
      fprintf(stderr, "[W::%s] %d @SQ lines provided with -H; %d sequences in the index. Continue anyway.\n", __func__, n_SQ, bns->n_seqs);
    
    if (hdr_line) { char buffer[200]; sprintf(buffer, "%s\n", hdr_line); out.append(buffer); } //err_printf("%s\n", hdr_line);
    //if (bwa_pg) { char buffer[100]; sprintf(buffer, "%s\n", bwa_pg); out.append(buffer); } // err_printf("%s\n", bwa_pg);
    
    return out;
  }
  
  bam_hdr_t* BWAWrapper::sam_hdr_read2(const std::string& hdr) const {
    kstring_t str;
    bam_hdr_t *h;
    str.l = str.m = 0; str.s = 0;
    
    std::istringstream iss(hdr);
    std::string line;
    while (std::getline(iss, line, '\n')) {
      //while (hts_getline(fp, KS_SEP_LINE, &fp->line) >= 0) {
      if (line.length() == 0 || line.at(0) != '@') break;
      
      //if (line.length() > 3 && line.substr(0,3) == "@SQ") has_SQ = 1;
      //if (fp->line.l > 3 && strncmp(fp->line.s,"@SQ",3) == 0) has_SQ = 1;
      //kputsn(fp->line.s, fp->line.l, &str);
      kputsn(line.c_str(), line.length(), &str);
      kputc('\n', &str);
    }
    /*
      if (! has_SQ && fp->fn_aux) {
      char line[2048];
      FILE *f = fopen(fp->fn_aux, "r");
      if (f == NULL) return NULL;
      while (fgets(line, sizeof line, f)) {
      const char *name = strtok(line, "\t");
      const char *length = strtok(NULL, "\t");
      ksprintf(&str, "@SQ\tSN:%s\tLN:%s\n", name, length);
      }
      fclose(f);
      }
    */
    if (str.l == 0) kputsn("", 0, &str);
    h = sam_hdr_parse(str.l, str.s);
    h->l_text = str.l; h->text = str.s;
    return h;
  }
  
  
  void BWAWrapper::constructIndex(const USeqVector& v) {

    if (!v.size())
      return;

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
    std::cerr << "ref seq length: " << tlen << std::endl;
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

  void BWAWrapper::alignSingleSequence(const std::string& seq, const std::string& name, BamReadVector& vec, 
				       double keep_sec_with_frac_of_primary_score, int max_secondary) {
    mem_alnreg_v ar;
    ar = mem_align1(memopt, idx->bwt, idx->bns, idx->pac, seq.length(), seq.c_str()); // get all the hits

#ifdef DEBUG_BWATOOLS
    std::cout << "num hits: " << ar.n << std::endl;
    //std::cout << __print_bns() << std::endl;
#endif    

    double primary_score = 0;
    //size_t num_secondary = 0;
    // loop through the hits
    for (size_t i = 0; i < ar.n; ++i) {

      if (ar.a[i].secondary >= 0 && (keep_sec_with_frac_of_primary_score < 0 || keep_sec_with_frac_of_primary_score > 1))
      	continue; // skip secondary alignments
      
      // get forward-strand position and CIGAR
      mem_aln_t a;
#ifdef DEBUG_BWATOOLS
      std::cerr << "awdlfkajsdf" << std::endl;
#endif
      a = mem_reg2aln(memopt, idx->bns, idx->pac, seq.length(), seq.c_str(), &ar.a[i]); 

      //if (name == "c_1_1453100_1473100_12")
      //std::cerr << name << " secondary " << ar.a[i].secondary << " primary_score " << primary_score << " a.score " << ar.a[i].score << " sub_N " << ar.a[i].sub_n << 
      //	" frac_rep " << ar.a[i].frac_rep << " flag " << a.flag << std::endl;

#ifdef DEBUG_BWATOOLS
      std::cerr << "DDDDDDDDDDDDDawdlfkajsdf" << std::endl;
#endif

      // if score not sufficient or past cap, continue
      bool sec_and_low_score =  ar.a[i].secondary >= 0 && (primary_score * keep_sec_with_frac_of_primary_score) > a.score;
      bool sec_and_cap_hit = ar.a[i].secondary >= 0 && i > max_secondary;
      if (sec_and_low_score || sec_and_cap_hit) {

#ifdef DEBUG_BWATOOLS
      std::cerr << "freeing cigar" << std::endl;
#endif

	free(a.cigar);
	continue;
      } else if (ar.a[i].secondary < 0) {
	primary_score = a.score;
	//num_secondary = 0;
      }

      //if (i == 0 && a.is_rev)
      //first_is_rev = true;

#ifdef DEBUG_BWATOOLS
      std::cerr << "allocing bamread" << std::endl;
#endif

      // instantiate the read
      BamRead b;
      b.init();

      b.b->core.tid = a.rid;
      b.b->core.pos = a.pos;
      b.b->core.qual = a.mapq;
      b.b->core.flag = a.flag;
      b.b->core.n_cigar = a.n_cigar;
      
      // set dumy mate
      b.b->core.mtid = -1;
      b.b->core.mpos = -1;
      b.b->core.isize = 0;

      // if alignment is reverse, set it
      if (a.is_rev) 
	b.b->core.flag |= BAM_FREVERSE;

      // allocate all the data
      b.b->core.l_qname = name.length() + 1;
      b.b->core.l_qseq = seq.length(); //(seq.length()>>1) + seq.length() % 2; // 4-bit encoding
      b.b->l_data = b.b->core.l_qname + (a.n_cigar<<2) + ((b.b->core.l_qseq+1)>>1) + (b.b->core.l_qseq);
      b.b.get()->data = (uint8_t*)malloc(b.b.get()->l_data);

#ifdef DEBUG_BWATOOLS
      std::cerr << "memcpy" << std::endl;
#endif

      // allocate the qname
      memcpy(b.b->data, name.c_str(), name.length() + 1);

      // allocate the cigar. Reverse if aligned to neg strand, since mem_aln_t stores
      // cigars relative to referemce string oreiatnion, not forward alignment
      memcpy(b.b->data + b.b->core.l_qname, (uint8_t*)a.cigar, a.n_cigar<<2);
      //std::cerr << "ORIGINAL CIGAR "  << b.CigarString() << " rev " << a.is_rev << std::endl;

#ifdef DEBUG_BWATOOLS
      std::cerr << "memcpy2" << std::endl;
#endif

      /*      if (a.is_rev != first_is_rev) {
	uint32_t temp;
	int start = 0; 
	int end = a.n_cigar - 1;
	uint32_t* arr = bam_get_cigar(b.b);
	while(start < end) {
	    temp = arr[start];   
	    arr[start] = arr[end];
	    arr[end] = temp;
	    ++start;
	    --end;
	  }  
	  }*/

      //std::cerr << "NEW CIGAR "  << b.CigarString() << " rev " << a.is_rev << std::endl;
      
      // convert N to S
      uint32_t * cigr = bam_get_cigar(b.b);
      for (int k = 0; k < b.b->core.n_cigar; ++k) {
	if ( (cigr[k] & BAM_CIGAR_MASK) == BAM_CREF_SKIP) {
	  cigr[k] &= ~BAM_CIGAR_MASK;
	  cigr[k] |= BAM_CSOFT_CLIP;
	}
      }
	
      // allocate the sequence
      uint8_t* m_bases = b.b->data + b.b->core.l_qname + (b.b->core.n_cigar<<2);

      // TODO move this out of bigger loop
      int slen = seq.length();
      int j = 0;
      if (a.is_rev/* && false*/) {
	for (int i = slen-1; i >= 0; --i) {
	  
	  // bad idea but works for now
	  // this is REV COMP things
	  uint8_t base = 15;
	  if (seq.at(i) == 'T')
	    base = 1;
	  else if (seq.at(i) == 'G')
	    base = 2;
	  else if (seq.at(i) == 'C')
	    base = 4;
	  else if (seq.at(i) == 'A')
	    base = 8;

	  m_bases[j >> 1] &= ~(0xF << ((~j & 1) << 2));   ///< zero out previous 4-bit base encoding
	  m_bases[j >> 1] |= base << ((~j & 1) << 2);  ///< insert new 4-bit base encoding
	  ++j;
	}
      } else {
	for (int i = 0; i < slen; ++i) {
	// bad idea but works for now
	  uint8_t base = 15;
	  if (seq.at(i) == 'A')
	    base = 1;
	  else if (seq.at(i) == 'C')
	    base = 2;
	  else if (seq.at(i) == 'G')
	    base = 4;
	  else if (seq.at(i) == 'T')
	    base = 8;
	  
	  m_bases[i >> 1] &= ~(0xF << ((~i & 1) << 2));   ///< zero out previous 4-bit base encoding
	  m_bases[i >> 1] |= base << ((~i & 1) << 2);  ///< insert new 4-bit base encoding

	}
      }

#ifdef DEBUG_BWATOOLS
      std::cerr << "memcpy3" << std::endl;
#endif

      // allocate the quality to NULL
      uint8_t* s = bam_get_qual(b.b);
      s[0] = 0xff;

      //b.SetQname(name);
      //b.SetSequence(seq);

      b.AddIntTag("NA", ar.n); // number of matches
      //std::cerr << a.NM << std::endl;
      b.AddIntTag("NM", a.NM);

      if (a.XA)
	b.AddZTag("XA", std::string(a.XA));

      // add num sub opt
      b.AddIntTag("SB", ar.a[i].sub_n);
      b.AddIntTag("AS", a.score);

      vec.push_back(b);

      //#ifdef DEBUG_BWATOOLS
      // print alignment
      //printf("\t%c\t%s\t%ld\t%d\t", /*ks->name.s,*/ "+-"[a.is_rev], idx->bns->anns[a.rid].name, (long)a.pos, a.mapq);
      //for (int k = 0; k < a.n_cigar; ++k) // print CIGAR
      //	printf("%d%c", a.cigar[k]>>4, "MIDSH"[a.cigar[k]&0xf]);
      // printf("\t%d\n", a.NM); // print edit distance
      //#endif
#ifdef DEBUG_BWATOOLS
      std::cerr << "final done" << std::endl;
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
    name->l = v[k].name.length() + 1;
    name->m = v[k].name.length() + 3;
    name->s = (char*)calloc(name->m, sizeof(char));
    memcpy(name->s, v[k].name.c_str(), v[k].name.length()+1);
    
    // make the sequence kstring
    kstring_t * t = (kstring_t*)malloc(sizeof(kstring_t));
    t->l = v[k].seq.length();
    t->m = v[k].seq.length() + 2;
    //t->s = (char*)calloc(v[k].seq.length(), sizeof(char));
    t->s = (char*)malloc(t->m);
    memcpy(t->s, v[k].seq.c_str(), v[k].seq.length());
    
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
      std::cerr << "...clearing old index" << std::endl;
      bwa_idx_destroy(idx);
      idx = 0;
    }
    
    // read in the bwa index
    std::cerr << "...loading in the index for BWA from location: " << file << std::endl;
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

