#include "SeqLib/BamRecord.h"

#include <cassert>
#include <bitset>
#include <cctype>

//#ifdef HAVE_BOOST
//#include <boost/algorithm/string.hpp>
//#endif

#include "SeqLib/ssw_cpp.h"

#define TAG_DELIMITER "^"
#define CTAG_DELIMITER '^'

namespace SeqLib {

  struct free_delete {
    void operator()(void* x) { bam_destroy1((bam1_t*)x); }
  };
  
  void BamRecord::init() {
    bam1_t* f = bam_init1();
    b = std::shared_ptr<bam1_t>(f, free_delete());
  }

  void BamRecord::assign(bam1_t* a) { 
    b = std::shared_ptr<bam1_t>(a, free_delete()); 
  }

  GenomicRegion BamRecord::asGenomicRegion() const {
    return GenomicRegion(b->core.tid, b->core.pos, PositionEnd());
  }

  GenomicRegion BamRecord::asGenomicRegionMate() const {
    return GenomicRegion(b->core.mtid, b->core.mpos, b->core.mpos + Length());
  }

  std::string BamRecord::Sequence() const {
    uint8_t * p = bam_get_seq(b);
    std::string out(b->core.l_qseq, 'N');
    for (int32_t i = 0; i < b->core.l_qseq; ++i) 
      out[i] = BASES[bam_seqi(p,i)];
    return out;
    
  }

  int32_t BamRecord::AlignmentLength() const {
    uint32_t* c = bam_get_cigar(b);
    int32_t len = 0;
    for (int k = 0; k < b->core.n_cigar; ++k) 
      if (bam_cigar_type(bam_cigar_op(c[k]))&1 || bam_cigar_opchr(c[k]) == 'H') // consumes query
	len += bam_cigar_oplen(c[k]);
    //cig.add(CigarField("MIDSSHP=XB"[c[k]&BAM_CIGAR_MASK], bam_cigar_oplen(c[k])));
    return len;
  }

  BamRecord::BamRecord(const std::string& name, const std::string& seq, const std::string& ref, const GenomicRegion * gr) {

    StripedSmithWaterman::Aligner aligner;
    // Declares a default filter
    StripedSmithWaterman::Filter filter;
    // Declares an alignment that stores the result
    StripedSmithWaterman::Alignment alignment;
    // Aligns the seq to the ref
    aligner.Align(seq.c_str(), ref.c_str(), ref.size(), filter, &alignment);

    init();
    b->core.tid = gr->chr;
    b->core.pos = gr->pos1 + alignment.ref_begin;
    b->core.qual = alignment.sw_score;
    b->core.flag = 0;
    b->core.n_cigar = alignment.cigar.size();
    
    // set dumy mate
    b->core.mtid = -1;
    b->core.mpos = -1;
    b->core.isize = 0;

    // allocate all the data
    b->core.l_qname = name.length() + 1;
    b->core.l_qseq = seq.length(); //(seq.length()>>1) + seq.length() % 2; // 4-bit encoding
    b->l_data = b->core.l_qname + (b->core.n_cigar<<2) + ((b->core.l_qseq+1)>>1) + (b->core.l_qseq);
    b.get()->data = (uint8_t*)malloc(b.get()->l_data);

    // allocate all the data
    b->core.l_qname = name.length() + 1;
    b->core.l_qseq = seq.length(); //(seq.length()>>1) + seq.length() % 2; // 4-bit encoding
    b->l_data = b->core.l_qname + (b->core.n_cigar<<2) + ((b->core.l_qseq+1)>>1) + (b->core.l_qseq);
    b.get()->data = (uint8_t*)malloc(b.get()->l_data);
    
    // allocate the qname
    memcpy(b->data, name.c_str(), name.length() + 1);

    // allocate the cigar. 32 bits per elem (4 type, 28 length)
    uint32_t * cigr = bam_get_cigar(b);
    for (size_t i = 0; i < alignment.cigar.size(); ++i)
      cigr[i] = alignment.cigar[i]; //Length << BAM_CIGAR_SHIFT | BAM_CMATCH;

    // allocate the sequence
    uint8_t* m_bases = b->data + b->core.l_qname + (b->core.n_cigar<<2);
    
    // TODO move this out of bigger loop
    int slen = seq.length();
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

  void BamRecord::SmartAddTag(const std::string& tag, const std::string& val)
  {
    // get the old tag
    assert(tag.length());
    assert(val.length());
    std::string tmp = GetZTag(tag);
    if (!tmp.length()) 
      {
	AddZTag(tag, val);
	return;
      }
    
    // check that we don't have the tag delimiter in the stirng
    if (val.find(TAG_DELIMITER) != std::string::npos)
      std::cerr << "BamRecord::SmartAddTag -- Tag delimiter " << TAG_DELIMITER << " is in the value to be added. Compile with diff tag delimiter or change val" << std::endl;

    // append the tag
    tmp += TAG_DELIMITER + val;
    
    // remove the old tag
    RemoveTag(tag.c_str());
    
    // add the new one
    assert(tmp.length());
    AddZTag(tag, tmp);
  }

  void BamRecord::clearSeqQualAndTags() {

    int new_size = b->core.l_qname + ((b)->core.n_cigar<<2);// + 1; ///* 0xff seq */ + 1 /* 0xff qual */;
    b->data = (uint8_t*)realloc(b->data, new_size);
    b->l_data = new_size;
    b->core.l_qseq = 0;
  }

  void BamRecord::SetSequence(const std::string& seq) {

    int new_size = b->l_data - ((b->core.l_qseq+1)>>1) - b->core.l_qseq + ((seq.length()+1)>>1) + seq.length();    
    int old_aux_spot = (b->core.n_cigar<<2) + b->core.l_qname + ((b->core.l_qseq + 1)>>1) + b->core.l_qseq;
    int old_aux_len = bam_get_l_aux(b); //(b->core.n_cigar<<2) + b->core.l_qname + ((b->core.l_qseq + 1)>>1) + b->core.l_qseq;

    // copy out all the old data
    uint8_t* oldd = (uint8_t*)malloc(b->l_data);
    memcpy(oldd, b->data, b->l_data);
    
    // clear out the old data and alloc the new amount
    free(b->data);
    b->data = (uint8_t*)calloc(new_size, sizeof(uint8_t)); 
    
    // add back the qname and cigar
    memcpy(b->data, oldd, b->core.l_qname + (b->core.n_cigar<<2));

    // update the sizes
    // >>1 shift is because only 4 bits needed per ATCGN base
    b->l_data = new_size; //b->l_data - ((b->core.l_qseq + 1)>>1) - b->core.l_qseq + ((seq.length()+1)>>1) + seq.length();
    b->core.l_qseq = seq.length();
    
    // allocate the sequence
    uint8_t* m_bases = b->data + b->core.l_qname + (b->core.n_cigar<<2);
    int slen = seq.length();

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

    // add in a NULL qual
    uint8_t* s = bam_get_qual(b);
    s[0] = 0xff;

    // add the aux data
    uint8_t* t = bam_get_aux(b);
    memcpy(t, oldd + old_aux_spot, old_aux_len);

    // reset the max size
    b->m_data = b->l_data;
    
  }
  
  void BamRecord::SetQname(const std::string& n)
  {
    // copy out the non-qname data
    size_t nonq_len = b->l_data - b->core.l_qname;
    uint8_t* nonq = (uint8_t*)malloc(nonq_len);
    memcpy(nonq, b->data + b->core.l_qname, nonq_len);

    // clear the old data and alloc the new amount 
    free(b->data);
    b->data = (uint8_t*)calloc(nonq_len + n.length() + 1, 1);
    
    // add in the new qnamev
    memcpy(b->data, (uint8_t*)n.c_str(), n.length() + 1); // +1 for \0

    // update the sizes
    b->l_data = b->l_data - b->core.l_qname + n.length() + 1;
    b->core.l_qname = n.length() + 1;    
    
    // copy over the old data
    memcpy(b->data + b->core.l_qname, nonq, nonq_len);
    free(nonq);

    // reset the max size
    b->m_data = b->l_data;
  }

  /*void BamRecord::SetSequence(std::string s)
  {

    // change the size to accomodate new sequence. Clear the quality string
    //std::cout << "osize " << b->l_data << " calcsize " << (b->core.l_qseq + b->core.l_qname + (b->core.n_cigar<<2)) << std::endl;
    b->data = (uint8_t*)realloc(b->data, b->core.l_qname + s.length() + (b->core.n_cigar<<2));
    
    // copy in the new sequence
    memcpy(b->data + b->core.l_qname + (b->core.n_cigar<<2), (uint8_t*)s.c_str(), s.length());

    }*/

  double BamRecord::MeanPhred() const {

    if (b->core.l_qseq <= 0)
      return -1;

    double s = 0;
    uint8_t* p = bam_get_qual(b);
    for (int32_t i = 0; i < b->core.l_qseq; ++i)
      s += p[i];
    return s / b->core.l_qseq;
  }

  std::string BamRecord::QualitySequence() const {
    std::string seq = GetZTag("GV");
    if (!seq.length()) 
      seq = Sequence();
    return seq;
  }

  std::ostream& operator<<(std::ostream& out, const BamRecord &r)
  {
    if (!r.b) {
      out << "empty read";
      return out;
    }

    out << bam_get_qname(r.b) << "\t" << r.b->core.flag
	<< "\t" << (r.b->core.tid+1) << "\t" << r.b->core.pos 
	<< "\t" << r.b->core.qual << "\t" << r.CigarString() 
	<< "\t" << (r.b->core.mtid+1) << "\t" << r.b->core.mpos << "\t" 
        << r.FullInsertSize() //r.b->core.isize 
	<< "\t" << r.Sequence() << "\t*" << 
      "\tAS:" << r.GetIntTag("AS") << 
      "\tDD:" << r.GetIntTag("DD");/* << "\t" << r.Qualities()*/;;/* << "\t" << r.Qualities()*/;
    return out;
      
    
  }

  int32_t BamRecord::CountSecondaryAlignments() const 
  {
    int xp_count = 0;
    
    // xa tag
    std::string xar_s = GetZTag("XA");
    //r_get_Z_tag(r, "XA", xar_s);
    if (xar_s.length()) {
      xp_count += std::count(xar_s.begin(), xar_s.end(), ';');
    }
    
    // xp tag
    std::string xpr_s = GetZTag("XP");
    //r_get_Z_tag(r, "XP", xpr_s);
    if (xpr_s.length()) {
      xp_count += std::count(xpr_s.begin(), xpr_s.end(), ';');
    }

    return xp_count;
    
  }

  int32_t BamRecord::CountNBases() const {
    uint8_t* p = bam_get_seq(b); 
    int32_t n = 0;
    for (int ww = 0; ww < b->core.l_qseq; ww++)
      if (bam_seqi(p,ww) == 15) 
	++n; 
    return n;
  }

  void BamRecord::QualityTrimmedSequence(int32_t qualTrim, int32_t& startpoint, int32_t& endpoint) const {

    endpoint = -1; //seq.length();
    startpoint = 0;
    int i = 0; 
    
    uint8_t * qual = bam_get_qual(b.get());
    
    // if there is no quality score, return whole thing
    if (qual[0] == 0xff) {
      startpoint = 0;
      return;
      
      //return Sequence();
    }
    
    // get the start point (loop forward)
    while(i < b->core.l_qseq) {
      int ps = qual[i];
      if (ps >= qualTrim) {
          startpoint = i;
          break;
	}
	++i;
    }

    // get the end point (loop backwards)
    i = b->core.l_qseq - 1; //seq.length() - 1;
    while(i >= 0) {

      int ps = qual[i];
      
      if (ps >= qualTrim) { //ps >= qualTrim) {
	endpoint = i + 1; // endpoint is one past edge
	break;
      }
      --i;
    }

    /*
    // check that they aren't all bad
    if (startpoint == 0 && endpoint == -1) 
      return;

    // if they're all good, mark and return
    if (startpoint == 0 && endpoint == b->core.l_qseq) {
      return; 
    }
    has_trim = true;
    
    std::string output = std::string(endpoint-startpoint, 'N');
    try { 
      uint8_t * p = bam_get_seq(b);
      for (int32_t i = startpoint; i < (endpoint - startpoint); ++i) 
	output[i] = BASES[bam_seqi(p,i)];
    } catch (...) {
      std::cerr << "Trying to subset string in BamRecord::QualityTrimRead out of bounds. String: " << Sequence() << " start " << startpoint << " length " << (endpoint - startpoint) << std::endl;
    }

    return output;
    */
  }

  void BamRecord::AddZTag(std::string tag, std::string val) {
    if (tag.empty() || val.empty())
      return;
    bam_aux_append(b.get(), tag.data(), 'Z', val.length()+1, (uint8_t*)val.c_str());
  }

  std::string BamRecord::GetZTag(const std::string& tag) const {
    uint8_t* p = bam_aux_get(b.get(),tag.c_str());
    if (!p)
      return std::string();
    char* pp = bam_aux2Z(p);
    if (!pp) 
      return std::string();
    return std::string(pp);
  }

  
  // get a string tag that might be separted by "x"
  std::vector<std::string> BamRecord::GetSmartStringTag(const std::string& tag) const {
    
    std::vector<std::string> out;
    std::string tmp = GetZTag(tag);

    if (tmp.empty())
      return std::vector<std::string>();
    
    if (tmp.find(TAG_DELIMITER) != std::string::npos) {
      std::istringstream iss(tmp);
      std::string line;
      while (std::getline(iss, line, CTAG_DELIMITER)) {
	out.push_back(line);
      }
    } else {
      out.push_back(tmp);
    }
    
    assert(out.size());
    return out;
    
  }
  
  
  std::vector<int> BamRecord::GetSmartIntTag(const std::string& tag) const {
    
    std::vector<int> out;
    std::string tmp;
    
    tmp = GetZTag(tag);
    //r_get_Z_tag(a, tag.c_str(), tmp);
    //assert(tmp.length());
    if (tmp.empty())
      return std::vector<int>();
    
    if (tmp.find(TAG_DELIMITER) != std::string::npos) {
      std::istringstream iss(tmp);
      std::string line;
      while (std::getline(iss, line, CTAG_DELIMITER)) {
	try { out.push_back(stoi(line)); } catch (...) { std::cerr << "Failed to read parsed int tag " << tag << " for value " << tmp << " with line " << line << std::endl; std::exit(EXIT_FAILURE); }
      }
    } else {
      try { out.push_back(stoi(tmp)); } catch (...) { std::cerr << "Failed to read int tag " << tag << " for value " << tmp << std::endl; std::exit(EXIT_FAILURE); }
    }
    
    assert(out.size());
    return out;
    
  }

  bool BamRecord::coveredBase(int pos) const {

    if (pos < 0) 
      return false;

    if (NumClip() == 0) 
      return true;
    
    Cigar cig = GetCigar();
    assert(cig.size() > 1); // are clips, so has to be at least two fields

    // get length of read (whether hardclipped or not)
    int len = 0;
    for (auto& c : cig) 
      if (c.ConsumesQuery())
	len += c.Length(); 

    int lbound  = 0, posr = 0, rbound = len - 1;

    // progress the left side
    for (auto& c : cig) {
      if (c.ConsumesQuery()) {
	lbound = posr;
	break;
      }
      posr += c.Length();
    }

    // progress the right side
    posr = Length() - 1;
    for (auto& c : cig) {
      if (c.ConsumesQuery()) {
	rbound = posr;
	break;
      }
      posr -= c.Length();
    }

    return pos >= lbound && pos <= rbound;
  }

  bool BamRecord::coveredMatchBase(int pos) const {

    if (pos < 0) 
      return false;

    if (NumClip() == 0) 
      return true;
    
    Cigar cig = GetCigar();
    assert(cig.size() > 1); // are clips, so has to be at least two fields

    // get length of read (whether hardclipped or not)
    int len = 0;
    for (auto& c : cig) 
      if (c.ConsumesQuery())
	len += c.Length(); 

    int lbound  = 0, posr = 0, rbound = len - 1;

    // progress the left side
    for (auto& c : cig) {
      if (c.Type() == 'M') {
	lbound = posr;
	break;
      }
      posr += c.Length();
    }

    // progress the right side
    posr = Length() - 1;
    for (auto& c : cig) {
      if (c.Type() == 'M') {
	rbound = posr;
	break;
      }
      posr -= c.Length();
    }

    return pos >= lbound && pos <= rbound;
  }

  
  BamRecord::BamRecord(const std::string& name, const std::string& seq, const GenomicRegion * gr, const Cigar& cig) {

    // make sure cigar fits with sequence
    size_t clen = 0;
    for (auto& i : cig)
      if (i.ConsumesQuery())
	//if (i.RawType() == BAM_CMATCH || i.RawType() == BAM_CSOFT_CLIP || i.RawType() == BAM_CINS)
	clen += i.Length();
    assert(seq.length() == clen);

    init();
    b->core.tid = gr->chr;
    b->core.pos = gr->pos1 + 1;
    b->core.qual = 60;
    b->core.flag = 0;
    b->core.n_cigar = cig.size();
    
    // set dumy mate
    b->core.mtid = -1;
    b->core.mpos = -1;
    b->core.isize = 0;
      
    // if alignment is reverse, set it
    if (gr->strand == '-') // just choose this convention to reverse
      b->core.flag |= BAM_FREVERSE;
    
    // allocate all the data
    b->core.l_qname = name.length() + 1;
    b->core.l_qseq = seq.length(); //(seq.length()>>1) + seq.length() % 2; // 4-bit encoding
    b->l_data = b->core.l_qname + (b->core.n_cigar<<2) + ((b->core.l_qseq+1)>>1) + (b->core.l_qseq);
    b.get()->data = (uint8_t*)malloc(b.get()->l_data);
    
    // allocate all the data
    b->core.l_qname = name.length() + 1;
    b->core.l_qseq = seq.length(); //(seq.length()>>1) + seq.length() % 2; // 4-bit encoding
    b->l_data = b->core.l_qname + (b->core.n_cigar<<2) + ((b->core.l_qseq+1)>>1) + (b->core.l_qseq);
    b.get()->data = (uint8_t*)malloc(b.get()->l_data);
    
    // allocate the qname
    memcpy(b->data, name.c_str(), name.length() + 1);
      
    // allocate the cigar. 32 bits per elem (4 type, 28 length)
    uint32_t * cigr = bam_get_cigar(b);
    for (size_t i = 0; i < cig.size(); ++i)
      cigr[i] = cig[i].raw(); //Length << BAM_CIGAR_SHIFT | BAM_CMATCH;
    
    // allocate the sequence
    uint8_t* m_bases = b->data + b->core.l_qname + (b->core.n_cigar<<2);

      // TODO move this out of bigger loop
      int slen = seq.length();
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
  

  CigarField::CigarField(char  t, uint32_t len) {
    
    int op = 0;

    // should use a table for this
    if (t == 'M')
      op = BAM_CMATCH;
    else if (t == 'D')
      op = BAM_CDEL;
    else if (t == 'I')
      op = BAM_CINS;
    else if (t == 'S')
      op = BAM_CSOFT_CLIP;
    else if (t == 'N')
      op = BAM_CREF_SKIP;
    else
      assert(false);
    data = len << BAM_CIGAR_SHIFT | op;

  }

  std::ostream& operator<<(std::ostream& out, const CigarField& c) { 
    out << bam_cigar_oplen(c.data) << bam_cigar_opchr(c.data); 
    return out; 
  }


  std::ostream& operator<<(std::ostream& out, const Cigar& c) { 
    for (auto& i : c)
      out << i;
    return out; 
  }


  Cigar cigarFromString(const std::string& cig) {

    Cigar tc;


    // get the ops (MIDSHPN)
    std::vector<char> ops;
    for (auto& c : cig)
      if (!isdigit(c)) {
	ops.push_back(c);
      }

    std::size_t prev = 0, pos;
    std::vector<std::string> lens;
    while ((pos = cig.find_first_of("MIDSHPNX", prev)) != std::string::npos) {
        if (pos > prev)
	  lens.push_back(cig.substr(prev, pos-prev));
        prev = pos+1;
      }
    if (prev < cig.length())
      lens.push_back(cig.substr(prev, std::string::npos));

    assert(ops.size() == lens.size());
    for (size_t i = 0; i < lens.size(); ++i) {
      //tc.push_back(CigarField(ops[i], std::stoi(lens[i])));
      tc.add(CigarField(ops[i], std::stoi(lens[i])));
    }
    
    return tc;


/*
#ifdef HAVE_BOOST
    std::vector<std::string> lens;
    boost::split(lens, cig, boost::is_any_of(cigar_delimiters));
    lens.pop_back(); // fills in empty at end for some reason

#else

    std::cerr << "!!!!!! BOOST REQUIRED FOR THIS FUNCTION !!!!!!" << std::endl << 
      "!!!!!! Returning EMPTY CIGAR !!!!!!!!!" << std::endl;
    return tc;
#endif
*/

  }
  
}
