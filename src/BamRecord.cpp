#include "SeqLib/BamRecord.h"

#include <tuple>
#include <cassert>
#include <bitset>
#include <cctype>
#include <sstream>
#include <stdexcept>
#include <numeric>
#include <regex>
#include <algorithm>
#include <cctype>

//inline constexpr char CTAG_DELIMITER = '^';
//inline constexpr std::string_view TAG_DELIMITER = "^";

namespace SeqLib {

  const int CigarCharToInt[128] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, //0-9
				   -1,-1,-1,-1,-1,-1,-1,-1,-1,-1, //10-19
				   -1,-1,-1,-1,-1,-1,-1,-1,-1,-1, //20
				   -1,-1,-1,-1,-1,-1,-1,-1,-1,-1, //30
				   -1,-1,-1,-1,-1,-1,-1,-1,-1,-1, //40
				   -1,-1,-1,-1,-1,-1,-1,-1,-1,-1, //50
				   -1,BAM_CEQUAL,-1,-1,-1,-1,BAM_CBACK,-1,BAM_CDEL,-1, //60-69
				   -1,-1,BAM_CHARD_CLIP,BAM_CINS,-1,-1,-1,BAM_CMATCH,BAM_CREF_SKIP,-1,
				   BAM_CPAD,-1,-1,BAM_CSOFT_CLIP,-1,-1,-1,-1,BAM_CDIFF,-1,
				   -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
				   -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
				   -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
				   -1,-1,-1,-1,-1,-1,-1,-1};
  
  CigarField::CigarField(char  opChr, uint32_t len) {
    int op = CigarCharToInt[static_cast<unsigned char>(opChr)];
    if (op < 0)
      throw std::invalid_argument("Cigar type must be one of MIDSHPN=X");      
    data = len << BAM_CIGAR_SHIFT;
    data = data | static_cast<uint32_t>(op);
  }
  
  /** Return the number of query-consumed bases */
  int Cigar::NumQueryConsumed() const noexcept {
    int total = 0;
    for (const auto& cf : m_data) {
      if (cf.ConsumesQuery()) {
	total += static_cast<int>(cf.Length());
      }
    }
    return total;
  }
  
  /** Return the number of reference-consumed bases */
  int Cigar::NumReferenceConsumed() const noexcept {
    int total = 0;
    for (const auto& cf : m_data) {
      if (cf.ConsumesReference()) {
	total += static_cast<int>(cf.Length());
      }
    }
    return total;
  }
  
  /** Add a new cigar op */
  void Cigar::add(const CigarField& c) { 
    m_data.push_back(c); 
  }
  
  std::ostream& operator<<(std::ostream& os, const CigarField& cf) noexcept {
    // print length then opchar
    return os << cf.Length() << bam_cigar_opchr(cf.data);
  }
  
  std::ostream& operator<<(std::ostream& os, const Cigar& c) noexcept {
    for (const auto& cf : c) 
      os << cf;
    return os;
  }

  bool Cigar::operator==(const Cigar& o) const noexcept {
    return m_data == o.m_data;
  }

  Cigar::Cigar(const std::string& cig) {
    static const std::regex token{R"((\d+)([MIDNSHP=X]))"};
    std::smatch m;
    auto s = cig;
    while (std::regex_search(s, m, token)) {
      int len = std::stoi(m[1].str());
      char op  = m[2].str()[0];
      add(CigarField(op, len));
      s = m.suffix().str();
    }
    assert(size() == 0 || // empty is OK
           std::all_of(begin(), end(), [](const CigarField& cf){
             return cf.Length() > 0;
           }));
  }
  
  BamRecord::BamRecord() {
    this->init_();
  }
  
  void BamRecord::init_() {
    bam1_t* f = bam_init1();
    b = SeqPointer<bam1_t>(f, Bam1Deleter());
  }

  BamRecord::BamRecord(std::string_view qname,
		       std::string_view seq,
		       const GenomicRegion& gr,
		       const Cigar& cig) 
  {
    // Validate lengths up front
    if (cig.NumQueryConsumed() != seq.size())
      throw std::invalid_argument("Sequence length mismatches CIGAR query consumption");
    if (cig.NumReferenceConsumed() != gr.Width())
      throw std::invalid_argument("GenomicRegion width mismatches CIGAR reference consumption");
    
    // Initialize an empty bam1_t
    init_();
    
    // Alias core & b->data for brevity
    auto& core = b->core;
    
    // Set up the basics
    core.tid     = gr.chr;
    core.pos     = gr.pos1;
    core.qual    = 60;
    core.flag    = (gr.strand == '-') ? (BAM_FREVERSE) : 0;
    core.n_cigar = static_cast<uint32_t>(cig.size());
    core.mtid    = -1;
    core.mpos    = -1;
    core.isize   =  0;
    
    // Compute sizes
    const size_t nameBytes = qname.size() + 1;            // include null
    const size_t cigarBytes = cig.size() * sizeof(uint32_t);
    const size_t seq4bitBytes = ((seq.size() + 1) / 2);    // 4-bit packing
    const size_t seqQualBytes = seq.size();               // 1 byte per base quality
    const size_t totalData = nameBytes + cigarBytes + seq4bitBytes + seqQualBytes;
    
    // Allocate data buffer
    core.l_qname = static_cast<uint32_t>(nameBytes);
    core.l_qseq  = static_cast<uint32_t>(seq.size());
    b->l_data    = static_cast<int>(totalData);
    b->m_data    = b->l_data;
    b->data      = static_cast<uint8_t*>(std::malloc(totalData));
    
    // 1) copy in the qname (with terminal '\0')
    std::memcpy(b->data, qname.data(), nameBytes);
    
    // 2) write the CIGAR
    auto* cigPtr = bam_get_cigar(b.get());
    for (size_t i = 0; i < cig.size(); ++i)
      cigPtr[i] = cig[i].raw();
    
    // 3) encode the sequence in 4-bit
    uint8_t* seqPtr = b->data + nameBytes + cigarBytes;
    std::fill_n(seqPtr, seq4bitBytes, 0);
    for (size_t i = 0; i < seq.size(); ++i) {
      uint8_t code = 15;
      switch (std::toupper(seq[i])) {
      case 'A': code = 1; break;
      case 'C': code = 2; break;
      case 'G': code = 4; break;
      case 'T': code = 8; break;
      }
      size_t idx   = i >> 1;
      int    shift = ((~i & 1) << 2);
      seqPtr[idx] = (seqPtr[idx] & ~(0xFu << shift)) | (code << shift);
    }
    
    // 4) set a NULL quality (0xff)
    uint8_t* qualPtr = bam_get_qual(b.get());
    qualPtr[0] = 0xff;
    
    // any remaining aux-tags are uninitialized (e.g. you can call AddZTag etc.)
}
  
  // struct free_delete {
  //   void operator()(void* x) { bam_destroy1((bam1_t*)x); }
  // };
  
  // void BamRecord::init() {
  //   bam1_t* f = bam_init1();
  //   b = SeqPointer<bam1_t>(f, Bam1Deleter());
  // }

  //  void BamRecord::assign(bam1_t* a) { 
  //   b = SeqPointer<bam1_t>(a, Bam1Deleter()); 
  //}

  // int32_t BamRecord::PositionWithSClips() const {
  //   if(!b) return -1; // to be consistent with BamRecord::Position()
    
  //   uint32_t* cig = bam_get_cigar(b);
  //   return ((*cig) & 0xF) == BAM_CSOFT_CLIP ? b->core.pos - ((*cig) >> 4) : b->core.pos;
  // }
  
  int32_t BamRecord::PositionEnd() const { 
    return b ? (b->core.l_qseq > 0 ? bam_endpos(b.get()) : b->core.pos + GetCigar().NumQueryConsumed()) : -1;
  }

  // int32_t BamRecord::PositionEndWithSClips() const {
  //   if(!b) return -1; // to be consistent with BamRecord::PositionEnd()

  //   uint32_t* cig_last = bam_get_cigar(b) + b->core.n_cigar - 1;
  //   if(b->core.l_qseq > 0) {
  //     return ((*cig_last) & 0xF) == BAM_CSOFT_CLIP ? bam_endpos(b.get()) + ((*cig_last) >> 4) :
  //                                                    bam_endpos(b.get());
  //   } else {
  //     return b->core.pos + GetCigar().NumQueryConsumed();
  //   }
  // }

  int32_t BamRecord::PositionEndMate() const {
    if (!b) 
      return -1;
    
    // figure out how many query bases the mate consumes
    // (either the flagged read length or, if missing, the sum of CIGAR ops)
    int32_t queryLen = b->core.l_qseq > 0
      ? static_cast<int32_t>(b->core.l_qseq)
      : GetCigar().NumQueryConsumed();
    
    // and return its end-position
    return b->core.mpos + queryLen;
  }

  GenomicRegion BamRecord::AsGenomicRegion() const {
    char s = '*';
    if (MappedFlag())
      s = ReverseFlag() ? '-' : '+';
    return GenomicRegion(
			 b->core.tid,
			 static_cast<int32_t>(b->core.pos),
			 PositionEnd(),
			 s
			 );
  }
  
  GenomicRegion BamRecord::AsGenomicRegionMate() const {
    char s = '*';
    if (MateMappedFlag())
      s = MateReverseFlag() ? '-' : '+';
    return GenomicRegion(
      b->core.mtid,
      static_cast<int32_t>(b->core.mpos),
      PositionEndMate(),
      s
			 );
  }


  std::string BamRecord::Sequence() const {
    // if there's no record loaded, return an empty string
    if (!b)
      return {};
    
    // grab the encoded 4-bit sequence
    uint8_t* seq_enc = bam_get_seq(b.get());
    size_t len     = b->core.l_qseq;
    
    // we know exactly how many bases we'll emit
    std::string out;
    out.reserve(len);
    
    // decode one base at a time
    for (size_t i = 0; i < len; ++i) {
      out.push_back(BASES[bam_seqi(seq_enc, i)]);
    }
    
    return out;
  }

  void BamRecord::SetCigar(const Cigar& c) {
    if (!b) return;  // no record loaded
    
    const size_t oldNCigar = b->core.n_cigar;
    const size_t newNCigar = c.size();
    const size_t qNameBytes = b->core.l_qname;
    const size_t oldCigarBytes = oldNCigar * sizeof(uint32_t);
    const size_t newCigarBytes = newNCigar * sizeof(uint32_t);
    const size_t oldDataLen = b->l_data;
    const size_t auxOffset = qNameBytes + oldCigarBytes;
    const size_t auxLen = bam_get_l_aux(b.get())
      + ((b->core.l_qseq + 1) / 2)
      + b->core.l_qseq;
    
    // Case 1: same number of ops just overwrite in-place
    if (newNCigar == oldNCigar) {
      uint32_t* dst = bam_get_cigar(b.get());
      std::transform(
		     c.begin(), c.end(), dst,
		     [](const CigarField& cf){ return cf.raw(); }
		     );
      return;
    }
    
    // Case 2: different size rebuild data buffer
    size_t newDataLen = oldDataLen - oldCigarBytes + newCigarBytes;
    
    // Backup old raw bytes
    std::vector<uint8_t> oldData(b->data, b->data + oldDataLen);
    
    // Allocate new zeroed buffer
    std::unique_ptr<uint8_t, decltype(&std::free)> 
      newBuf(
	     static_cast<uint8_t*>(std::calloc(newDataLen, 1)),
	     &std::free
	     );
    
    // 1) copy qname
    std::memcpy(newBuf.get(), oldData.data(), qNameBytes);
    
    // 2) write new CIGAR
    uint32_t* cigDst = reinterpret_cast<uint32_t*>(newBuf.get() + qNameBytes);
    std::transform(
		   c.begin(), c.end(), cigDst,
		   [](const CigarField& cf){ return cf.raw(); }
		   );
    
    // 3) copy the rest (sequence+aux)
    std::memcpy(
		newBuf.get() + qNameBytes + newCigarBytes,
		oldData.data() + auxOffset,
		auxLen
		);
    
    // 4) swap into bam1_t
    b->core.n_cigar = static_cast<uint32_t>(newNCigar);
    b->l_data      = static_cast<int>(newDataLen);
    b->m_data      = b->l_data;
    
    // free old buffer and install the new one
    std::free(b->data);
    b->data = newBuf.release();
  }

  // void BamRecord::SmartAddTag(const std::string& tag, const std::string& val)
  // {
  //   // get the old tag
  //   assert(tag.length());
  //   assert(val.length());
  //   std::string tmp;
  //   GetZTag(tag, tmp);
  //   if (!tmp.length()) 
  //     {
  // 	AddZTag(tag, val);
  // 	return;
  //     }
    
  //   // check that we don't have the tag delimiter in the stirng
  //   if (val.find(TAG_DELIMITER) != std::string::npos)
  //     std::cerr << "BamRecord::SmartAddTag -- Tag delimiter " << TAG_DELIMITER << " is in the value to be added. Compile with diff tag delimiter or change val" << std::endl;

  //   // append the tag
  //   tmp.append(TAG_DELIMITER).append(val);
    
  //   // remove the old tag
  //   RemoveTag(tag.c_str());
    
  //   // add the new one
  //   assert(tmp.length());
  //   AddZTag(tag, tmp);
  // }

  void BamRecord::ClearSeqQualAndTags() {
    if (!b) return;
    
    // Compute the new buffer size: qname + all cigar ops (4 bytes each)
    const size_t newSize = 
      static_cast<size_t>(b->core.l_qname) 
      + static_cast<size_t>(b->core.n_cigar) * sizeof(uint32_t);
    
    // Shrink the data buffer
    uint8_t* newData = reinterpret_cast<uint8_t*>(
						  std::realloc(b->data, newSize)
						  );
    if (!newData && newSize > 0) {
      // realloc failed, original buffer is still valid but we can't proceed
      throw std::bad_alloc();
    }
    
    b->data   = newData;
    b->l_data = static_cast<int>(newSize);
    b->core.l_qseq = 0;
  }
  
  void BamRecord::SetSequence(std::string_view seq) {
    if (!b) return;
    
    // old sequence lengths and offsets
    size_t oldSeqBases = static_cast<size_t>(b->core.l_qseq);
    size_t oldSeqBytes = (oldSeqBases + 1) / 2 + oldSeqBases;
    size_t dataHeaderSize = static_cast<size_t>(b->core.l_qname) + b->core.n_cigar * sizeof(uint32_t);
    size_t oldAuxSpot = dataHeaderSize + oldSeqBytes;
    size_t oldAuxLen = static_cast<size_t>(bam_get_l_aux(b)) + oldSeqBytes;
    
    // new sequence lengths
    size_t newSeqBases = seq.size();
    size_t newSeqBytes = (newSeqBases + 1) / 2 + newSeqBases;
    
    // compute new overall buffer size
    size_t newSize = dataHeaderSize + newSeqBytes + (b->l_data - (dataHeaderSize + oldSeqBytes));
    
    // backup old data
    uint8_t* oldData = reinterpret_cast<uint8_t*>(std::malloc(b->l_data));
    std::memcpy(oldData, b->data, b->l_data);
    
    // reallocate data buffer
    std::free(b->data);
    b->data = reinterpret_cast<uint8_t*>(std::calloc(newSize, 1));
    
    // copy header (qname + CIGAR)
    std::memcpy(b->data, oldData, dataHeaderSize);
    
    // encode new sequence in 4-bit packed format
    uint8_t* seqPtr = b->data + dataHeaderSize;
    for (size_t i = 0; i < newSeqBases; ++i) {
      uint8_t code = 15;
      switch (seq[i]) {
      case 'A': code = 1; break;
      case 'C': code = 2; break;
      case 'G': code = 4; break;
      case 'T': code = 8; break;
      }
      size_t byteIdx = i / 2;
      size_t shift = ((~i & 1) << 2);
      seqPtr[byteIdx] &= static_cast<uint8_t>(~(0xF << shift));
      seqPtr[byteIdx] |= static_cast<uint8_t>(code << shift);
    }
    
    // update header fields
    b->core.l_qseq = static_cast<int>(newSeqBases);
    b->l_data = static_cast<int>(newSize);
    b->m_data = b->l_data;
    
    // null-quality marker
    if (uint8_t* qual = bam_get_qual(b)) qual[0] = 0xFF;
    
    // copy auxiliary tags
    if (uint8_t* aux = bam_get_aux(b))
      std::memcpy(aux, oldData + oldAuxSpot, oldAuxLen);
    
    std::free(oldData);
  }

  void BamRecord::SetQname(std::string_view name) {
    if (!b) return;
    
    // backup non-qname data
    size_t oldQnameLen = static_cast<size_t>(b->core.l_qname);
    size_t oldDataLen  = static_cast<size_t>(b->l_data);
    size_t nonQNameLen = oldDataLen - oldQnameLen;
    
    auto* nonQNameData = reinterpret_cast<uint8_t*>(std::malloc(nonQNameLen));
    std::memcpy(nonQNameData, b->data + oldQnameLen, nonQNameLen);
    
    // reallocate with new qname length
    std::free(b->data);
    size_t newQnameLen = name.size() + 1;  // include room for '\0'
    size_t newSize     = newQnameLen + nonQNameLen;
    b->data = reinterpret_cast<uint8_t*>(std::calloc(newSize, 1));
    
    // copy new qname (not NUL-terminated by data()), then add '\0'
    std::memcpy(b->data, name.data(), name.size());
    b->data[name.size()] = '\0';
    
    // copy the rest of the old data back in
    std::memcpy(b->data + newQnameLen, nonQNameData, nonQNameLen);
    
    // update header fields
    b->core.l_qname = static_cast<int>(newQnameLen);
    b->l_data       = static_cast<int>(newSize);
    b->m_data       = b->l_data;
    
    std::free(nonQNameData);
  }

  void BamRecord::SetQualities(std::string_view quals, int offset) {
    // If non empty, must match sequence length
    auto qlen = static_cast<size_t>(b->core.l_qseq);
    if (!quals.empty() && quals.size() != qlen) {
      throw std::invalid_argument("New quality string must match sequence length");
    }
    
    uint8_t* out = bam_get_qual(b);
    if (quals.empty()) {
      // signal no qualities in BAM by writing a 0-qual byte
      out[0] = 0;
      return;
    }
    
    // apply offset and write in place
    for (size_t i = 0; i < qlen; ++i) {
      out[i] = static_cast<uint8_t>(quals[i] - offset);
    }
  }
  
  // double BamRecord::MeanPhred() const {

  //   if (b->core.l_qseq <= 0)
  //     return -1;

  //   double s = 0;
  //   uint8_t* p = bam_get_qual(b);
  //   for (int32_t i = 0; i < b->core.l_qseq; ++i)
  //     s += p[i];
  //   return s / b->core.l_qseq;
  // }

  // std::string BamRecord::QualitySequence() const {
  //   std::string seq;
  //   GetZTag("GV", seq);
  //   if (!seq.length()) 
  //     seq = Sequence();
  //   return seq;
  // }
  
  std::ostream& operator<<(std::ostream& out, const BamRecord& r) {
    if (r.b == nullptr) {
      return out << "empty read\n";
    }
    
    // Cache a reference to the SAM core for brevity
    const auto& core = r.b->core;
    
    out
      << bam_get_qname(r.b)       << '\t'
      << core.flag                << '\t'
      << (core.tid   + 1)         << '\t'
      << core.pos                 << '\t'
      << static_cast<int>(core.qual) << '\t'
      << r.CigarString()          << '\t'
      << (core.mtid  + 1)         << '\t'
      << core.mpos                << '\t'
      << r.FullInsertSize()       << '\t'
      << r.Sequence()             << '\t'
      << '*'
      << '\n';  // single newline, no flush
    
    return out;
  }

  // int32_t BamRecord::CountBWASecondaryAlignments() const 
  // {
  //   int xp_count = 0;
    
  //   // xa tag
  //   std::string xar_s;
  //   GetZTag("XA", xar_s);
  //   if (xar_s.length()) {
  //     xp_count += std::count(xar_s.begin(), xar_s.end(), ';');
  //   }

  //   return xp_count;
    
  // }

  // int32_t BamRecord::CountBWAChimericAlignments() const 
  // {
  //   int xp_count = 0;
    
  //   // sa tag (post bwa mem v0.7.5)
  //   std::string xar_s;
  //   GetZTag("SA", xar_s);
  //   if (xar_s.length()) 
  //     xp_count += std::count(xar_s.begin(), xar_s.end(), ';');

  //   // xp tag (pre bwa mem v0.7.5)
  //   std::string xpr_s;
  //   GetZTag("XP", xpr_s);
  //   if (xpr_s.length()) 
  //     xp_count += std::count(xpr_s.begin(), xpr_s.end(), ';');

  //   return xp_count;
    
  // }

  int32_t BamRecord::CountNBases() const {
    const auto* seq_data = bam_get_seq(b);
    const size_t len = b->core.l_qseq;
    int32_t count = 0;
    for (size_t i = 0; i < len; ++i) {
      // bam_seqi pulls out the 4-bit code; 15 is the N sentinel
      count += (bam_seqi(seq_data, i) == 15);
    }
    return count;
  }

  void BamRecord::QualityTrimmedSequence(int32_t qualTrim,
					 int32_t& startpoint,
					 int32_t& endpoint) const 
  {
    const auto* qual = bam_get_qual(b.get());
    const int32_t len = b->core.l_qseq;
    
    // defaults
    startpoint = 0;
    endpoint   = -1;
    
    // no real qualities whole read
    if (len <= 0 || qual[0] == 0xff) 
      return;
    
    auto* begin = qual;
    auto* endp  = qual + len;
    
    // find first base qualTrim
    auto it1 = std::find_if(begin, endp,
                            [qualTrim](uint8_t qv){ return qv >= qualTrim; });
    if (it1 == endp) {
      // no base meets the threshold
      startpoint = len;
        return;
    }
    startpoint = static_cast<int32_t>(it1 - begin);
    
    // find last base qualTrim via reverse_iterator
    auto rit = std::find_if(std::make_reverse_iterator(endp),
                            std::make_reverse_iterator(begin),
                            [qualTrim](uint8_t qv){ return qv >= qualTrim; });
    endpoint = static_cast<int32_t>(std::distance(begin, rit.base()));
  }
  
  bool BamRecord::GetTag(std::string_view tag, std::string& out) const {
    // Try string (Z) tags first
    if (GetZTag(tag, out))
      return true;
    
    // Try integer tags
    if (int32_t ti; GetIntTag(tag, ti)) {
      out = std::to_string(ti);
      return true;
    }
    
    // Try float tags
    if (float tf; GetFloatTag(tag, tf)) {
      out = std::to_string(tf);
      return true;
    }
    
    return false;
  }
  
  void BamRecord::AddZTag(std::string_view tag, std::string_view val) {
    if (tag.empty() || val.empty()) return;
    // bam_aux_append expects a C-string with trailing '\0'
    bam_aux_append(b.get(),
		   tag.data(),
		   'Z',
		   static_cast<int>(val.size() + 1),
		   reinterpret_cast<const uint8_t*>(val.data()));
  }
  
  bool BamRecord::GetZTag(std::string_view tag, std::string& out) const {
    if (auto* p = bam_aux_get(b.get(), tag.data())) {
      if (*p == 'Z') {
	if (char* z = bam_aux2Z(p)) {
	  out.assign(z);
	  return true;
	}
      }
    }
    return false;
  } 
  
  // get a string tag that might be separted by "x"
  // std::vector<std::string> BamRecord::GetSmartStringTag(const std::string& tag) const {
    
  //   std::vector<std::string> out;
  //   std::string tmp;
  //   GetZTag(tag, tmp);

  //   if (tmp.empty())
  //     return std::vector<std::string>();
    
  //   if (tmp.find(TAG_DELIMITER) != std::string::npos) {
  //     std::istringstream iss(tmp);
  //     std::string line;
  //     while (std::getline(iss, line, CTAG_DELIMITER)) {
  // 	out.push_back(line);
  //     }
  //   } else {
  //     out.push_back(tmp);
  //   }
    
  //   assert(out.size());
  //   return out;
    
  // }
  
  
  // std::vector<int> BamRecord::GetSmartIntTag(const std::string& tag) const {
    
  //   std::vector<int> out;
  //   std::string tmp;
    
  //   GetZTag(tag, tmp);
  //   if (tmp.empty())
  //     return std::vector<int>();
    
  //   if (tmp.find(TAG_DELIMITER) != std::string::npos) {
  //     std::istringstream iss(tmp);
  //     std::string line;
  //     while (std::getline(iss, line, CTAG_DELIMITER))
  // 	out.push_back(atoi(line.c_str())); 
  //   } else {
  //     out.push_back(atoi(tmp.c_str())); 
  //   }
    
  //   assert(out.size());
  //   return out;
    
  // }

  // std::vector<double> BamRecord::GetSmartDoubleTag(const std::string& tag) const {
    
  //   std::vector<double> out;
  //   std::string tmp;
    
  //   GetZTag(tag, tmp);
  //   if (tmp.empty())
  //     return std::vector<double>();
    
  //   if (tmp.find(TAG_DELIMITER) != std::string::npos) {
  //     std::istringstream iss(tmp);
  //     std::string line;
  //     while (std::getline(iss, line, CTAG_DELIMITER))
  // 	out.push_back(std::atof(line.c_str())); 
  //   } else { // single entry
  //     out.push_back(std::atof(tmp.c_str())); 
  //   }
    
  //   assert(out.size());
  //   return out;
    
  // }
 
  bool BamRecord::operator<(const BamRecord& other) const {
    return std::tie(b->core.tid, b->core.pos)
      < std::tie(other.b->core.tid, other.b->core.pos);
  }
  
  bool BamRecord::operator==(const BamRecord& other) const {
    return std::tie(b->core.tid, b->core.pos)
      == std::tie(other.b->core.tid, other.b->core.pos);
  }
  
  // Cigar::Cigar(const std::string& cig) {

  //   // get the ops (MIDSHPN)
  //   std::vector<char> ops;
  //   for (size_t i = 0; i < cig.length(); ++i)
  //     if (!isdigit(cig.at(i))) {
  // 	ops.push_back(cig.at(i));
  //     }
    
  //   std::size_t prev = 0, pos;
  //   std::vector<std::string> lens;
  //   while ((pos = cig.find_first_of("MIDSHPNX", prev)) != std::string::npos) {
  //       if (pos > prev)
  // 	  lens.push_back(cig.substr(prev, pos-prev));
  //       prev = pos+1;
  //     }
  //   if (prev < cig.length())
  //     lens.push_back(cig.substr(prev, std::string::npos));

  //   assert(ops.size() == lens.size());
  //   for (size_t i = 0; i < lens.size(); ++i) {
  //     add(CigarField(ops[i], std::atoi(lens[i].c_str())));
  //   }
    
  //   //return tc;

  // }

  // bool Cigar::operator==(const Cigar& c) const { 
  //    if (m_data.size() != c.size())
  //      return false;
  //    if (!m_data.size()) // both empty
  //      return true;
  //    for (size_t i = 0; i < m_data.size(); ++i)
  //      if (m_data[i].Type() != c[i].Type() || m_data[i].Length() != c[i].Length())
  // 	 return false;
  //    return true;
  // }

  int BamRecord::OverlappingCoverage(const BamRecord& r) const {
    // 1) Build a coverage mask for *this* read
    const auto len1 = GetCigar().NumQueryConsumed();
    std::vector<uint8_t> cov1(len1, 0);
    
    size_t pos = 0;
    for (const auto& f : GetCigar()) {
      if (f.Type() == 'M') {
	for (size_t j = 0; j < f.Length(); ++j) {
	  cov1[pos + j] = 1;
	}
      }
      if (f.ConsumesQuery()) {
	pos += f.Length();
      }
    }
    
    // 2) Now count overlaps in r
    size_t ocov = 0;
    pos = 0;
    for (const auto& f : r.GetCigar()) {
      if (f.Type() == 'M') {
	for (size_t j = 0; j < f.Length(); ++j) {
	  if (cov1[pos + j]) {
	    ++ocov;
	  }
	}
      }
      if (f.ConsumesQuery()) {
	pos += f.Length();
      }
    }
    
    return static_cast<int>(ocov);
  }
  
  bool BamRecord::GetFloatTag(std::string_view tag, float& out) const {
    if (!b) return false;                           // guard if record is empty
    // bam_aux_get expects a null-terminated C-string; string_view::data()
    // is fine if the view comes from a null-terminated string literal
    // or std::string, otherwise you may need to materialize it.
    uint8_t* aux = bam_aux_get(b.get(), tag.data());
    if (!aux) return false;                          // tag not present
    
    char type = *aux++;
    if (type != 'f' && type != 'd')                  // not a float/double tag
      return false;
    
    out = bam_aux2f(aux);                            // extract as float
    return true;
  }

  bool BamRecord::GetIntTag(std::string_view tag, int32_t& out) const {
    if (!b) return false;                                 // no record
    uint8_t* aux = bam_aux_get(b.get(), tag.data());
    if (!aux) return false;                                // tag not present
    
    char type = *aux++;
    // only integer-like tags:
    if (type!='i' && type!='I' && type!='c' && type!='C' &&
        type!='s' && type!='S')
      return false;
    
    out = bam_aux2i(aux);
    return true;
  }

  std::string BamRecord::CigarString() const {
    if (!b) return {};  // or throw if you prefer

    const uint32_t n = b->core.n_cigar;
    auto* cig = bam_get_cigar(b.get());
    static constexpr char ops[] = "MIDNSHP=XB";

    // pre-reserve a bit (avg ~34 chars per op)
    std::string out;
    out.reserve(n * 4);

    for (uint32_t i = 0; i < n; ++i) {
        uint32_t len = bam_cigar_oplen(cig[i]);
        char op   = ops[cig[i] & BAM_CIGAR_MASK];
        out += std::to_string(len);
        out.push_back(op);
    }
    return out;
  }

  // -----------------------------------------------------------------------------
  // --- getters ---------------------------------------------------------------
  
  int32_t BamRecord::ChrID() const {
    if (!b)
      return -1;
    return b->core.tid;
  }
  
  int32_t BamRecord::MateChrID() const {
    return b ? b->core.mtid : -1;
  }
  
  int32_t BamRecord::MapQuality() const {
    return b ? b->core.qual : -1;
  }
  
  // -----------------------------------------------------------------------------
  // --- setters / flags --------------------------------------------------------

  void BamRecord::SetMapQuality(int32_t m) {
    if (b)
      b->core.qual = m;
  }
  
  void BamRecord::SetChrID(int32_t i) {
    if (b)
      b->core.tid = i;
  }
  
  void BamRecord::SetChrIDMate(int32_t i) {
    if (b)
      b->core.mtid = i;
  }
  
  void BamRecord::SetPositionMate(int32_t i) {
    if (b)
      b->core.mpos = i;
  }

  void BamRecord::SetQCFail(bool f) {
    if (!b) return;
    if (f)  b->core.flag |= BAM_FQCFAIL;
    else    b->core.flag &= ~BAM_FQCFAIL;
  }
  
  void BamRecord::SetPairMappedFlag(bool f) {
    if (!b) return;
    if (f)  b->core.flag |= BAM_FPAIRED;
    else    b->core.flag &= ~BAM_FPAIRED;
  }
  
  void BamRecord::SetMateReverseFlag(bool f) {
    if (!b) return;
    if (f)  b->core.flag |= BAM_FMREVERSE;
    else    b->core.flag &= ~BAM_FMREVERSE;
  }

  void BamRecord::AddIntTag(std::string_view tag, int32_t val) {
    if (!b) return;
    // bam_aux_append wants: (bam1_t*, const char* tag, char type, int len, const uint8_t* data)
    bam_aux_append(
		   b.get(),
		   tag.data(),         // no longer require a null-terminated std::string
		   'i',
		   sizeof(val),
		   reinterpret_cast<const uint8_t*>(&val)
		   );
  }
  
  void BamRecord::SetID(int32_t id) {
    if (b) b->core.tid = id;
  }
  
  void BamRecord::SetPosition(int32_t pos) {
    if (b) b->core.pos = pos;
  }
  
  // --------------
  // --------------
  
  std::string BamRecord::ParseReadGroup() const {
    // 1) look for explicit RG tag
    std::string rg;
    if (GetZTag("RG", rg))
      return rg;
    
    // 2) else fall back to the prefix of the QNAME
    auto qn = Qname();
    if (auto pos = qn.find(':'); pos != std::string::npos)
      return qn.substr(0, pos);
    
    // 3) neither found  NA
    return "NA";
  }

  int BamRecord::NumAlignedBases() const {
  int total = 0;
  for (const auto& f : GetCigar()) {
    switch (f.Type()) {
      case 'M': case 'I': case '=': case 'X': case 'D':
        total += f.Length();
        break;
      default:
        break;
    }
  }
  return total;
  }
  
  uint32_t BamRecord::MaxInsertionBases() const {
    uint32_t best = 0;
    for (const auto& f : GetCigar()) {
      if (f.Type() == 'I')
	best = std::max(best, f.Length());
    }
    return best;
  }
  
  uint32_t BamRecord::MaxDeletionBases() const {
    uint32_t best = 0;
    for (const auto& f : GetCigar()) {
      if (f.Type() == 'D')
	best = std::max(best, f.Length());
    }
    return best;
  }
  
  uint32_t BamRecord::NumMatchBases() const {
    uint32_t total = 0;
    for (const auto& f : GetCigar()) {
      if (f.Type() == 'M')
	total += f.Length();
    }
    return total;
  }

  Cigar BamRecord::GetCigar() const {
    Cigar cig;
    const auto n = b->core.n_cigar;
    cig.reserve(n);                         // avoid reallocations if supported
    auto* raw = bam_get_cigar(b.get());     // aligned, safe pointer
    
    for (size_t i = 0; i < n; ++i) {
      cig.add(CigarField{raw[i]});
    }
    return cig;
  }
  
  Cigar BamRecord::GetReverseCigar() const {
    Cigar cig;
    const auto n = b->core.n_cigar;
    cig.reserve(n);
    auto* raw = bam_get_cigar(b.get());
    
    // reverse-iterate from raw[n-1] down to raw[0]
    for (size_t i = n; i-- > 0; ) {
      cig.add(CigarField{raw[i]});
    }
    return cig;
  }
  
  std::string BamRecord::Qualities(int offset) const { 
    if (!b) return {};
    auto* qual = bam_get_qual(b.get());
    if (!qual) return {};
    
    std::string out;
    out.reserve(b->core.l_qseq);
    for (int i = 0; i < b->core.l_qseq; ++i) 
      out.push_back(static_cast<char>(qual[i] + offset));
    return out;
  }
  
  int32_t BamRecord::AlignmentPositionReverse() const {
    if (!b) return -1;
    auto* cig = bam_get_cigar(b.get());
    const size_t n = b->core.n_cigar;
    int32_t pos = 0;
    for (size_t i = n; i-- > 0; ) {
      char op = bam_cigar_opchr(cig[i]);
      if (op == 'S' || op == 'H') 
	pos += bam_cigar_oplen(cig[i]);
      else 
	break;
    }
    return pos;
  }
  
  int32_t BamRecord::AlignmentEndPositionReverse() const {
    if (!b) return -1;
    auto* cig = bam_get_cigar(b.get());
    const size_t n = b->core.n_cigar;
    int32_t clip = 0;
    // count leading S/H in forward orientation
    for (size_t i = 0; i < n; ++i) {
      char op = bam_cigar_opchr(cig[i]);
      if (op == 'S' || op == 'H')
	clip += bam_cigar_oplen(cig[i]);
      else
	break;
    }
    return b->core.l_qseq - clip;
  }
  
  int32_t BamRecord::AlignmentPosition() const {
    if (!b) return -1;
    auto* cig = bam_get_cigar(b.get());
    const size_t n = b->core.n_cigar;
    int32_t pos = 0;
    // skip initial H then S
    for (size_t i = 0; i < n; ++i) {
      char op = bam_cigar_opchr(cig[i]);
      if (op == 'H') 
	continue;
      if (op == 'S') 
	pos += bam_cigar_oplen(cig[i]);
      else 
	break;
    }
    return pos;
  }
  
  int32_t BamRecord::AlignmentEndPosition() const {
    if (!b) return -1;
    auto* cig = bam_get_cigar(b.get());
    const size_t n = b->core.n_cigar;
    int32_t clip = 0;
    // skip trailing S/H in reverse orientation
    for (size_t i = n; i-- > 0; ) {
      char op = bam_cigar_opchr(cig[i]);
      if (op == 'S' || op == 'H')
	clip += bam_cigar_oplen(cig[i]);
      else
	break;
    }
    return b->core.l_qseq - clip;
  }

  int32_t BamRecord::NumSoftClip() const {
    if (!b) return 0;
    auto* cig = bam_get_cigar(b.get());
    int32_t total = 0;
    uint32_t n = b->core.n_cigar;
    for (uint32_t i = 0; i < n; ++i) {
      if (bam_cigar_opchr(cig[i]) == 'S') 
	total += bam_cigar_oplen(cig[i]);
    }
    return total;
  }
  
  int32_t BamRecord::NumHardClip() const {
    if (!b) return 0;
    auto* cig = bam_get_cigar(b.get());
    int32_t total = 0;
    uint32_t n = b->core.n_cigar;
    for (uint32_t i = 0; i < n; ++i) {
      if (bam_cigar_opchr(cig[i]) == 'H') 
	total += bam_cigar_oplen(cig[i]);
    }
    return total;
  }
  
  int32_t BamRecord::NumClip() const {
    if (!b) return 0;
    auto* cig = bam_get_cigar(b.get());
    int32_t total = 0;
    uint32_t n = b->core.n_cigar;
    for (uint32_t i = 0; i < n; ++i) {
      char op = bam_cigar_opchr(cig[i]);
      if (op == 'S' || op == 'H') 
	total += bam_cigar_oplen(cig[i]);
    }
    return total;
  }
  
  std::string BamRecord::ChrName(const BamHeader& hdr) const {
    if (!b || b->core.tid < 0) return {};
    if (!hdr.isEmpty())            // named contigs
      return hdr.IDtoName(b->core.tid);
    // fallback to stringified index
    return std::to_string(b->core.tid);
  }
  
  std::string BamRecord::Brief() const {
    if (!b) return {};
    // 1-based chrom, comma-formatted pos, strand
    char strand = (b->core.flag & BAM_FREVERSE) ? '-' : '+';
    return std::to_string(b->core.tid + 1)
      + ':' + AddCommas(b->core.pos)
      + '(' + strand + ')';
  }
  
  std::string BamRecord::BriefMate() const {
    if (!b) return {};
    char strand = (b->core.flag & BAM_FMREVERSE) ? '-' : '+';
    return std::to_string(b->core.mtid + 1)
      + ':' + AddCommas(b->core.mpos)
      + '(' + strand + ')';
  }

  int BamRecord::PairOrientation() const {
    if (!PairMappedFlag())
      return UDORIENTATION;
    
    const bool rev      = ReverseFlag();
    const bool mate_rev = MateReverseFlag();
    const auto pos      = Position();
    const auto mpos     = MatePosition();
    
    // FR: first forward, mate reverse
    if ((!rev && pos <= mpos &&  mate_rev) ||
	( rev && pos >= mpos && !mate_rev))
      return FRORIENTATION;
    
    // FF: both forward
    if (!rev && !mate_rev)
      return FFORIENTATION;
    
    // RR: both reverse
    if (rev && mate_rev)
      return RRORIENTATION;
    
    // RF: read reverse, mate forward
    if (( rev && pos <  mpos && !mate_rev) ||
	(!rev && pos >  mpos &&  mate_rev))
      return RFORIENTATION;
    
    // shouldn't happen, but fall back to unoriented
    return UDORIENTATION;
  }
  
  bool BamRecord::ProperOrientation() const {
    if (!b) 
      return false;
    
    // must be on the same chromosome
    if (b->core.tid != b->core.mtid) 
      return false;
    
    const bool rev      = ReverseFlag();
    const bool mate_rev = MateReverseFlag();
    const auto pos      = Position();
    const auto mpos     = MatePosition();
    
    if (pos < mpos) {
      // "FR" orientation: this read forward, mate reverse
      return (!rev && mate_rev);
    } else {
      // "RF" orientation: this read reverse, mate forward
      return (rev && !mate_rev);
    }
  }

  
}
