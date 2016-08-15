#include "SeqLib/BamHeader.h"
#include "SeqLib/BamWalker.h"

#include <sstream>
#include <stdexcept>
#include <iostream>

#include "htslib/khash.h"

namespace SeqLib {

BamHeader::BamHeader(const std::string& hdr)  {

  h = std::shared_ptr<bam_hdr_t>(sam_hdr_read2(hdr)); 

  ConstructName2IDTable();

}

  void BamHeader::WriteToStdout() const {

    std::shared_ptr<htsFile> f = std::shared_ptr<htsFile>(sam_open("-", "w"), sam_write_delete());
    f->format.format == text_format;
    sam_hdr_write(f.get(), h.get());
  }

  BamHeader::BamHeader(const bam_hdr_t * hdr) {

    h = std::shared_ptr<bam_hdr_t>(bam_hdr_dup(hdr));

    ConstructName2IDTable();

  }

  int BamHeader::ConstructName2IDTable() {

    // create the lookup table if not already made
    if (!n2i) {
      n2i = std::shared_ptr<std::unordered_map<std::string, int>>(new std::unordered_map<std::string, int>());
      for (int i = 0; i < h->n_targets; ++i)
	n2i->insert(std::pair<std::string, int>(std::string(h->target_name[i]), i));
    }

  }

  int BamHeader::Name2ID(const std::string& name) const {

    std::unordered_map<std::string, int>::const_iterator ff = n2i->find(name);
    if (ff != n2i->end())
      return ff->second;
    else
      return -1;

  }

bam_hdr_t* BamHeader::sam_hdr_read2(const std::string& hdr) const {

  kstring_t str;
  bam_hdr_t *hhh;
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
  hhh = sam_hdr_parse(str.l, str.s);
  hhh->l_text = str.l; hhh->text = str.s;
  return hhh;
}

int BamHeader::NumSequences() const {
  
  if (!h)
    return 0;
  return h->n_targets;

}

std::string BamHeader::IDtoName(int id) const {

  if (id < 0)
    throw std::invalid_argument("BamHeader::IDtoName - ID must be >= 0");

  if (!h)
    throw std::out_of_range("BamHeader::IDtoName - Header is uninitialized");

  if (id >= h->n_targets)
    throw std::out_of_range("BamHeader::IDtoName - Requested ID is higher than number of sequences");
  
  return std::string(h->target_name[id]);

}

}
