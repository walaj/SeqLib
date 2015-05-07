#include "SnowTools/BamRead.h"

namespace SnowTools {

  BamRead::BamRead() {
    b = bam_init1();
  }

  BamRead::~BamRead() {
    bam_destroy1(b);
  }

}
