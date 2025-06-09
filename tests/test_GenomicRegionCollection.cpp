#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>
#include "SeqLib/GenomicRegionCollection.h"

TEST_CASE("Shuffle randomizes but preserves size", "[GenomicRegionCollection]") {
  SeqLib::GenomicRegionCollection grc;
  grc.add(SeqLib::GenomicRegion(0, 100, 200, '+'));
  grc.add(SeqLib::GenomicRegion(1, 300, 400, '-'));
  size_t before = grc.size();

  grc.Shuffle();
  REQUIRE(grc.size() == before);

}
