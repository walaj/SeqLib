#include <catch2/catch.hpp>
#include "SeqLib/BamRecord.h"
#include <sstream>

using SeqLib::Cigar;
using SeqLib::CigarField;

TEST_CASE("CigarField prints correct string", "[CigarField]") {
    std::stringstream ss;

    CigarField f1('M', 10);
    ss << f1;
    REQUIRE(ss.str() == "10M");

    ss.str(""); ss.clear();
    CigarField f2('I', 3);
    ss << f2;
    REQUIRE(ss.str() == "3I");

    ss.str(""); ss.clear();
    CigarField f3('D', 7);
    ss << f3;
    REQUIRE(ss.str() == "7D");
}

TEST_CASE("Cigar parses a CIGAR string correctly", "[Cigar]") {
    Cigar c("5M2I3D4S");
    std::vector<std::pair<char,int>> expected = {
        {'M', 5},
        {'I', 2},
        {'D', 3},
        {'S', 4}
    };

    size_t idx = 0;
    for (const auto& cf : c) {
        std::stringstream ss;
        ss << cf;
        REQUIRE(ss.str() == std::to_string(expected[idx].second) + expected[idx].first);
        ++idx;
    }
    REQUIRE(idx == expected.size());
}

TEST_CASE("Cigar equality operator", "[Cigar]") {
    Cigar a("10M5I");
    Cigar b;
    b.add(CigarField('M', 10));
    b.add(CigarField('I', 5));
    REQUIRE(a == b);

    Cigar c("10M4I");
    REQUIRE_FALSE(a == c);
}

TEST_CASE("NumQueryConsumed sums M, I, and S correctly", "[Cigar]") {
    // M=5, I=2, S=4 => total query consumed = 5+2+4 = 11
    Cigar c("5M2I3D4S");
    REQUIRE(c.NumQueryConsumed() == 11);
}

TEST_CASE("NumReferenceConsumed sums M and D correctly", "[Cigar]") {
    // M=5, D=3 => total reference consumed = 5+3 = 8
    Cigar c("5M2I3D4S");
    REQUIRE(c.NumReferenceConsumed() == 8);
}
