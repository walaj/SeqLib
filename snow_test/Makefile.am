bin_PROGRAMS = snow_test

#noinst_LIBRARIES = libsnowtest.a

snow_test_CPPFLAGS = \
     -I../src \
     -I../htslib \
     -I.. \
     -I../blat/jkOwnLib -I../blat/inc \
     -I../multifast-v1.4.2 --coverage -Wno-write-strings

snow_test_LDADD = \
	../src/libsnowtools.a \
	../bwa/libbwa.a \
	../multifast-v1.4.2/ahocorasick/build/libahocorasick.a \
	../htslib/libhts.a ../blat/lib/jkOwnLib.a ../blat/lib/jkweb.a \
	@boost_lib@/libboost_regex.a @boost_lib@/libboost_unit_test_framework.a \
	@boost_lib@/libboost_system.a #libsnowtest

snow_test_LDFLAGS = --coverage ##-BOOST_TEST_DYN_LINK

snow_test_SOURCES = snow_test.cpp snow_test_main.cpp \
	../src/GenomicRegion.cpp \
	../src/Histogram.cpp ../src/BamWalker.cpp ../src/BamStats.cpp \
	../src/MiniRules.cpp ../src/SnowUtils.cpp ../src/BamRead.cpp \
	../src/BWAWrapper.cpp ../src/STCoverage.cpp \
        ../src/Fractions.cpp ../src/AlignedContig.cpp \
	../src/BreakPoint.cpp ../src/RefGenome.cpp \
	../src/BLATWrapper.cpp ../src/SVVCF.cpp
