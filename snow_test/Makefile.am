bin_PROGRAMS = snow_test

#noinst_LIBRARIES = libsnowtest.a

snow_test_CPPFLAGS = \
     -I../src \
     -I../htslib \
     -I.. \
     -I$(top_srcdir)/../multifast-v1.4.2 --coverage

snow_test_LDADD = \
	$(top_builddir)/../src/libsnowtools.a \
	$(top_builddir)/../bwa/libbwa.a \
	$(top_builddir)/../htslib/libhts.a \
	@boost_lib@/libboost_regex.a @boost_lib@/libboost_unit_test_framework.a \
	@boost_lib@/libboost_system.a #libsnowtest

snow_test_LDFLAGS = --coverage ##-BOOST_TEST_DYN_LINK

snow_test_SOURCES = snow_test.cpp snow_test_main.cpp \
	$(top_builddir)/../src/GenomicRegion.cpp \
	$(top_builddir)/../src/Histogram.cpp $(top_builddir)/../src/BamWalker.cpp $(top_builddir)/../src/BamStats.cpp \
	$(top_builddir)/../src/MiniRules.cpp $(top_builddir)/../src/SnowUtils.cpp $(top_builddir)/../src/BamRead.cpp \
	$(top_builddir)/../src/BWAWrapper.cpp $(top_builddir)/../src/STCoverage.cpp \
        $(top_builddir)/../src/Fractions.cpp
