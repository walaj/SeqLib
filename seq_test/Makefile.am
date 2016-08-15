bin_PROGRAMS = seq_test

#noinst_LIBRARIES = libseqtest.a

seq_test_CPPFLAGS = \
     -I../src \
     -I../htslib \
     -I../bwa \
     -I../fermi-lite \
     -I.. --coverage

seq_test_LDADD = \
	../src/libseqkit.a \
	../bwa/libbwa.a \
	../htslib/libhts.a \
	-lboost_unit_test_framework -lboost_system -lboost_timer -lboost_chrono

seq_test_LDFLAGS = --coverage

seq_test_SOURCES = seq_test.cpp \
	../src/GenomicRegion.cpp \
	../src/Histogram.cpp ../src/BamWriter.cpp ../src/BamReader.cpp ../src/BamStats.cpp \
	../src/ReadFilter.cpp ../src/SeqKitUtils.cpp ../src/BamRecord.cpp \
	../src/BWAWrapper.cpp ../src/STCoverage.cpp \
        ../src/Fractions.cpp ../src/RefGenome.cpp ../src/SeqPlot.cpp ../src/BamHeader.cpp \
	../src/FermiAssembler.cpp
