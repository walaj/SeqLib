bin_PROGRAMS = seq_test read_write 

seq_test_CPPFLAGS = \
     -I../src \
     -I../htslib \
     -I.. \
     -I../fermi-lite --coverage

seq_test_LDADD = \
	../fermi-lite/libfml.a \
	../bwa/libbwa.a \
	../htslib/libhts.a \
	-lboost_unit_test_framework -lboost_system -lboost_timer -lboost_chrono

seq_test_LDFLAGS = --coverage

seq_test_SOURCES = seq_test.cpp \
	../src/GenomicRegion.cpp \
	../src/BamWriter.cpp ../src/BamReader.cpp \
	../src/ReadFilter.cpp ../src/SeqLibUtils.cpp ../src/BamRecord.cpp \
	../src/BWAWrapper.cpp \
        ../src/RefGenome.cpp ../src/SeqPlot.cpp ../src/BamHeader.cpp \
	../src/FermiAssembler.cpp ../src/ssw_cpp.cpp ../src/ssw.c ../src/jsoncpp.cpp

read_write_CPPFLAGS = \
     -I../src \
     -I../htslib \
     -I.. \
     -I../fermi-lite --coverage

read_write_LDADD = \
	../fermi-lite/libfml.a \
	../bwa/libbwa.a \
	../htslib/libhts.a \
	-lboost_unit_test_framework -lboost_system -lboost_timer -lboost_chrono

read_write_LDFLAGS = --coverage

read_write_SOURCES = read_write.cpp \
	../src/GenomicRegion.cpp \
	../src/BamWriter.cpp ../src/BamReader.cpp \
	../src/ReadFilter.cpp ../src/SeqLibUtils.cpp ../src/BamRecord.cpp \
	../src/BWAWrapper.cpp ../src/BamPolyReader.cpp \
        ../src/RefGenome.cpp ../src/SeqPlot.cpp ../src/BamHeader.cpp \
	../src/FermiAssembler.cpp ../src/ssw_cpp.cpp ../src/ssw.c ../src/jsoncpp.cpp 
