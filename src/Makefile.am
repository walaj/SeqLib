noinst_LIBRARIES = libseqlib.a

libseqlib_a_CPPFLAGS = -I../ -I../htslib

libseqlib_a_SOURCES =   FastqReader.cpp BFC.cpp ReadFilter.cpp SeqPlot.cpp jsoncpp.cpp ssw_cpp.cpp ssw.c \
			GenomicRegion.cpp RefGenome.cpp BamWriter.cpp BamReader.cpp \
			BWAWrapper.cpp BamRecord.cpp FermiAssembler.cpp BamHeader.cpp

bin_PROGRAMS = seqtools
seqtools_SOURCES = seqtools.cpp
seqtools_LDADD=libseqlib.a ../htslib/libhts.a ../bwa/libbwa.a ../fermi-lite/libfml.a
seqtools_CPPFLAGS=-I../htslib -I..

#seqtools.o: seqtools.cpp
#	g++ -c seqtools.cpp $(INCLUDES) $(CFLAGS) $(CXXFLAGS)

#seqtools: seqtools.o libseqlib.a
#	g++ seqtools.o -o seqtools $(LIBS) -lrt -lpthread -lz -lm
