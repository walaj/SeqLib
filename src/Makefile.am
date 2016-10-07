noinst_LIBRARIES = libseqlib.a

libseqlib_a_CPPFLAGS = -I../ -I../htslib -I../blat/inc -I../blat/jkOwnLib -Wwrite-strings

libseqlib_a_SOURCES =   FastqReader.cpp BFC.cpp ReadFilter.cpp SeqPlot.cpp jsoncpp.cpp ssw_cpp.cpp ssw.c \
			GenomicRegion.cpp RefGenome.cpp BamWriter.cpp BamReader.cpp \
			BLATWrapper.cpp BWAWrapper.cpp BamRecord.cpp FermiAssembler.cpp BamHeader.cpp

LIBS=libseqlib.a ../htslib/libhts.a ../bwa/libbwa.a ../fermi-lite/libfml.a
INCLUDES=-I../htslib -I.. 

seqtools.o: seqtools.cpp
	g++ -c seqtools.cpp $(INCLUDES) $(CFLAGS) $(CXXFLAGS)

seqtools: seqtools.o libseqlib.a
	g++ seqtools.o -o seqtools $(LIBS) -lrt -lpthread -lz -lm
