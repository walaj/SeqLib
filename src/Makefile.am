noinst_LIBRARIES = libseqlib.a

<<<<<<< HEAD
libseqlib_a_CPPFLAGS = \
        -I.. -I../htslib -I../bwa -I../fermi-lite -I../blat/inc -I../blat/jkOwnLib -Wno-write-strings
=======
libseqlib_a_CPPFLAGS = -I../ -I..//htslib
>>>>>>> master

libseqlib_a_SOURCES =   FastqReader.cpp BFC.cpp ReadFilter.cpp SeqPlot.cpp jsoncpp.cpp ssw_cpp.cpp ssw.c \
			GenomicRegion.cpp RefGenome.cpp BamWriter.cpp BamReader.cpp \
<<<<<<< HEAD
			BLATWrapper.cpp BWAWrapper.cpp BamRecord.cpp FermiAssembler.cpp BamHeader.cpp
=======
			BWAWrapper.cpp BamRecord.cpp FermiAssembler.cpp BamHeader.cpp

bin_PROGRAMS = seqtools

seqtools_CPPFLAGS = -I.. -I../htslib

seqtools_SOURCES = seqtools.cpp

seqtools_LDADD = libseqlib.a ../bwa/libbwa.a ../fermi-lite/libfml.a ../htslib/libhts.a
>>>>>>> master
