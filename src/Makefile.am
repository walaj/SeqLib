noinst_LIBRARIES = libseqlib.a

libseqlib_a_CPPFLAGS = \
        -I.. -I../htslib -I../bwa -I../fermi-lite -I../blat/inc -I../blat/jkOwnLib -Wno-write-strings

libseqlib_a_SOURCES = ReadFilter.cpp SeqPlot.cpp jsoncpp.cpp ssw_cpp.cpp ssw.c \
			GenomicRegion.cpp RefGenome.cpp BamWriter.cpp BamReader.cpp \
			BLATWrapper.cpp BWAWrapper.cpp BamRecord.cpp FermiAssembler.cpp BamHeader.cpp
