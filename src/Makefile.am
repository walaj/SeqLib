SUBDIRS = bwa

noinst_LIBRARIES=libsnowtools.a

libsnowtools_a_SOURCES = GenomicRegion.cpp Histogram.cpp \
			BamWalker.cpp BamStats.cpp MiniRules.cpp gzstream.C SnowUtils.cpp \
			gChain.cpp HTSTools.cpp BamRead.cpp gzstream.C BWAWrapper.cpp \
			bwa/is.c bwa/bwtindex.c \
			bwa/bwt_gen.c bwa/QSufSort.c bwa/bwt.c bwa/bwamem_extra.c

bin_PROGRAMS = snowtools

snowtools_LDADD = $(top_builddir)/libsnowtools.a $(top_builddir)/bwa/libbwa.a
snowtools_LDFLAGS = -pthread -std=c++11
snowtools_SOURCES = SnowTools.cpp CoverageTool.cpp RealignTool.cpp
