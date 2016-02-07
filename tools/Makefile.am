bin_PROGRAMS = snow

snow_CPPFLAGS = \
        -I.. \
	-I../multifast-v1.4.2 \
        -I../htslib \
	-I../SGA/Util -I../SGA/Bigraph -I../SGA/SuffixTools \
	-I../SGA/Algorithm -I../SGA/SQG -I../blat/inc -I../src \
	-I../blat/jkOwnLib -Wno-write-strings

snow_LDADD = ../src/libsnowtools.a ../bwa/libbwa.a ../htslib/libhts.a ../blat/lib/jkOwnLib.a ../blat/lib/jkweb.a

snow_SOURCES = snowtools.cpp 

## these are the executable versions of snowtools
#snowtools_LDADD = $(top_builddir)/libsnowtools.a 
#snowtools_LDFLAGS = -pthread -std=c++11
#snowtools_SOURCES = SnowTools.cpp CoverageTool.cpp RealignTool.cpp
