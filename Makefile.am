AUTOMAKE_OPTIONS = foreign
SUBDIRS = bwa htslib fermi-lite src

install:  
	mkdir -p lib && cp src/libseqlib.a fermi-lite/libfml.a bwa/libbwa.a htslib/libhts.a lib

seqtools:
	mkdir -p bin && cd src/seqtools && make && mv seqtools ../../bin
