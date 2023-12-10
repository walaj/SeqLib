AUTOMAKE_OPTIONS = foreign
SUBDIRS = bwa fermi-lite src

install:  
	mkdir -p lib && cp src/libseqlib.a fermi-lite/libfml.a bwa/libbwa.a lib

seqtools:
	mkdir -p bin && cd src/seqtools && make && mv seqtools ../../bin
