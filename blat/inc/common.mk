CC=gcc
COPT=-O -g
HG_DEFS=-D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -D_GNU_SOURCE -DMACHTYPE_${MACHTYPE}
HG_INC=-I../inc

#global external libraries 
L=

# pthreads is required
L+=-pthread

ifeq (${BINDIR},)
    BINDIR = ${HOME}/bin/${MACHTYPE}
endif

MKDIR=mkdir -p
ifeq (${STRIP},)
   STRIP=strip
endif

# portable naming of compiled executables: add ".exe" if compiled on 
# Windows (with cygwin).
ifeq (${OS}, Windows_NT)
  AOUT=a
  EXE=.exe
else
  AOUT=a.out
  EXE=
endif

%.o: %.c
	${CC} ${COPT} ${CFLAGS} ${HG_DEFS} ${LOWELAB_DEFS} ${HG_WARN} ${HG_INC} ${XINC} -o $@ -c $<

