#!/usr/bin/env bash

#set -x
set -e
set -o pipefail

echo "CXX: $CXX TRAVIS_BRANCH $TRAVIS_BRANCH CC $CC TRAVIS_OS_NAME $TRAVIS_OS_NAME"

## only build for one compiler
if [ "$COMPILER" == "g++-4.9" ] && [ "$TRAVIS_BRANCH" == "master" ]; 
then

    if [ "${TRAVIS_OS_NAME}" = "osx" ]; then
	brew update
    fi
    
    if [ -n "${BOOST_VERSION}" ]; then
	mkdir -p $BOOST_ROOT
	wget --no-verbose --output-document=- \
	    http://sourceforge.net/projects/boost/files/boost/${BOOST_VERSION}/boost_${BOOST_VERSION//./_}.tar.bz2/download \
	    | tar jxf - --strip-components=1 -C "${BOOST_ROOT}"
    fi

    ## install htslib
    mkdir htslib
    wget --no-verbose --output-document=- https://github.com/samtools/htslib/releases/download/1.12/htslib-1.12.tar.bz2 | tar xfj - --strip-components=1 -C htslib
    cd htslib && ./configure
    sudo make
    sudo make install

    #git clone https://github.com/samtools/htslib.git
    #sudo make -C htslib

fi
