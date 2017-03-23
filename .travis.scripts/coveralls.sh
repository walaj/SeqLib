#!/bin/bash

echo "...running unit tests and code coverage"

if [ "$COMPILER" == "g++-4.9" ] && [ "$TRAVIS_BRANCH" == "master" ] && [ "$TRAVIS_PULL_REQUEST" == "false" ]; 
then
  ## only install if not on home environment (eg travis ci)
    if [ -z "$REFHG19" ]; then
	sudo pip install cpp-coveralls
    fi
    cd seq_test

    ## download the latest matched gcov
    sudo update-alternatives --install /usr/bin/gcov gcov /usr/bin/gcov-4.9 90
    sudo ln -sf /usr/bin/gcov-4.9 /usr/bin/gcov
    GCOV_VERSION=`gcov --version`
    echo "GCOV version $GCOV_VERSION"

    ## download the test data
    mkdir test_data
    cd test_data
    wget -r -nH -nd -np -R index.html* https://data.broadinstitute.org/snowman/SeqLib/
    cd ..
    
    export LD_LIBRARY_PATH=${BOOST_ROOT}/lib:${LD_LIBRARY_PATH}
    echo "LD_LIBRARY_PATH: $LD_LIBRARY_PATH"
    ./configure --with-boost=${BOOST_ROOT}
    make CXXFLAGS='-DHAVE_C11=1 -std=c++11' CXX=$COMPILER
    ./seq_test 1> stdout.log 2> stderr.log
    tail -n stderr.log 
    
    EXCL="-e src/non_api -e seq_test/seq_test.cpp -e htslib -e bwa -e fermi-lite -e config.h -e seq_test/config.h -e seq_test/config.h -e src/jsoncpp.cpp -e src/json -e src/SeqLib/ssw.h -e src/SeqLib/ssw_cpp.h -e src/ssw.c -e src/ssw_cpp.cpp -e SeqLib/aho_corasick.hpp -e json/json.h -e SeqLib/ssw.h -e SeqLib/ssw_cpp.h"
    cpp-coveralls -r ../ -t ${COVERALLS_TOKEN} ${EXCL} ##--dryrun
    cd ..
fi
