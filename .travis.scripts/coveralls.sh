#!/bin/bash

echo "...running unit tests and code coverage"

if [ "$BUILDDOX" == "1" ] && [ "$TRAVIS_BRANCH" == "master" ] && [ "$TRAVIS_PULL_REQUEST" == "false" ]; 
then
  ## only install if not on home environment (eg travis ci)
    if [ -z "$REFHG19" ]; then
	sudo pip install cpp-coveralls
    fi
    cd seq_test

    ## download the test data
    mkdir test_data
    cd test_data
    wget -r -nH -nd -np -R index.html* https://data.broadinstitute.org/snowman/SeqLibTest/
    cd ..
    
    export LD_LIBRARY_PATH=${BOOST_ROOT}/lib:${LD_LIBRARY_PATH}
    echo "LD_LIBRARY_PATH: $LD_LIBRARY_PATH"
    ./configure --with-boost=${BOOST_ROOT}
    make CXXFLAGS='-DHAVE_C11=1 -std=c++11'
    ./seq_test
    
    EXCL="-e src/non_api -e seq_test/seq_test.cpp -e htslib -e bwa -e fermi-lite -e config.h -e seq_test/config.h -e seq_test/config.h -e src/jsoncpp.cpp -e src/json -e src/SeqLib/ssw.h -e src/SeqLib/ssw_cpp.h -e src/ssw.c -e src/ssw_cpp.cpp"
    cpp-coveralls -r ../ -t ${COVERALLS_TOKEN} ${EXCL} ##--dryrun
    cd ..
fi
