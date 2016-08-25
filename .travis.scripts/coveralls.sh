#!/bin/bash

echo "...running unit tests and code coverage"
#reuse -q .python-3.5.1
#.zlib-1.2.8

# Note that this only works if the tests were built using --coverage for
# compile and link flags!
#if [ "$CXX" == "g++" ];
#then
  ## only install if not on home environment (eg travis ci)
  if [ -z "$REFHG19" ]; then
    sudo pip install cpp-coveralls
  fi
  cd seq_test
  
  export LD_LIBRARY_PATH=${BOOST_ROOT}/lib:${LD_LIBRARY_PATH}
  echo "LD_LIBRARY_PATH: $LD_LIBRARY_PATH"
  ./configure --with-boost=${BOOST_ROOT}
  make
  ./seq_test

  EXCL="-e src/non_api -e seq_test/seq_test.cpp -e htslib -e bwa -e fermi-lite -e config.h -e seq_test/config.h -e seq_test/config.h -e src/jsoncpp.cpp -e src/json -e src/SeqLib/ssw.h -e src/SeqLib/ssw_cpp.h -e src/ssw.c -e src/ssw_cpp.cpp"
  cpp-coveralls -r ../ -t ${COVERALLS_TOKEN} ${EXCL} ##--dryrun
  cd ..
#fi
