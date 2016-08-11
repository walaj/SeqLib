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
  cd snow_test
  
  #./configure --with-boost=${BOOST_ROOT}
  #make
  rm snow_test-snow-test.gcno ## dont get cov on test prog

  ./snow_test
  EXCL="-e snow_test/snow_test.cpp -e htslib -e bwa -e tools -e snow_test/snow_test_main.cpp -e R -e examples -e doxy -e src/deprecated.h -e config.h -e snow_test/config.h -e src/SnowTools/gzstream.h -e snow_test/config.h"
  cpp-coveralls -r ../ -t ${COVERALLS_TOKEN} ${EXCL} ##--dryrun
  #find ./ -type f -regex ".*/[a-z_]+.*gcov" -delete
  #cpp-coveralls -r ../ -t ${COVERALLS_TOKEN} ${EXCL} -n
  cd ..
#fi
