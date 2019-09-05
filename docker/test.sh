#!/bin/bash
echo running in ${PWD}
docker run -v ${PWD}:/mnt carlocreator/varsim /opt/varsim/tests/quickstart_test/quickstart.sh /mnt
