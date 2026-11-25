#!/bin/sh
#
# installs additional Ubuntu packages needed by the Intel build
#
# get packages
#
sudo apt-get update -y -o Acquire::Retries=3
sudo apt-get install -y -o Acquire::Retries=3 --no-install-recommends libfftw3-dev
