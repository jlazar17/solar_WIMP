#!/bin/bash
CH=$1
M=$2
XS=$3

eval `/cvmfs/icecube.opensciencegrid.org/py3-v4/setup.sh` python3 /data/user/jlazar/solar_WIMP/asimov_sens.py $CH $M $XS
