#!/bin/bash
CH=$1
MCHI=$2

source /home/jlazar/golemspace.sh
eval `/cvmfs/icecube.opensciencegrid.org/py2-v3/setup.sh`
python /data/user/jlazar/solar_WIMP/calc_detector_flux.py --ch $CH -m $MCHI
