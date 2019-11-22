#!/bin/bash
CH=$1
MCHI=$2
OUTFILE=$3

source /home/jlazar/golemspace.sh
echo "words"
eval `/cvmfs/icecube.opensciencegrid.org/py2-v3/setup.sh`
echo "words2"
python /data/user/jlazar/solar_WIMP/calc_detector_flux.py ${RUNN} ${NUTYPE} ${CH} ${MCHI} ${MCFILE} | tail -1 >> ${OUTFILE}
