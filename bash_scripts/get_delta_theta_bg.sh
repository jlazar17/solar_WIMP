#!/bin/bash
RUNN=$1
BGPARTICLE=$2
NUTYPE=$3
MCFILE=$4
OUTFILE=$5

eval `/cvmfs/icecube.opensciencegrid.org/py3-v4/setup.sh` python /data/user/jlazar/solar_WIMP/calc_delta_theta_bg.py ${RUNN} ${BGPARTICLE} ${NUTYPE} ${MCFILE} | tail -1 >> ${OUTFILE}
