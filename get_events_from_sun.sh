#!/bin/bash
RUNN=$1
NUTYPE=$2
CH=$3
MCHI=$4
MCFILE=$5
OUTFILE=$6

eval `/cvmfs/icecube.opensciencegrid.org/py3-v4/setup.sh` python /data/user/jlazar/solar_WIMP/calc_events_from_sun.py ${RUNN} ${NUTYPE} ${CH} ${MCHI} ${MCFILE} | tail -1 >> ${OUTFILE}
