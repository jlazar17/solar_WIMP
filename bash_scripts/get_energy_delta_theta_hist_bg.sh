#!/bin/bash
RUNN=$1
BGPARTICLE=$2
NUTYPE=$3
BINNING=$4

eval `/cvmfs/icecube.opensciencegrid.org/py3-v4/setup.sh` python /data/user/jlazar/solar_WIMP/calc_energy_delta_theta_hist_bg.py -n ${RUNN} --bgparticle ${BGPARTICLE} --nt ${NUTYPE} --binning ${BINNING}
