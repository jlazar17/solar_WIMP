#!/bin/bash
RUNN=$1
NUTYPE=$2
CH=$3
MCHI=$4
FACTOR=$5
BINNING=$6

#eval `/cvmfs/icecube.opensciencegrid.org/py3-v4/setup.sh` python /data/user/jlazar/solar_WIMP/calc_energy_delta_theta_hist_signal.py ${RUNN} ${NUTYPE} ${CH} ${MCHI} ${MCFILE} | tail -1 >> ${OUTFILE}
eval `/cvmfs/icecube.opensciencegrid.org/py3-v4/setup.sh` python /data/user/jlazar/solar_WIMP/calc_energy_delta_theta_hist_signal.py -n ${RUNN} --nt ${NUTYPE} --ch ${CH} -m ${MCHI} -f ${FACTOR} --binning ${BINNING}
