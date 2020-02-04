#!/bin/bash
CH=$1
MCHI=$2
FACTOR=$3
BINNING=$4

#eval `/cvmfs/icecube.opensciencegrid.org/py3-v4/setup.sh` python /data/user/jlazar/solar_WIMP/calc_energy_delta_theta_hist_signal.py ${RUNN} ${NUTYPE} ${CH} ${MCHI} ${MCFILE} | tail -1 >> ${OUTFILE}
eval `/cvmfs/icecube.opensciencegrid.org/py3-v4/setup.sh` python /data/user/jlazar/solar_WIMP/calc_energy_delta_theta_hist_signal.py --ch ${CH} -m ${MCHI} -f ${FACTOR} --binning ${BINNING}
