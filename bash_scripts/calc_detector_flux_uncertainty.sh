#!/bin/bash

CH=$1
M=$2
THETA12=$3
THETA23=$4
THETA13=$5
DELTAM12=$6
DELTAM13=$7
DELTA=$8

source /home/jlazar/golemsolarwimp_py2.sh; python /data/user/jlazar/solar_WIMP/calc_detector_flux_uncertainty.py --ch $CH -m $M --theta_12 $THETA12 --theta_23 $THETA23 --theta_13 $THETA13 --delta_m_12 $DELTAM12 --delta_m_13 $DELTAM13 --delta $DELTA
