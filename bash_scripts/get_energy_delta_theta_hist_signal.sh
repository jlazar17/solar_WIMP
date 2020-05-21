#!/bin/bash
MCFILE=$1
FLUX=$2


source /home/jlazar/golemsolarwimp.sh; python /data/user/jlazar/solar_WIMP/calc_energy_delta_theta_hist_signal.py --mcfile $MCFILE --flux $FLUX
