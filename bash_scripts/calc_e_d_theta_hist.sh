#!/bin/bash
MCFILE=$1
FLUX=$2


source /home/jlazar/golemsolarwimp.sh; python /data/user/jlazar/solar_WIMP/calc_e_d_theta_hist.py $MCFILE $FLUX
