#!/bin/bash

FLUXFILE=$1
MCFILE=$2
M=$3
source /home/jlazar/golemsolarwimp.sh; python /data/user/jlazar/solar_WIMP/dn_dz_interpolation.py --fluxfile ${FLUXFILE} --mcfile ${MCFILE} -m ${M}
