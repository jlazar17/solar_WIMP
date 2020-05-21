#!/bin/bash
MCFILE=$1


source /home/jlazar/golemsolarwimp.sh; python /data/user/jlazar/solar_WIMP/conventional_flux_interp.py $MCFILE
