#!/bin/bash

FLUXFILE=$1
MCFILE=$2
LICFILE=$3

source /home/jlazar/golemsolarwimp.sh; python /data/user/jlazar/solar_WIMP/do_everything.py --fluxfile $FLUXFILE --mcfile $MCFILE --LICFile $LICFILE
