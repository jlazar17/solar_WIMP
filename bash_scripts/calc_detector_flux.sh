#!/bin/bash

source /home/jlazar/programs/charon/environment.sh; python /data/user/jlazar/solar_WIMP/calc_detector_flux.py --ch $1 -m $2 --whichgen $3
