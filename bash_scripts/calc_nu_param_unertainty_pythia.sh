#!/bin/bash

source /home/jlazar/golemsolarwimp_py2.sh; python /data/user/jlazar/solar_WIMP/calc_nu_param_unertainty.py --ch $1 -m $2 --ordering $3 -n 2000 --pythia True
