#!/bin/bash

source /home/jlazar/golemsolarwimp_py2.sh; python /data/user/jlazar/solar_WIMP/calc_xs_uncertainty.py --ch $1 -m $2 --xs error --ordering $3

