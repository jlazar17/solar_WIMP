#!/bin/bash
source /home/jlazar/golemsolarwimp_py2.sh; python /data/user/jlazar/solar_WIMP/calc_total_uncertainty.py --ch $1 -m $2 --ordering $3 -n 2000
