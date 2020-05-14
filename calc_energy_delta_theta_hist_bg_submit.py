import pycondor
import numpy as np
import sys

import argparse
##### SET UP COMMAND LINE ARGUMENTS #####
parser = parser = argparse.ArgumentParser()
parser.add_argument("--nt",
                    type=str,
                    help="neutrino type (nu or nuBar)"
                   )
parser.add_argument("--binning",
                    type=str,
                    help="bin width (f=0.5 degree, f=1.6 degrees)"
                   )
parser.add_argument("--bgparticle",
                    type=str,
                    help="bin width (f=0.5 degree, f=1.6 degrees)"
                   )
args = parser.parse_args()
bgParticle = args.bgparticle
nuType     = args.nt
binning    = args.binning
#args       = sys.argv
#bgParticle = args[1]
#nuType     = args[2]
#binning    = args[3]

path   = "/home/jlazar/condor_logs/energy_delta_theta_hist_bg"
error  = "%s/error" % path
output = "%s/output" % path
log    = "%s/log" % path
submit = "%s/submit" % path

runNs  = [ int(i) for i in np.linspace(0,9,10) ]
mcFile = "/data/user/jlazar/data/solar_WIMP/data/mcRecarray.npy"

outfile = "energy_delta_theta_hist_bg_%s_%s.out" % (bgParticle, nuType)
run     = pycondor.Job("energy_delta_theta_hist_bg_%s_%s" % (bgParticle, nuType), 
                       "/data/user/jlazar/solar_WIMP/get_energy_delta_theta_hist_bg.sh", 
                       error=error, 
                       output=output, 
                       log=log, 
                       submit=submit, 
                       universe="vanilla", 
                       notification="never"
                      )

for n in runNs:	
    run.add_arg("%s %s %s %s %s" % (n, bgParticle, nuType, binning, mcFile))
run.build()
#run.build_submit()
