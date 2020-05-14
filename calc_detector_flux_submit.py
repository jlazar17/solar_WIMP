import pycondor
import numpy as np
import argparse

##### SET UP COMMAND LINE ARGUMENTS #####
parser = argparse.ArgumentParser()
#parser.add_argument("--mc_file",
#                    type=str,
#                    help="Monte Carlo file"
#                   )
parser.add_argument("-m",
                    type=int,
                    help="Dark matter mass"
                   )

args    = parser.parse_args()
m       = args.m 
chs     = [5, 8, 11]

path   = "/home/jlazar/condor_logs/detector_flux"
error  = "%s/error" % path
output = "%s/output" % path
log    = "%s/log" % path
submit = "%s/submit" % path

for ch in chs:
    outfile = "ch%s_m%s_detector_flux.out" % (ch, m)
    run     = pycondor.Job("ch%s_m%s_detector_flux" % (ch, m), 
               "/data/user/jlazar/solar_WIMP/get_detector_flux.sh", 
               error=error, 
               output=output, 
               log=log, 
               submit=submit, 
               universe="vanilla", 
               notification="never"
              )
    run.add_arg("%s %s" % (ch, m))
    run.build()
#run.build_submit()
