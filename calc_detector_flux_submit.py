import pycondor
import numpy as np
from sys import argv as args

m = int(args[1])
chs      = ["bb", "WW", "tautau"]

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
    run.add_arg("%s %s %s" % (ch, m, outfile))
    run.build()
#run.build_submit()
