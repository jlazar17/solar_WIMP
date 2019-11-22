import pycondor
import numpy as np
import sys

args       = sys.argv
bgParticle = args[1]
nuType     = args[2]

path   = "/home/jlazar/condor_logs/energy_delta_theta_hist_bg_course"
error  = "%s/error" % path
output = "%s/output" % path
log    = "%s/log" % path
submit = "%s/submit" % path

runNs  = [ int(i) for i in np.linspace(0,9,10) ]
mcFile = "/data/user/jlazar/data/solarWIMP/data/mcRecarray.npy"

outfile = "energy_delta_theta_hist_bg_course_%s_%s.out" % (bgParticle, nuType)
run     = pycondor.Job("energy_delta_theta_hist_bg_course_%s_%s" % (bgParticle, nuType), 
                       "/data/user/jlazar/solar_WIMP/get_energy_delta_theta_hist_bg_course.sh", 
                       error=error, 
                       output=output, 
                       log=log, 
                       submit=submit, 
                       universe="vanilla", 
                       notification="never"
                      )

for n in runNs:	
    run.add_arg("%s %s %s %s %s" % (n, bgParticle, nuType, mcFile, outfile))
#run.build()
run.build_submit()
