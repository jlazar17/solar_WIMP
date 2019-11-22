import pycondor
import numpy as np
import sys

args   = sys.argv
ch     = args[1]
mChi   = args[2]
nuType = args[3]

path   = "/home/jlazar/condor_logs/num_events_from_sun"
error  = "%s/error" % path
output = "%s/output" % path
log    = "%s/log" % path
submit = "%s/submit" % path

runNs = [ int(i) for i in np.linspace(0,9,10) ]
mcFile = "/data/user/jlazar/data/solarWIMP/data/mcRecarray.npy"

outfile = "exp_num_events_from_sun_ch%s_m%s_%s.txt" % (ch, mChi, nuType)
run     = pycondor.Job("events_from_sun_ch%s_m%s_%s" % (ch, mChi, nuType), 
                       "/data/user/jlazar/solar_WIMP/get_events_from_sun.sh", 
                       error=error, 
                       output=output, 
                       log=log, 
                       submit=submit, 
                       universe="vanilla", 
                       notification="never"
                      )

for n in runNs:	
    run.add_arg("%s %s %s %s %s %s" % (n, nuType, ch, mChi, mcFile, outfile))
run.build()
#run.build_submit()
