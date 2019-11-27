import pycondor
import numpy as np
import sys

args    = sys.argv
nuType  = args[1]
ch      = args[2]
mChi    = args[3]
binning = args[4]

path   = "/home/jlazar/condor_logs/energy_delta_theta_hist_signal"
error  = "%s/error" % path
output = "%s/output" % path
log    = "%s/log" % path
submit = "%s/submit" % path

runNs = [ int(i) for i in np.linspace(0,99,100) ]
mcFile = "/data/user/jlazar/data/solar_WIMP/data/mcRecarray.npy"

outfile = "ch%s_m%s_%s_energy_delta_theta_hist_signal.out" % (ch, mChi, nuType)
run     = pycondor.Job("ch%s_m%s_%s_energy_delta_theta_hist_signal" % (ch, mChi, nuType), 
                       "/data/user/jlazar/solar_WIMP/get_energy_delta_theta_hist_signal.sh", 
                       error=error, 
                       output=output, 
                       log=log, 
                       submit=submit, 
                       universe="vanilla", 
                       notification="never"
                      )

for n in runNs:	
    run.add_arg("%s %s %s %s %s %s %s" % (n, nuType, ch, mChi, binning, mcFile, outfile))
#run.build()
run.build_submit()
