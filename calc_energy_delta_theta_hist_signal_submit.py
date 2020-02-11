import pycondor
import numpy as np
from sys import argv as args

ch       = int(args[1])
m        = int(args[2])
factor   = float(args[3])
binning  = args[4] 

path   = "/home/jlazar/condor_logs/energy_delta_theta_hist_signal"
error  = "%s/error" % path
output = "%s/output" % path
log    = "%s/log" % path
submit = "%s/submit" % path

run     = pycondor.Job("ch%s_m%s_energy_delta_theta_hist_signal" % (ch, m), 
           "/data/user/jlazar/solar_WIMP/get_energy_delta_theta_hist_signal.sh", 
           error=error, 
           output=output, 
           log=log, 
           submit=submit, 
           universe="vanilla", 
           notification="never"
          )
run.add_arg("%s %s %s %s" % (ch, m, factor, binning))
run.build()

#run_ns = [ int(i) for i in np.linspace(0,24,25) ]
#
#
#for nt in nu_types:
#    run     = pycondor.Job("ch%s_m%s_%s_energy_delta_theta_hist_signal" % (ch, m, nt), 
#               "/data/user/jlazar/solar_WIMP/get_energy_delta_theta_hist_signal.sh", 
#               error=error, 
#               output=output, 
#               log=log, 
#               submit=submit, 
#               universe="vanilla", 
#               notification="never"
#              )
#    for n in run_ns:
#        run.add_arg("%s %s %s %s %s %s" % (n, nt, ch, m, factor, binning))
#    run.build()
