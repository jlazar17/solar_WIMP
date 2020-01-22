import pycondor
import numpy as np
from sys import argv as args

ch = int(args[1])
m  = int(args[2])
nu_types = ["nu", "nuBar"]

path   = "/home/jlazar/condor_logs/e_d_theta_signal_rescale"
error  = "%s/error" % path
output = "%s/output" % path
log    = "%s/log" % path
submit = "%s/submit" % path

run_ns = [ int(i) for i in np.linspace(0,99,100) ]


for nt in nu_types:
    outfile = "ch%s_m%s_%s_energy_delta_theta_hist_signal.out" % (ch, m, nt)
    run     = pycondor.Job("ch%s_m%s_%s_energy_delta_theta_hist_signal" % (ch, m, nt), 
               "/data/user/jlazar/solar_WIMP/test.sh", 
               error=error, 
               output=output, 
               log=log, 
               submit=submit, 
               universe="vanilla", 
               notification="never"
              )
    for n in run_ns:
        run.add_arg("%s %s %s %s" % (n, nt, ch, m))
    run.build()
