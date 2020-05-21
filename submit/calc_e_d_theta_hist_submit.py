import pycondor
import numpy as np
from sys import argv as args
from controls import chs, ms, mcfiles

path   = "/home/jlazar/condor_logs/e_d_theta_hist"
error  = "%s/error" % path
output = "%s/output" % path
log    = "%s/log" % path
submit = "%s/submit" % path

xlines = ["request_memory = (NumJobStarts is undefined) ? 2 * pow(2, 10) : 2048 * pow(2, NumJobStarts + 1)",
          "periodic_release = (HoldReasonCode =?= 21 && HoldReasonSubCode =?= 1001) || HoldReasonCode =?= 21",
          "periodic_remove = (JobStatus =?= 5 && (HoldReasonCode =!= 34 && HoldReasonCode =!= 21)) || (RequestMemory > 13192)"
         ]

dagman = pycondor.Dagman("e_d_theta_hist_dag", submit=submit, verbose=2)
run     = pycondor.Job("e_d_theta_hist", 
           "/data/user/jlazar/solar_WIMP/calc_e_d_theta_hist.sh", 
           error=error, 
           output=output, 
           log=log, 
           submit=submit, 
           universe="vanilla", 
           notification="never",
           dag=dagman,
           verbose=2,
           extra_lines=xlines
          )

fluxes = [f"ch{ch}-m{m}" for ch in chs for m in ms]
print(fluxes)

for mcfile in mcfiles:
    for flux in fluxes:
       run.add_arg(f"{mcfile} {flux}") 

dagman.build()
