import pycondor
import numpy as np
import sys
import glob


path   = "/home/jlazar/condor_logs/conventional_flux_interpolation/"
error  = "%s/error" % path
output = "%s/output" % path
log    = "%s/log" % path
submit = "%s/submit" % path

xlines = ["request_memory = (NumJobStarts is undefined) ? 2 * pow(2, 10) : 2048 * pow(2, NumJobStarts + 1)",
          "periodic_release = (HoldReasonCode =?= 21 && HoldReasonSubCode =?= 1001) || HoldReasonCode =?= 21",
          "periodic_remove = (JobStatus =?= 5 && (HoldReasonCode =!= 34 && HoldReasonCode =!= 21)) || (RequestMemory > 13192)"
         ]
dagman = pycondor.Dagman("conventional_flux_interpolation", submit=submit, verbose=2)

run     = pycondor.Job("conv_interp", 
                       "/data/user/jlazar/solar_WIMP/conventional_flux_interpolation.sh", 
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

for mcf in np.genfromtxt("/data/user/jlazar/solar_WIMP/mc_paths.txt", dtype=str):
    run.add_arg("%s" % (mcf))
dagman.build()
