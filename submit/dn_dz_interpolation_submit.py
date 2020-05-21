import pycondor
import numpy as np
import sys
import glob

ms = [500,800,1000,2000,3000,5000,10000]
chs = [5,8,11]


path   = "/home/jlazar/condor_logs/dn_dz_interpolation/"
error  = "%s/error" % path
output = "%s/output" % path
log    = "%s/log" % path
submit = "%s/submit" % path

xlines = ["request_memory = (NumJobStarts is undefined) ? 2 * pow(2, 10) : 2048 * pow(2, NumJobStarts + 1)",
          "periodic_release = (HoldReasonCode =?= 21 && HoldReasonSubCode =?= 1001) || HoldReasonCode =?= 21",
          "periodic_remove = (JobStatus =?= 5 && (HoldReasonCode =!= 34 && HoldReasonCode =!= 21)) || (RequestMemory > 13192)"
         ]
dagman = pycondor.Dagman("dn_dz_interpolation", submit=submit, verbose=2)

run     = pycondor.Job("interp", 
                       "/data/user/jlazar/solar_WIMP/dn_dz_interpolation.sh", 
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

for m in ms:
    fluxfiles = ["/data/user/jlazar/solar_WIMP/data/qr_dn_dz/ch%d-m%d_dn_dz.npy" % (ch, m) for ch in chs]
    for fluxf in fluxfiles:
        for mcf in np.genfromtxt("/data/user/jlazar/solar_WIMP/mc_paths.txt", dtype=str):
            run.add_arg("%s %s %d" % (fluxf, mcf, m))
dagman.build()
