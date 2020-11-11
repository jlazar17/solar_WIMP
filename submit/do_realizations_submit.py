import pycondor
import numpy as np
import argparse

chs = [5,8,11]
ms  = [300,500,800,1000,3000,5000,8000,10000]

path   = "/home/jlazar/condor_logs/realizations/"
error  = "%s/error" % path
output = "%s/output" % path
log    = "%s/log" % path
submit = "%s/submit" % path
        
run     = pycondor.Job("realizations", 
               "/data/user/jlazar/solar_WIMP/bash_scripts/data_realizations.sh", 
               error=error, 
               output=output, 
               log=log, 
               submit=submit, 
               universe="vanilla", 
               notification="never"
              )

for ch in chs:
    for m in ms:
        for opt in ['00', '01']:
            run.add_arg("%s %s %s" % (ch, m, opt))
run.build()
