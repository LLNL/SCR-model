#!/usr/bin/env python

##  SCR Model Code Version 0.1
##  Copyright (c) 2020, Lawrence Livermore National Security, LLC.
##  Produced at the Lawrence Livermore National Laboratory.
##  LLNL-CODE-817266. All rights reserved.
##  Please see the file LICENSE for the license for this code.
##

import os,sys,popen2
from glob import glob

driverName = "./modelSCR.py"
regressionFilesDir = "regressionFiles"
testDir = "regressTestingDir"

machine = "made_up"

L1low = 0
L1high = 1000
L1step = 200

L2low = 0
L2high = 1000
L2step = 200

tlow = 0
thigh=40000
tstep=500

Os = [1,2,50]
Fs = [1,2,50]

def my_exec2(cmd1,cmd2):
     """Executes cmd1 and gives it arguments in cmd2 (cmd2 is just a string
        containing the arguments)
        Returns stdout, stderr and the return status of cmd1
     """
     Cmd = cmd1.strip() + ' ' + cmd2
     #print "my_exec2 cmd is " + Cmd
     x = popen2.Popen3(Cmd, True)
     x.tochild.close()
     childout = x.fromchild.read().strip()
     childerr = x.childerr.read().strip()
     retval = x.wait()
     status = os.WEXITSTATUS(retval)
     #if (status != 0):
        #print "exec failed for " + Cmd
     return (childout, childerr, retval)


for o in Os:
  for f in Fs:
     cmd = "%s --L1Low %s --L1High %s --L1Step %s --L2Low %s --L2High %s --L2Step %s -t %s -T %s -s %s -f %s -o %s --machine %s --expDir %s --fullResults" % (driverName, L1low, L1high, L1step, L2low, L2high, L2step, tlow, thigh, tstep, f, o, machine, testDir)
     print cmd
     #os.system(cmd)
     out,err,retval = my_exec2(cmd,"")
     if retval != 0:
        print "FAIL. run failed for O:%s F:%s" % (o,f)
        print out
        print err
        sys.exit()


files = glob("%s/*" % testDir)

for file in files:
   rootname = file.lstrip(testDir)
   rootname = "%s%s" % (regressionFilesDir,rootname)

   cmd = "diff %s %s" % (file, rootname)
   #print cmd
   #os.system(cmd)
   out,err,retval = my_exec2(cmd,"")

   if retval != 0:
      print "FAIL. diff failed for %s" % file
      print out
      sys.exit()

print "SUCCESS! All tests PASSED!"
print ""

