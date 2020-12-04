#!/usr/bin/env python

##  SCR Model Code Version 0.1
##  Copyright (c) 2020, Lawrence Livermore National Security, LLC.
##  Produced at the Lawrence Livermore National Laboratory.
##  LLNL-CODE-817266. All rights reserved.
##  Please see the file LICENSE for the license for this code.
##

from math import * 
import sys,os
from optparse import OptionParser

from modelSCR_core import *


#import psyco
#psyco.full()

class options:
    def __init__(self, L, compTime, checkWriteT, checkRecoverT, lambdas, lam, V, k, c):
       self.L = L
       self.t = compTime
       self.checkWriteTimes = checkWriteT 
       self.checkRecoverTimes = checkRecoverT
       self.lambdas = lambdas
       self.lam = lam
       self.V = V
       self.k = k
       self.c = c
       self.debug = False
       # to hold the precomputed values for recovery states. 
       # They don't change over 
       # the course of the execution, so no reason to recompute 
       # them each X state.
       self.recoveryVals = None
    def setRvals(self, rVals):
       self.recoveryVals = rVals
    def setLevelRatio(self,lr):
       self.levelRatio = lr
    def setOverheadFactor(self, o):
       self.overheadFactor = o
    def setFailureFactor(self, f):
       self.failureFactor = f
    def setDebug(self, debug):
       self.debug  = debug
    def setEpsilon(self, eps):
       self.epsilon = eps
       
# NOTE: all arrays contain one extra box of storage at index zero. It's a little wasteful, but it makes the code cleaner since all the math is done with indices starting at 1

def setupRun(L,c,k,V,t,cws,crs,frs,lam, debug, epsilon):
    opts = options(int(L),float(t), cws, crs, frs, lam, V, int(k), int(c))
    opts.setEpsilon(epsilon)
    opts.setDebug(debug)

    R = [None]
    for l in range(1,opts.L+1):
       r = RbaseCase(l, opts)
       R.append(r)
    opts.setRvals(R)
    
    return opts

def computeMetrics(Xvals, opts):
    mets = {}
    V = opts.V
    # number of intervals executed (application timesteps completed)
    # is the product of the number of times each checkpoint level
    # was taken, all checkpoint levels < K are taken V[i] + 1 times 
    # for each V[i+1] level checkpoint, e.g. if V = [0, 2, 3, 1]
    # the number of intervals executed = (2+1) * (3+1) * 1 
    # (remember index zero is ignored because we start counting levels at 1)
    # T:   T0 T1 T2 T3 T4 T5 T6 T7 T8 T9 T10 T11| T12 T13
    # L1:     C1 C1    C1 C1    C1 C2    C1  C1 |
    # L2:           C2       C2       C2        |
    # L3:  C3                                   | C3
    mets["time"]= Xvals.tX0
    l = len(V)-1
    count = 1.0
    for i in range(1,l):
       count *= (V[i] + 1) 
    count *= V[l]
    mets["intervals"] = count
    saved = count * opts.t
    mets["saved compute time"]= saved
    if(Xvals.tX0 != 0):
       efficiency = saved/Xvals.tX0
    else:
       efficiency = 0.0
    mets["efficiency"] = efficiency
    timePerInterval = Xvals.tX0/count
    mets["levels"] = opts.L
    mets["levelRatios"] = opts.V[1] #TODO there should be a better way to get this
    mets["intervalTime"] = opts.t
    mets["L-Level Overhead"] = opts.checkWriteTimes[opts.L] # only the highest level overhead
                                                            # changes, e.g. lustre
    mets["average interval time (t+cC)"]= timePerInterval
    
    return mets


def parseArgs():
    parser = OptionParser()
 
    # KM: for now, only support L = 3. Later, we'll add the ability to model
    # other numbers of levels.
    #parser.add_option("-l", "--levelsLow", dest="levelsLow", help="specify the lower bound on the number of checkpoint levels", default=3, type="int")
    #parser.add_option("--levelsStep", dest="levelsStep", help="specify the step value for the number of checkpoint levels", default=1, type="int")
    #parser.add_option("-L", "--levelsHigh", dest="levelsHigh", help="specify the upper bound on th number of checkpoint levels", default=3, type = "int")

    parser.add_option("--L1Low", dest="L1Low", help="specify the lower bound on the number of L1 checkpoint levels", default=2, type="int")
    parser.add_option("--L1Step", dest="L1Step", help="specify the step value for the number of L1 checkpoint levels", default=1, type="int")
    parser.add_option( "--L1High", dest="L1High", help="specify the upper bound on th number of L1 checkpoint levels", default=2, type = "int")
    parser.add_option("--L2Low", dest="L2Low", help="specify the lower bound on the number of L2 checkpoint levels", default=2, type="int")
    parser.add_option("--L2Step", dest="L2Step", help="specify the step value for the number of L2 checkpoint levels", default=1, type="int")
    parser.add_option( "--L2High", dest="L2High", help="specify the upper bound on th number of L2 checkpoint levels", default=2, type = "int")

    parser.add_option("-t", "--checkIntervalLow",dest="checkIntervalLow", help="specify the lower bound on checkpoint interval, i.e. application compute time", default=10, type="int")
    parser.add_option("-s", "--checkIntervalStep",dest="checkIntervalStep", help="specify the step value for the increase in checkpoint interval length", default=1, type="int")
    parser.add_option("-T", "--checkIntervalHigh",dest="checkIntervalHigh", help="specify the highest value for checkpoint interval length", default=10, type="int")
#
    parser.add_option("-f", "--failureFactor", dest="failureFactor", help="specify the factor to increase or decrease the failure rates (multiplicative)", default=1.0, type="float")
    # KM: removed the options to iterate over failure factors for a couple 
    # reasons
    # 1) because we weren't stepping evenly over the space, i.e. we used
    # factors of 1, 2, 10, 50
    # 2) in practice, to keep the experiments neater (and so they would 
    # finish running in a timely fashion), we just did one 
    # factor at a time 

    # KM: removed baseFailureRate because we were using real data intead of 
    # generating some made up numbers based on a base rate
    #parser.add_option("--baseFailureRate", dest="baseFailureRate", help="specify  the base failure rate ", default=1.0, type="float")

    #parser.add_option("-f", "--failureFactorLow", dest="failureFactorLow", help="specify lower bound on factor to increase or decrease the failure rates (multiplicative)", default=1.0, type="float")
    #parser.add_option("--failureFactorStep", dest="failureFactorStep", help="specify the step value for the increase in the factor to increase or decrease the failure rates (multiplicative)", default=1.0, type="float")
    #parser.add_option("-F", "--failureFactorHigh", dest="failureFactorHigh", help="specify upper bound on factor to increase or decrease the failure rates (multiplicative)", default=1.0, type="float")

    parser.add_option("-o", "--overheadFactor", dest="overheadFactor", help="specify the factor to increase or decrease the recover/write checkpoint overhead (multiplicative). Checkpoint recover and write times have the same value at this time.", default=1.0, type="float")
    # KM: removed the options to iterate over overhead factors for a couple
    # reasons
    # 1) because we weren't stepping evenly over the space, i.e. we used
    # factors of 1, 2, 10, 50
    # 2) in practice, to keep the experiments neater (and so they would
    # finish running in a timely fashion), we just did one
    # factor at a time

    # KM: removed baseOverheadRate because we were using real data intead of
    # generating some made up numbers based on a base rate
    #parser.add_option("--baseOverheadTime", dest="baseOverheadTime", help="specify the base value for recover/write checkpoint overhead. Checkpoint recover and write times have the same value at this time.", default=1.0, type="float")

    #parser.add_option("-o", "--overheadFactorLow", dest="overheadFactorLow", help="specify lower bound on factor to increase or decrease the recover/write checkpoint overhead (multiplicative). Checkpoint recover and write times have the same value at this time.", default=1.0, type="float")
    #parser.add_option("--overheadFactorStep", dest="overheadFactorStep", help="specify the step value for the increase in the factor to increase or decrease the checkpoint overhead (multiplicative)", default=1.0, type="float")
    #parser.add_option("-O", "--overheadFactorHigh", dest="overheadFactorHigh", help="specify upper bound on factor to increase or decrease the recover/write checkpoint overhead (multiplicative). Checkpoint recover and write times have the same value at this time.", default=1.0, type="float")
#
    parser.add_option("--machine", dest="machine", help="specify what machine to emulate. choices: atlas, coastal, made_up (made_up means not real data)")
    parser.add_option("-e", "--epsilon", dest="epsilon", help="specify the epsilon used in debugging checks", default=0.000001, type="float")
    parser.add_option("-d", "--debug", dest="debug", help="turn on printing of debugging statements", default=False, action="store_true")

    parser.add_option("--expDir", dest="expDir", help="specify the directory to write the data files to")
    parser.add_option("--expPrefix", dest="expPrefix", help="specify the a tag to be prepended to the data files", default="test")
    parser.add_option("--fullResults", dest="writeFullResults", help="turn on writing of full results (note: this can generate a lot of data, depending on the experiment)", default=False, action="store_true")
    (argOptions, args) = parser.parse_args()
    if not argOptions.expDir:
       print "you need to specify the output directory with the --expDir flag"
       sys.exit()
    else:
       cmd = "mkdir -p %s" % argOptions.expDir
       os.system(cmd)

    return argOptions

def buildV(L, ratios):
   # all arrays have a prepended zero, because ignoring index 0
   # allowing each level to have a different v[i], so pass in a 
   # list of counts for each level, e.g. ratios = [10,2,1], then 
   # V = [0, 10, 2, 1]
   V = [0]  + ratios
   return V

def buildOverheadTimes(overheadFactor, overheadTimes):
   # this simply takes an array of overhead times and multiplies each entry
   # by the overhead factor. 
   # KM: I don't believe this is used anymore, because we decided to model
   # increasing only the highest level overhead (e.g. Lustre)
   o = overheadFactor
   ovs = [i*o for i in overheadTimes]
   return ovs

def buildFailureRates(failFactor, failureRates):
   # this simply takes an array of failure rates and multiplies each entry
   # by the failure factor. 
   f = failFactor
   failRates = [i*f for i in failureRates]
   return failRates

def doExperiment(whatExp, expVector, level, lastCheckLevel, levelVector, intTime, overhdTimes,  failureRates, debug, epsilon):
 
   # only allows one parameter to vary at a time, specified by "whatExp"
   # all other values are used as they are passed in

   # KM: most used mode is "intervalTime", so it is the most tested

   t = intTime
   failRates = failureRates
   lam = sum(failRates)
   L = level
   k = level
   c = lastCheckLevel
   V = levelVector
   overheadTimes = overhdTimes
   checkWriteTimes = overheadTimes
   checkReadTimes = checkWriteTimes
 
   results = []
   #print expVector
   for v in expVector:
      if whatExp == "intervalTime":
         t = v
      elif whatExp == "levels":
         L = v
         k = v
      elif whatExp == "levelRatios":
         V = buildV(L, v)
      elif whatExp == "overheads":
         checkWriteTimes = buildOverheadTimes(v, overheadTimes)
         checkReadTimes = checkWriteTimes
      elif whatExp == "failRates":
         failRates = buildFailureRates(v, failureRates)
         lam = sum(failRates)
      else:
         print " I don't recognize experiment type: %s. FAIL." % whatExp
         sys.exit() 
      #print "running with L:%s c:%s k:%s V:%s t:%s " % (L,c,k,V,t)
      opts =  setupRun(L,c,k,V,t,checkWriteTimes,checkReadTimes, failRates, lam, debug, epsilon)
          
      Xr = X(opts.k, opts.c, opts)
      mets = computeMetrics(Xr,opts)
      results.append(mets)

   return results

def writeExperimentResults(results,dir, prefix, whatExp, low, high, step, L, t, l1,l2, overheads, o, failrates, f):
   filename = "%s/%s_%s_%s_%s_%s_%s_%s_%s_%s_%s_%s.data" % (dir,prefix, whatExp, low, high, step, L, t, l1,l2, o, f)
   print filename
   f = open(filename, 'w')
   f.write("# Experiment Type: %s\n" % whatExp)
   f.write("# Low: %s High: %s Step: %s\n" % (low, high, step)) 
   f.write("# Level: %s checkInterval: %s levelRatios: %s %s 1\n" % (L, t, l1,l2))
   f.write("# overheads: %s\n" % overheads)
   f.write("# failure rates: %s\n" % failrates)

   metrics = results[0].keys()
   metrics.sort()
   for m in metrics:
      f.write("%s" % m)
      if metrics.index(m) != (len(metrics)-1):
         f.write("\t")
   f.write("\n")

   for mets in results:
      for m in metrics:
           f.write("%s" % mets[m])
           if metrics.index(m) != (len(metrics)-1):
               f.write("\t")
      f.write("\n")

   f.close()

def getMaxEfficiencyAndTime(effs, times):
       maxeff = max(effs)
       loweff = maxeff - maxeff * 0.05
       higheff = maxeff + maxeff * 0.05
       loweffindx = 0
       maxeffindx = 0
       higheffindx = 0
       i = 0
       while i < len(effs)-1 and effs[i] < loweff:
          i += 1
       loweffindx = i
       while i < len(effs)-1 and effs[i] != maxeff:
          i += 1
       maxeffindx = i
       while i < len(effs)-1 and effs[i] < higheff:
          i+= 1
       higheffindx = i

       avgeff = sum(effs[loweffindx:higheffindx+1])/(higheffindx+1-loweffindx)
       avgtime = sum(times[loweffindx:higheffindx+1])/(higheffindx+1-loweffindx)
       return (maxeff, loweff, higheff, avgeff, times[maxeffindx],times[loweffindx], times[higheffindx], avgtime)


def main():
   go()

def go():
   opts = parseArgs()
   outDir = opts.expDir
   

   #KM: for now only suports L = 2 or 3, leaving this in for future work 
   # will have to build L1, L2, ..., Ln vectors for each L
   #Ls = [i for i in range(opts.levelsLow,opts.levelsHigh+1,opts.levelsStep)] 

   L1 = [i for i in range(opts.L1Low,opts.L1High+1,opts.L1Step)]
   L2 = [i for i in range(opts.L2Low,opts.L2High+1,opts.L2Step)]
   # L3 always equals 1
   
   checkIntervalLow = opts.checkIntervalLow
   checkIntervalHigh = opts.checkIntervalHigh
   checkIntervalStep = opts.checkIntervalStep
   timeVector = [i for i in range(checkIntervalLow,checkIntervalHigh+1,checkIntervalStep)]
   
   t = 1.0
   
   L = 3
   Ls = [L]
   # KM: for the real machines, we only have 2 levels of checkpointing, so
   # we'll model that by having 0 L3 checkpoints
   # and there's just 1 level 2 checkpoint 
   # the "made_up" machine has 3 levels
   if opts.machine == "atlas":
     baseFailRts = [0, 7.78e-6, 5.19e-7] # atlas
     oTimes = [0, 9, 438]
     L = 2
     L2 = [1]
   elif opts.machine == "coastal":
     baseFailRts = [0, 8.54e-7, 2.01e-7] # coastal
     oTimes = [0, 15, 1835 ]
     L = 2
     L2 = [1]
   elif opts.machine == "made_up":
     baseFailRts = [0, 7.78e-6, 5.19e-7, 3.0e-8] # made up 
     oTimes = [0, 9, 50, 438] # made up
   #elif opts.machine == "hera": # don't have hera data yet
     #baseFailRts = [0, 0, ] #hera
   else:
       print "Unknown machine specified: %s. Please choose one of: atlas, coastal" % opts.machine
       sys.exit()
   
   Ls = [L]
   
   for l in Ls: 
       o = opts.overheadFactor
       ovrTimes = oTimes
       ovrTimes[L] = o * oTimes[L] # only change highest level ovehread, e.g. Lustre
       f = opts.failureFactor
       failRts = buildFailureRates(f, baseFailRts)
       datas = []
       for l1 in L1:
         for l2 in L2:
           if opts.machine == "made_up":
              Levs = [l1,l2,1]
           else:
              Levs = [l1,1]
           levelVec = buildV(L, Levs)
           results = doExperiment("intervalTime", timeVector, L, L, levelVec, timeVector[0], ovrTimes, failRts, opts.debug, opts.epsilon)
           if opts.writeFullResults:
               print "writing"
               writeExperimentResults(results, outDir, opts.expPrefix, "intervalTime", checkIntervalLow, checkIntervalHigh, checkIntervalStep, L, t, l1, l2, ovrTimes, o, failRts, f)
           effs = [m['efficiency'] for m in results]
           times = [m['intervalTime'] for m in results]
           data = getMaxEfficiencyAndTime(effs,times)
           datas.append((l1,l2,data))
   

       filename = "%s/maxes.L%s.o%s.f%s.data" % (outDir,L,o,f)
       lines = "Max Efficency\tLow Efficiency\tHigh Efficiency\tAvg Efficiency\tMax Time\tLow Time\tHigh Time\tAvg Time\tL1\tL2\n"
       for l1,l2,(me,le,he,ae,mt,lt,ht,at) in datas:
           lines +=  "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (me, le, he, ae, mt,lt,ht,at, l1, l2)
       f = open(filename,'w')
     
       f.write(lines)
       f.close()

   

if __name__ == "__main__":
    main()

