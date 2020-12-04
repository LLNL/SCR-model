
#import psyco
#psyco.full()

##  SCR Model Code Version 0.1
##  Copyright (c) 2020, Lawrence Livermore National Security, LLC.
##  Produced at the Lawrence Livermore National Laboratory.
##  LLNL-CODE-817266. All rights reserved.
##  Please see the file LICENSE for the license for this code.
##

from math import * 



#Y(1, _) . base case
#X(L, c) . Y(L, c)
#Y(L, c) . if V[L-1] > 0:
#                Y(L-1, L-1) + Z(L-1, c)
#             else:
#                Y(L-1, c)
#Z(L, c) . (V[L]-1)*X(L, L) + X(L, c)

def valDebugPrint(val, p0, t0, pis, tis):
      str = "%s: p%s0=%4f t%s0=%4f \np%sis=[" % (val,val,p0,val, t0,val)
      for i in range(len(pis)):
        str += "%4f " % (pis[i])
      str += "]\nt%sis=[" % val
      for i in range(len(tis)):
        str += "%4f " % ( tis[i])
      str += "]"
      print str

XValCache = {}
YbaseCache = {}   
YValCache = {}   
ZValCache = {}   
     
class XValues:
    def __init__(self,pX0, tX0, pXis, tXis):
      self.pX0 = pX0
      self.tX0 = tX0
      self.pXis = pXis
      self.tXis = tXis
    def debugPrint(self):
      valDebugPrint("X", self.pX0, self.tX0, self.pXis, self.tXis)

class YValues:
    def __init__(self,pY0, tY0, pYis, tYis):
      self.pY0 = pY0
      self.tY0 = tY0
      self.pYis = pYis
      self.tYis = tYis
    def debugPrint(self):
      valDebugPrint("Y", self.pY0, self.tY0, self.pYis, self.tYis)


class ZValues:
    def __init__(self,pZ0, tZ0, pZis, tZis):
      self.pZ0 = pZ0
      self.tZ0 = tZ0
      self.pZis = pZis
      self.tZis = tZis
    def debugPrint(self):
      valDebugPrint("Z", self.pZ0, self.tZ0, self.pZis, self.tZis)

class RValues:
    def __init__(self,pR0, tR0, pRis, tRis):
      self.pR0 = pR0
      self.tR0 = tR0
      self.pRis = pRis
      self.tRis = tRis
    def debugPrint(self):
      valDebugPrint("R", self.pR0, self.tR0, self.pRis, self.tRis)

def p0(T, opts): 
    lam = opts.lam
    return  exp(-lam * T)

def t0(T, opts):
    return T

def pi(T, i, opts):
    lam = opts.lam
    lami = opts.lambdas[i]
    a =  (1.0 - exp(-lam * T))
    return ((lami*a ))/lam 


def ti(T, opts):
    lam = opts.lam
    top = 1.0 - (lam * T + 1.0)* exp(-lam * T)
    bottom = lam * (1.0 - exp(-lam * T))
    #print "T: %f lam: %f top: %f bottom %f" % (T,lam,top, bottom)
    return  top/bottom

def checkYbaseCache(c, opts):
    try:
       t = opts.t
       L = opts.L
       cC = opts.checkWriteTimes[c]
       lam = opts.lam
       Yvals = YbaseCache[(t, c, L, cC, lam)] 
       print "ybase match"
       return Yvals 
    except:
       return None

def addYbaseCache(c, opts, Yvals):
    t = opts.t
    L = opts.L
    cC = opts.checkWriteTimes[c]
    lam = opts.lam
    YbaseCache[(t, c, L, cC, lam)] = Yvals

def YbaseCase(c, opts):
    t = opts.t
    cws = opts.checkWriteTimes
    cC = cws[c]
    L = opts.L

    ret = checkYbaseCache(c, opts)
    if ret:
        return ret
    # compute probability  and time of exiting Y with no failures
    pY0 = p0(t + cC, opts)
    tY0 = t0(t + cC, opts)
   
    # compute probabilities of exiting Y at each failure level
    pYis = [0]# put a dummy in for index 0
    tYis = [0]
    for i in range(1,L+1):
        pYi = pi(t + cC, i, opts)
        tYi = ti(t + cC, opts)
        pYis.append(pYi)
        tYis.append(tYi)

    if pY0 == 0.0 :
       ret = YValues(pY0,tY0,pYis,tYis)
       addYbaseCache(c, opts, ret)
       return ret 

    if opts.debug:
       #print pY0
       #print pYis
       # the sum of all probabilities leaving a state should be 1.0
       total = (sum([pY0] + pYis))
       assert (fabs(total - 1.0) < opts.epsilon), ("%4f %4f" % (total, total - 1.0))
       # each individual ti should be less than t0
       for i in range(1,L+1):
           assert tY0 > tYis[i]

    ret = YValues(pY0,tY0,pYis,tYis)
    addYbaseCache(c, opts, ret)
    return ret

def RbaseCase(k, opts):
    # compute probabilities of exiting R  with no failures
    rs = opts.checkRecoverTimes
    rK = rs[k]
    L = opts.L
    pR0 = p0(rK,opts)
    tR0 = t0(rK, opts)

    # compute probabilities of exiting R  at each failure level
    pRis = [0] # put a dummy in for index 0
    tRis = [0]
    for i in range(1,L+1):
        pRi = pi(rK, i, opts)
        tRi = ti(rK, opts)
        pRis.append(pRi)
        tRis.append(tRi)
    #print "%s %s" % (pRis, tRis)

    if opts.debug:
       #print pR0
       #print pRis
       total = (sum([pR0] + pRis))
       assert (fabs(total - 1.0) < opts.epsilon), ("%4f %4f" % (total, total - 1.0))
       #print "Rbase %s %s" % (total, total - 1.0)
    ret = RValues(pR0,tR0,pRis,tRis)
    return ret



def checkZvalCache(k, c, v, X1r, X2r, opts):
    try:
       x1 = X1r
       pis =  X1r.pXis
       tis =  X1r.tXis
       zeros = True
       vals = zip(pis,tis)
       for p,t in vals:
           if p != 0 or t != 0:
              zeros = False
              break
       if zeros:
           x1 = 0
       Zvals = ZValCache[(k, c, v, x1, X2r)]
       print "zmatch"
       return Zvals
    except:
       return None

def addZvalCache(k, c, v, X1r, X2r, opts, Zvals):
    x1 = X1r
    pis =  X1r.pXis
    tis =  X1r.tXis
    zeros = True
    vals = zip(pis,tis)
    for p,t in vals:
        if p != 0 or t != 0:
           zeros = False
           break
    if zeros:
        x1 = 0
    ZValCache[(k, c, v, x1, X2r)] = Zvals

#Z(k, c) . (V[k]-1)*X(k, k) + X(k, c)
def Z(k, c, opts):
    if opts.debug:
          print "ENTER Z(%s, %s)" % (k, c)
    v = opts.V[k]
   
     
    L = opts.L
    # initialize to 0's for all values 
    X1r = XValues(0,0,[0 for i in range(L+1)],[0 for i in range(L+1)])
    # if there is more than one X state in this Z, first compute
    # the values for X(k, k)
    if (v - 1 > 0):
       X1r = X(k, k, opts)
       if opts.debug:
          X1r.debugPrint()
    # compute the values for X(k, c)
    X2r = X(k, c, opts) 
    if opts.debug:
       X2r.debugPrint()

    ret = checkZvalCache(k, c, v, X1r, X2r, opts)
    if ret:
        return ret

    # compute the probability and time of exiting Z with no failures
    pZ0 = pow(X1r.pX0,v-1) * X2r.pX0 
    tZ0 = (v-1) * X1r.tX0 + X2r.tX0

    # compute the probabilities and times of exiting Z at each failure level
    pZis = [0]
    tZis = [0]
    
    for i in range(1,L+1): # index starts with 1
        pZi = 0.0
        tZi = 0.0
        if i > k:
           # what should default values for pXis be?
           pZi = ((1.0 - pow(X1r.pX0,v-1))/(1.0 - X1r.pX0))*(X1r.pXis[i]) + pow(X1r.pX0,v-1) * X2r.pXis[i]
           
           B1 = (1.0 - pow(X1r.pX0,v-1))/(1 - X1r.pX0) *X1r.pXis[i]*X1r.tXis[i] 
           B2 = (X1r.pX0 - (v-1) * pow(X1r.pX0,v-1) + (v-2) * pow(X1r.pX0,v))
           B2 = B2/(pow((1-X1r.pX0),2)) * X1r.pXis[i] * X1r.tX0
           B = B1 + B2
           A = B + pow(X1r.pX0,v-1) * X2r.pXis[i]* ((v-1) * X1r.tX0 + X2r.tXis[i])
           tZi = A/pZi
        pZis.append(pZi)
        tZis.append(tZi)

    if pZ0 == 0.0:
       ret = ZValues(pZ0, tZ0, pZis, tZis)
       addZvalCache(k, c, v, X1r, X2r, opts, ret)
       return ret
    if opts.debug:
       # the sum of all probabilities leaving a state should be 1.0
       assert fabs(sum([pZ0] + pZis) - 1.0) < opts.epsilon
       # each individual ti should be less than t0
       for i in range(1,L+1):
           assert tZ0 > tZis[i]

    if opts.debug:
       print "LEAVE Z(%s, %s)" % (k, c)
    ret = ZValues(pZ0, tZ0, pZis, tZis)
    addZvalCache(k, c, v, X1r, X2r, opts, ret)
    return ret

def checkYvalCache(k, c, Zr, Yr, opts):
    try:
       L = opts.L
       print "Y: k, c, Zr, Yr" % (k, c, Zr, Yr)
       Yvals = YValCache[(k, c, Zr, Yr)]
       print "ymatch"
       return Yvals
    except:
       return None

def addYvalCache(k, c, Zr, Yr, opts, Yvals):
    YValCache[(k, c, Zr, Yr)] = Yvals


#Y(k, c) . if V[k-1] > 0:
#                Y(k-1, k-1) + Z(k-1, c)
#             else:
#                Y(k-1, c)

def Y(k, c, opts):
    if opts.debug:
         print "ENTER Y(%s, %s)" % (k, c)
 

    # if k == 1, we compute the values for the base case
    if(k == 1):
       Yr = YbaseCase(c, opts)
       if opts.debug:
           Yr.debugPrint()
           print "BASE LEAVE Y(%s, %s)" % (k, c)
       return  Yr
    # if there are no recovery states at level k-1, there is no Z state
    if(opts.V[k-1] == 0): 
       Yr = Y(k-1, c, opts)
       if opts.debug:
          Yr.debugPrint()
          print "LEAVE Y(%s, %s)" % (k, c)
       return Yr
    # otherwise, compute values for a Y and Z state
    L = opts.L
   
    Yr = Y(k-1, k-1, opts)
    if opts.debug:
       Yr.debugPrint()
    Zr = Z(k-1, c, opts)
    if opts.debug:
       Zr.debugPrint()

    ret = checkYvalCache(k, c, Zr, Yr, opts)
    if ret:
        return ret

    # compute the probability and time of exiting Y with no failures
    pY0 = Yr.pY0 * Zr.pZ0
    tY0 = Yr.tY0 + Zr.tZ0

    # compute the probabilities and times of exiting Y at each failure level
    pYis = [0]
    tYis = [0]
    for i in range(1,L+1):
        pYi = Yr.pYis[i] + Yr.pY0 * Zr.pZis[i]
        top = (Yr.pYis[i]*Yr.tYis[i] + Yr.pY0*Zr.pZis[i]*(Yr.tY0 + Zr.tZis[i]))
        tYi = top/pYi
        pYis.append(pYi)
        tYis.append(tYi)

    #if pY0 == 0.0:
    if pY0  < opts.epsilon:
       ret = YValues(pY0, tY0, pYis, tYis)
       addYvalCache(k, c, Zr, Yr, opts, ret)
       return ret

    if opts.debug:
       # the sum of all probabilities leaving a state should be 1.0
       print "sum: %s eps:%s" % (sum([pY0] + pYis), opts.epsilon)
       assert fabs(sum([pY0] + pYis) - 1.0) < opts.epsilon
       # each individual ti should be less than t0
       for i in range(1,L+1):
           assert tY0 > tYis[i]
       


    if opts.debug: 
        print "LEAVE Y(%s, %s)" % (k, c)
    ret = YValues(pY0, tY0, pYis, tYis)
    addYvalCache(k, c, Zr, Yr, opts, ret)
    return ret

def checkXvalCache(k, c, Yr, opts):
    try:
       Xvals = XValCache[(k, c, Yr)]
       print "xmatch"
       return Xvals
    except:
       return None

def addXvalCache(k, c, Yr, opts, Xvals):
    XValCache[(k, c, Yr)] = Xvals

#X(k, c) . Y(k, c)
def X(k, c, opts):
    if opts.debug:
       print "ENTER X(%s, %s)" % (k, c)

    L = opts.L
    # has Y and recovery
    Yr = Y(k, c, opts)
    if opts.debug:
       Yr.debugPrint()

    ret = checkXvalCache(k, c, Yr, opts)
    if ret:
       return ret

    r = opts.recoveryVals
    pRis = r[k].pRis
    pR0 = r[k].pR0
    tRis = r[k].tRis
    tR0 = r[k].tR0
    #print "pR0: %4f tR0: %4f pRis: %s tRis: %s" %(pR0,tR0, pRis, tRis)
    pYR = 0.0
    tYR = 0.0
    pRR = 0.0
    tRR = 0.0
    for i in range(1,k+1): # range(1,5+1) makes 1,2,3,4,5
      pYR += Yr.pYis[i] 
      tYR += Yr.pYis[i] * Yr.tYis[i]
 
    #print "before tYR %4f pYR %4f" % (tYR, pYR)
    tYR = tYR/pYR

    M = k-1
    if k == opts.L:
      M = k
    for i in range(1,M+1): # range(1,5+1) makes 1,2,3,4,5
      pRR += pRis[i]
      tRR += pRis[i] * tRis[i] 

    #print "%f %f" % (tYR, pYR)
    #print "before tRR %4f pRR %4f" % (tRR, pRR)

    if (k != 1):
       tRR = tRR/pRR

    # if too close to 0, just return 0
    if (pRR - 1.0 > opts.epsilon):
       pX0 = 0.0
       tX0 = 0.0
       pXis = [0 for i in range(1,L+1)]
       tXis = [0 for i in range(1,L+1)]
       ret = XValues(pX0, tX0, pXis, tXis)
       addXvalCache(k, c, Yr, opts, ret)
       return ret

    pRY = pR0/(1.0-pRR)
    tRY = tR0 + (pRR/(1.0-pRR)) * tRR

    #print "pRR: %4f tRR: %4f" % (pRR,tRR)
    #print "pYR: %4f tYR: %4f" % (pYR,tYR)
    #print "pRY: %4f tRY: %4f" % (pRY,tRY)
    PRis = [0]
    TRis = [0]
    for i in range(1,L+1):
       pRi = 0.0
       tRi = 0.0
       if i == k + 1:
          pRi = (pRis[k] + pRis[i])/(1.0 - pRR) 
          tRi = (pRis[k]*tRis[k] + pRis[i]*tRis[i])/(pRis[k]+pRis[i]) + pRR/(1.0-pRR) * tRR
       elif i > k + 1:
          pRi = pRis[i]/(1.0 - pRR)
          tRi = tRis[i] + pRR/(1.0 - pRR) * tRR
       PRis.append(pRi)
       TRis.append(tRi)
         

    # error condition, returns  all zeros
    if (pYR - 1.0 > opts.epsilon) or (pYR+Yr.pY0 - 1.0 > opts.epsilon) or (1.0 - pYR*pRY < opts.epsilon) :
       pX0 = 0.0
       tX0 = 0.0
       pXis = [0 for i in range(1,L+1)]
       tXis = [0 for i in range(1,L+1)]
       ret = XValues(pX0, tX0, pXis, tXis)
       addXvalCache(k, c, Yr, opts, ret)
       return ret

    pX0 = Yr.pY0/(1.0 - pYR*pRY)
    tX0 = Yr.tY0 + ((pYR*pRY)/(1.0 - pYR*pRY)) * (tYR + tRY)
    #print "Yr.pY0:%s  pYR:%s pRy:%s prod:%s pX0:%s" % (Yr.pY0, pYR, pRY, pYR*pRY, pX0)

    pXis = [0]
    tXis = [0]
    for i in range(1,L+1): 
       pXi = 0.0     
       tXi = 0.0     
       if i > k:
           pXi = (Yr.pYis[i] + pYR * PRis[i])/(1.0 - pYR*pRY)
           tXi_1 = Yr.pYis[i]*Yr.tYis[i] + pYR*PRis[i] * (tYR + TRis[i])
           tXi_1 = tXi_1/(Yr.pYis[i] + pYR*PRis[i])
           tXi_2 = (pYR*pRY)/(1.0 - pYR*pRY) * (tYR + tRY)
           tXi = tXi_1 + tXi_2
       pXis.append(pXi)
       tXis.append(tXi) 

    #print "pX0: %s" % pX0
    if pX0 < opts.epsilon:
       ret = XValues(pX0, tX0, pXis, tXis)
       addXvalCache(k, c, Yr, opts, ret)
       return ret
     

    if opts.debug:

       #print pX0
       #_print pXis
       # the sum of all probabilities leaving a state should be 1.0
       total = fabs(sum([pX0] + pXis))
       assert ((total - 1.0) < opts.epsilon), ("%4f %4f" % (total, total - 1.0))
       # each individual ti should be less than t0
       for i in range(1,L+1):
           assert tX0 > tXis[i]


    if opts.debug:
        print "LEAVE X(%s, %s)" % (k, c)
    ret = XValues(pX0, tX0, pXis, tXis)
    addXvalCache(k, c, Yr, opts, ret)
    return ret

