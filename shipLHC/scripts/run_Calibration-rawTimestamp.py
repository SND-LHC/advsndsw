#!/usr/bin/env python
import ROOT,os,sys,getopt,subprocess,atexit,time
import Monitor
import vetoTimeCalibration

def pyExit():
       print("Make suicide until solution found for freezing")
       os.system('kill '+str(os.getpid()))
atexit.register(pyExit)

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("-f", "--inputFile", dest="fname", help="file name", type=str,default=None,required=False)
parser.add_argument("--server", dest="server", help="xrootd server",default=os.environ["EOSSHIP"])
parser.add_argument("-r", "--runNumber", dest="runNumber", help="run number", type=int)
parser.add_argument("-p", "--path", dest="path", help="run number",required=False,default="")
parser.add_argument("-P", "--partition", dest="partition", help="partition of data", type=int,required=False,default=-1)
parser.add_argument("-g", "--geoFile", dest="geoFile", help="geofile", required=True)
parser.add_argument("-b", "--heartBeat", dest="heartBeat", help="heart beat", default=10000,type=int)
parser.add_argument("-c", "--command", dest="command", help="command", default="")
parser.add_argument("-n", "--nEvents", dest="nEvents", help="number of events", default=-1,type=int)
parser.add_argument("-s", "--nStart", dest="nStart", help="first event", default=0,type=int)
parser.add_argument("-t", "--trackType", dest="trackType", help="DS or Scifi", default="DS")

options = parser.parse_args()
options.online = False
options.xCheck = False

import ctypes
from array import array
def FCN(npar, gin, f, par, iflag):
#calculate chisquare
       chisq  = 0
       for bar in range(nbars):
           for bar1 in [bar-1,bar,bar+1]:
               if bar1 < 0 or bar1 > nbars: continue
               if X[bar][bar1][0]==0: continue
               d = X[bar][bar1][0] - par[bar] - par[bar1+nbars]
               if X[bar][bar1][1]>0: chisq += d**2/X[bar][bar1][1]**2
       f.value = chisq
       return

import pickle
def minimize():
  s=1
  M.h['tdcCalibPlane'] = {10:{'L':{},'R':{}},11:{'L':{},'R':{}}}
  for side in ['L','R']:
     global X
     X = M.h['matrix'][side]
     global nbars
     nbars  = M.systemAndBars[s]
     npar = nbars*2
     ierflg    = ctypes.c_int(0)
     vstart  = array('d',[0]*npar)
     gMinuit = ROOT.TMinuit(npar)
     gMinuit.SetFCN(FCN)
     gMinuit.SetErrorDef(1.0)
     gMinuit.SetMaxIterations(10000)
     p = 0
     err = 1E-3
     for l in range(2):
        for o in range(nbars):
             name = "o"+str(10*l+o)
             gMinuit.mnparm(p, name, vstart[p], err, 0.,0.,ierflg)
             p+=1
     gMinuit.FixParameter(4)
     strat = array('d',[0])
     gMinuit.mnexcm("SET STR",strat,1,ierflg) # 0 faster, 2 more reliable
     gMinuit.mnexcm("SIMPLEX",vstart,npar,ierflg)
     gMinuit.mnexcm("MIGRAD",vstart,npar,ierflg)
     cor0 = ctypes.c_double(0)
     ecor0 = ctypes.c_double(0)
     cor1 = ctypes.c_double(0)
     ecor1 = ctypes.c_double(0)

     for bar in range(npar):
        if bar<7: l=0
        else:   l=1
        rc = gMinuit.GetParameter(bar,cor0,ecor0)
        M.h['tdcCalibPlane'][10+l][side][bar-7*l] = cor0.value

     for bar in range(nbars):
        for bar1 in [bar-1,bar,bar+1]:
            if bar1 < 0 or bar1 > nbars: continue
            if X[bar][bar1][0]==0: continue
            rc = gMinuit.GetParameter(bar,cor0,ecor0)
            rc = gMinuit.GetParameter(bar1+nbars,cor1,ecor1)
            d = X[bar][bar1][0] - (cor1.value - cor0.value)
            print("relation between veto0 %i and veto1 %i, before correction %6.4F and after %6.4F"%(bar,bar1,X[bar][bar1][0],d))

  with open('tdcVetoPlaneCalibration', 'wb') as fh:
     pickle.dump(M.h['tdcCalibPlane'], fh)


def mergeResiduals():
   M.h['exBar'] = M.h['exBar_04'].Clone('exBar')
   M.h['exBar'].Reset()
   for l in range(2):
      for bar in range(7):
            M.h['exBar'].Add(M.h['exBar_'+str(l)+str(bar)])

def mergeTimeDiffs():
   s = 1
   for tag in ['','X']:
    for side in ['L','R']:
      M.h['dtBar_Veto'+side+tag] = M.h['dtBar_Veto4_4LX'].Clone('dtBar_Veto'+side+tag)
      M.h['dtBar_Veto'+side+tag].Reset()
      for bar in range(nbars):
          for bar1 in [bar-1,bar,bar+1]:
                    if bar1<0 or bar1>nbars-1:  continue
                    key = M.sdict[s]+str(bar)+'_'+str(bar1)+side+tag
                    M.h['dtBar_Veto'+side+tag].Add(M.h['dtBar_'+key])

FairTasks = {}
M = Monitor.Monitoring(options,FairTasks)

Tasks = {}
Tasks['vetoTDCchannelcalibration']   = vetoTDCcalibration.vetoTDCchannelCalibration()
Tasks['vetoTDCplanecalibration']       = vetoTDCcalibration.vetoTDCplaneCalibration()
Tasks['vetoTDCchannelcalibrationX']   = vetoTDCcalibration.vetoTDCchannelCalibration()
Tasks['vetoTDCplanecalibrationX']       = vetoTDCcalibration.vetoTDCplaneCalibration()

if options.nEvents < 0 : 
    options.nEvents = M.eventTree.GetEntries()

for t in ['vetoTDCchannelcalibration']:
   options.xCheck = False
   Tasks[t].Init(options,M)
   for n in range(options.nStart,options.nStart+options.nEvents):
     event = M.GetEvent(n)
     Tasks[t].ExecuteEvent(M.eventTree)
   Tasks[t].Finalize()

for t in ['vetoTDCchannelcalibrationX']:
   options.xCheck = True
   Tasks[t].Init(options,M)
   for n in range(options.nStart,options.nStart+options.nEvents):
     event = M.GetEvent(n)
     Tasks[t].ExecuteEvent(M.eventTree)
   Tasks[t].Finalize()

for t in ['vetoTDCplanecalibration']:
   options.xCheck = False
   Tasks[t].Init(options,M)
   for n in range(options.nStart,options.nStart+options.nEvents):
     event = M.GetEvent(n)
     Tasks[t].ExecuteEvent(M.eventTree)
   Tasks[t].Finalize()

minimize()

for t in ['vetoTDCplanecalibrationX']:
   options.xCheck = True
   Tasks[t].Init(options,M)
   for n in range(options.nStart,options.nStart+options.nEvents):
     event = M.GetEvent(n)
     Tasks[t].ExecuteEvent(M.eventTree)
   Tasks[t].Finalize()



mergeResiduals()
mergeTimeDiffs()

