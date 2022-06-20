#!/usr/bin/env python
import ROOT,os,sys,subprocess,atexit,time
import Monitor
import vetoTimeCalibration

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("-f", "--inputFile", dest="fname", help="file name", type=str,default=None,required=False)
parser.add_argument("--server", dest="server", help="xrootd server",default=os.environ["EOSSHIP"])
parser.add_argument("-r", "--runNumber", dest="runNumber", help="run number", type=int)
parser.add_argument("-p", "--path", dest="path", help="run number",required=False,default="")
parser.add_argument("-P", "--partition", dest="partition", help="partition of data", type=int,required=False,default=-1)
parser.add_argument("-g", "--geoFile", dest="geoFile", help="geofile", required=True)
parser.add_argument("-b", "--heartBeat", dest="heartBeat", help="heart beat", default=10000,type=int)
parser.add_argument("-c", "--command", dest="command", help="Time calibration / Time walk", default="TC")
parser.add_argument("-n", "--nEvents", dest="nEvents", help="number of events", default=-1,type=int)
parser.add_argument("-s", "--nStart", dest="nStart", help="first event", default=0,type=int)
parser.add_argument("-t", "--trackType", dest="trackType", help="DS or Scifi", default="DS")

options = parser.parse_args()
options.online = False
options.xCheck = False


def mergeResiduals():
   M.h['exBar'] = M.h['exBar_04'].Clone('exBar')
   M.h['exBar'].Reset()
   for l in range(2):
      for bar in range(7):
            M.h['exBar'].Add(M.h['exBar_'+str(l)+str(bar)])

def mergeTimeDiffs():
   s = 1
   nbars = M.systemAndBars[s]
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
if options.nEvents < 0 : 
    options.nEvents = M.eventTree.GetEntries()

Tasks = {}

if options.command=="TC":
   Tasks['vetoTimechannelcalibration']   = vetoTimeCalibration.vetoTDCchannelCalibration()
   Tasks['vetoTimeplanecalibration']       = vetoTimeCalibration.vetoTDCplaneCalibration()
   Tasks['vetoTimechannelcalibrationX']   = vetoTimeCalibration.vetoTDCchannelCalibration()
   Tasks['vetoTimeplanecalibrationX']       = vetoTimeCalibration.vetoTDCplaneCalibration()

   for t in ['vetoTimechannelcalibration']:
      print('start',t)
      options.xCheck = False
      Tasks[t].Init(options,M)
      for n in range(options.nStart,options.nStart+options.nEvents):
        event = M.GetEvent(n)
        Tasks[t].ExecuteEvent(event)
        event.Reco_MuonTracks.Delete()
      Tasks[t].Finalize()
      print('end',t)

   for t in ['vetoTimechannelcalibrationX']:
      print('start',t)
      options.xCheck = True
      Tasks[t].Init(options,M)
      for n in range(options.nStart,options.nStart+options.nEvents):
        event = M.GetEvent(n)
        Tasks[t].ExecuteEvent(event)
        event.Reco_MuonTracks.Delete()
      Tasks[t].Finalize()
      print('end',t)
   for t in ['vetoTimeplanecalibration']:
      print('start',t)
      options.xCheck = False
      Tasks[t].Init(options,M)
      for n in range(options.nStart,options.nStart+options.nEvents):
        event = M.GetEvent(n)
        Tasks[t].ExecuteEvent(event)
        event.Reco_MuonTracks.Delete()

      Tasks[t].Finalize()
      print('end',t)
      print('start minimization')
      Tasks[t].minimize()
      print('end minimization')

   for t in ['vetoTimeplanecalibrationX']:
      print('start',t)
      options.xCheck = True
      Tasks[t].Init(options,M)
      for n in range(options.nStart,options.nStart+options.nEvents):
        event = M.GetEvent(n)
        Tasks[t].ExecuteEvent(event)
        event.Reco_MuonTracks.Delete()
      Tasks[t].Finalize()
      print('end',t)

   mergeResiduals()
   mergeTimeDiffs()

if options.command=="TW":
   Tasks['vetoTimeWalk']       = vetoTimeCalibration.vetoTimeWalk()
   for t in Tasks:
      Tasks[t].Init(options,M)
      for n in range(options.nStart,options.nStart+options.nEvents):
        event = M.GetEvent(n)
        Tasks[t].ExecuteEvent(event)
        event.Reco_MuonTracks.Delete()
      Tasks[t].Finalize()

