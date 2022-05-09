#!/usr/bin/env python
import ROOT,os,sys,getopt,subprocess,atexit,time
import Monitor
import Scifi_monitoring
import Mufi_monitoring
import DAQ_monitoring
import EventDisplay_Task
import SndlhcMuonReco

def pyExit():
       print("Make suicide until solution found for freezing")
       os.system('kill '+str(os.getpid()))
atexit.register(pyExit)

from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("-A", "--auto", dest="auto", help="run in auto mode online monitoring",default=False,action='store_true')
parser.add_argument("--Nupdate", dest="Nupdate", help="frequence of updating online plots",default=100,type=int)
parser.add_argument("--Nlast",      dest="Nlast", help="last N events to analyze on file",default=10,type=int)
parser.add_argument("--sudo", dest="sudo", help="update files on EOS",default=False,action='store_true')

parser.add_argument("-M", "--online", dest="online", help="online mode",default=False,action='store_true')
parser.add_argument("--batch", dest="batch", help="batch mode",default=False,action='store_true')
parser.add_argument("--server", dest="server", help="xrootd server",default=os.environ["EOSSHIP"])
parser.add_argument("-r", "--runNumber", dest="runNumber", help="run number", type=int,default=-1)
parser.add_argument("-p", "--path", dest="path", help="run number",required=False,default="")
parser.add_argument("-P", "--partition", dest="partition", help="partition of data", type=int,required=False,default=-1)
parser.add_argument("-d", "--Debug", dest="debug", help="debug", default=False)
parser.add_argument("-cpp", "--convRawCPP", action='store_true', dest="FairTask_convRaw", help="convert raw data using ConvRawData FairTask", default=False)
parser.add_argument( "--withCalibration", action='store_true', dest="makeCalibration", help="make QDC and TDC calibration, not taking from raw data", default=False)

parser.add_argument("-f", "--inputFile", dest="fname", help="file name for MC", type=str,default=None,required=False)
parser.add_argument("-g", "--geoFile", dest="geoFile", help="geofile", required=True)
parser.add_argument("-b", "--heartBeat", dest="heartBeat", help="heart beat", default=10000,type=int)
parser.add_argument("-c", "--command", dest="command", help="command", default="")
parser.add_argument("-n", "--nEvents", dest="nEvents", help="number of events", default=-1,type=int)
parser.add_argument("-s", "--nStart", dest="nStart", help="first event", default=0,type=int)
parser.add_argument("-t", "--trackType", dest="trackType", help="DS or Scifi", default="DS")

parser.add_argument("--ScifiNbinsRes", dest="ScifiNbinsRes", default=100)
parser.add_argument("--Scifixmin", dest="Scifixmin", default=-2000.)
parser.add_argument("--ScifialignPar", dest="ScifialignPar", default=False)

parser.add_argument("--goodEvents", dest="goodEvents", action='store_true',default=False)
parser.add_argument("--withTrack", dest="withTrack", action='store_true',default=False)
parser.add_argument("--nTracks", dest="nTracks",default=0,type=int)
parser.add_argument("--save", dest="save", action='store_true',default=False)
parser.add_argument("--interactive", dest="interactive", action='store_true',default=False)

options = parser.parse_args()

options.dashboard = "currently_processed_file.txt"
if options.auto or options.batch: ROOT.gROOT.SetBatch(True)

def currentRun():
      with client.File() as f:
            f.open(options.server+options.path+options.dashboard)
            status, L = f.read()
            Lcrun = L.decode().split('\n')
            f.close()
      for l in Lcrun:
            if not l.find('FINISHED')<0:
               print("DAQ not running. Don't know which file to open.")
               print(Lcrun)
               curRun,curPart,start ="","",""
               break
            if not l.find('/home/snd/snd/') < 0:
                 tmp = l.split('/')
                 curRun = tmp[len(tmp)-2]
                 curPart = tmp[len(tmp)-1]
                 start = Lcrun[1]
                 break
      return curRun,curPart,start

if options.auto:
   options.online = True
   from XRootD import client
   from XRootD.client.flags import DirListFlags, OpenFlags, MkDirFlags, QueryCode
# search for current run
   if options.runNumber < 0:
        curRun = ""
        while curRun.find('run') < 0:
               curRun,curPart,options.startTime =  currentRun()
               if curRun.find('run') < 0:
                   print("sleep 300sec.",time.ctime())
                   time.sleep(300)
        options.runNumber = int(curRun.split('_')[1])
        options.partition = int(curPart.split('_')[1].split('.')[0])
else:
   if options.runNumber < 0:
       print("run number required for non-auto mode")
       os._exit(1)

# prepare tasks:
FairTasks = []
houghTransform = False # under construction, not yet tested
if houghTransform:
   muon_reco_task = SndlhcMuonReco.MuonReco()
   muon_reco_task.SetName("houghTransform")
   FairTasks.append(muon_reco_task)
else:
   import SndlhcTracking
   trackTask = SndlhcTracking.Tracking() 
   trackTask.SetName('simpleTracking')
   FairTasks.append(trackTask)

M = Monitor.Monitoring(options,FairTasks)
if options.nEvents < 0 : 
    if options.online: options.nEvents = M.converter.fiN.event.GetEntries()
    else:    options.nEvents = M.eventTree.GetEntries()

monitorTasks = {}
monitorTasks['daq']   = DAQ_monitoring.DAQ_boards()
monitorTasks['Scifi_hitMaps']   = Scifi_monitoring.Scifi_hitMaps()
monitorTasks['Mufi_hitMaps']   = Mufi_monitoring.Mufi_hitMaps()
monitorTasks['Mufi_QDCcorellations']   = Mufi_monitoring.Mufi_largeVSsmall()
monitorTasks['Scifi_residuals'] = Scifi_monitoring.Scifi_residuals()   # time consuming
if options.interactive:  monitorTasks['EventDisplay']   = EventDisplay_Task.twod()

for m in monitorTasks:
    monitorTasks[m].Init(options,M)

if not options.auto:   # default online/offline mode
   for n in range(options.nStart,options.nStart+options.nEvents):
     event = M.GetEvent(n)
     for m in monitorTasks:
        monitorTasks[m].ExecuteEvent(M.eventTree)

   if options.nEvents>0:
       for m in monitorTasks:
          monitorTasks[m].Plot()
   M.publishRootFile()
else: 
   """ auto mode
       check for open data file on the online machine
       analyze N events 
       re-open file and check for new events, 
       if none available, 
         check every 5 seconds: if no new file re-open again
         if new file, finish run, publish histograms, and restart with new file
   """
   N0 = 0
   lastFile = M.converter.fiN.GetName()
   tmp = lastFile.split('/')
   lastRun  = tmp[len(tmp)-2]
   lastPart = tmp[len(tmp)-1]
   nLast = options.nEvents
   nStart = nLast-options.Nlast
   M.updateHtml()
   M.updateSource(lastFile)
   while 1>0:
      for n in range(nStart,nLast):
        event = M.GetEvent(n)
        N0+=1
        # trackTask.event = event
        for m in monitorTasks:
           monitorTasks[m].ExecuteEvent(M.eventTree)
# update plots
        if N0%options.Nupdate==0:
           for m in monitorTasks:
               monitorTasks[m].Plot()
           if options.sudo: M.publishRootFile()

      M.updateSource(lastFile)
      newEntries = M.converter.fiN.event.GetEntries()
      if newEntries>nLast:
         nStart = max(nLast,newEntries-options.Nlast)
         nLast = newEntries
      else:  
      # check if file has changed
         curRun,curPart,options.startTime  =  currentRun()
         while curRun.find('run') < 0:
               curRun,curPart,options.startTime  =  currentRun()
               if curRun.find('run') < 0:  
                   print("sleep 300sec.",time.ctime())
                   time.sleep(300)
         if not curRun == lastRun:
            for m in monitorTasks:
               monitorTasks[m].Plot()
            print("run ",lastRun," has finished.")
            quit()  # reinitialize everything with new run number
         if not curPart == lastPart:
            lastPart = curPart
            lastFile = options.server+options.path+lastRun+"/"+ lastPart
            M.converter.fiN = ROOT.TFile.Open(lastFile)
         else:
            time.sleep(30) # sleep 30 seconds and check for new events
            print('DAQ inactive for 30sec. Last event = ',M.converter.fiN.event.GetEntries(), curRun,curPart,N0)
            nStart = nLast



