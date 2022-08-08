
#!/usr/bin/env python
import ROOT,os,sys,subprocess,atexit,time
from XRootD import client
from XRootD.client.flags import DirListFlags, OpenFlags, MkDirFlags, QueryCode
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
parser.add_argument("--ScifiResUnbiased", dest="ScifiResUnbiased", default=False)
parser.add_argument("--Mufixmin", dest="Mufixmin", default=-10.)

parser.add_argument("--goodEvents", dest="goodEvents", action='store_true',default=False)
parser.add_argument("--withTrack", dest="withTrack", action='store_true',default=False)
parser.add_argument("--nTracks", dest="nTracks",default=0,type=int)
parser.add_argument("--save", dest="save", action='store_true',default=False)
parser.add_argument("--interactive", dest="interactive", action='store_true',default=False)

options = parser.parse_args()
options.slowStream = True
options.startTime = ""
options.dashboard = "/mnt/raid1/data_online/currently_processed_file.txt"
if (options.auto and not options.interactive) or options.batch: ROOT.gROOT.SetBatch(True)

def currentRun():
      with client.File() as f:
            f.open(options.server+options.dashboard)
            status, L = f.read()
            Lcrun = L.decode().split('\n')
            f.close()
      curRun,curPart,start ="","",""
      for l in Lcrun:
            if not l.find('FINISHED')<0:
               print("DAQ not running. Don't know which file to open.")
               print(Lcrun)
               break
            if not l.find('.root') < 0:
                 tmp = l.split('/')
                 curRun = tmp[len(tmp)-2]
                 curPart = tmp[len(tmp)-1]
                 start = Lcrun[1]
                 break
      return curRun,curPart,start

if options.auto:
   options.online = True
# search for current run
   if options.runNumber < 0:
        curRun = ""
        while curRun.find('run') < 0:
               curRun,curPart,options.startTime =  currentRun()
               if curRun.find('run') < 0:
                   print("no run info, sleep 30sec.",time.ctime())
                   time.sleep(30)
        options.runNumber = int(curRun.split('_')[1])
        lastRun = curRun
        options.partition = 0   #   monitoring file set to data_0000.root   int(curPart.split('_')[1].split('.')[0])
else:
   if options.runNumber < 0:
       print("run number required for non-auto mode")
       os._exit(1)
# works only for runs on EOS
   if not options.server.find('eos')<0:
      runDir = "/eos/experiment/sndlhc/raw_data/commissioning/TI18/data/run_"+str(options.runNumber).zfill(6)
      jname = "run_timestamps.json"
      dirlist  = str( subprocess.check_output("xrdfs "+os.environ['EOSSHIP']+" ls "+runDir,shell=True) ) 
      if jname in dirlist:
         with client.File() as f:          
            f.open(options.server+runDir+"/run_timestamps.json")
            status, jsonStr = f.read()
         exec("date = "+jsonStr.decode())
         options.startTime = date['start_time']
         if 'stop_time' in date:
             options.startTime += " - "+ date['stop_time']
         options.startTime.replace('Z','')

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
if not options.fname:
   monitorTasks['daq']     = DAQ_monitoring.DAQ_boards()
   monitorTasks['rates']   = DAQ_monitoring.Time_evolution()
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
     if not options.online:
        if n%options.heartBeat == 0:
            print("--> run/event nr: %i %i %5.2F%%"%(M.eventTree.EventHeader.GetRunId(),n,n/options.nEvents*100))
# assume for the moment file does not contain fitted tracks
     for m in monitorTasks:
        monitorTasks[m].ExecuteEvent(M.eventTree)

   if options.nEvents>0:
       for m in monitorTasks:
          monitorTasks[m].Plot()
   M.publishRootFile()
   if options.sudo:
       print(options.runNumber,options.startTime)
       options.startTime += " #events="+str(options.nEvents)
       M.updateHtml()
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
   lastPart = 0   #   reading from second low speed DAQ stream    tmp[len(tmp)-1]
   nLast = options.nEvents
   nStart = nLast-options.Nlast
   if not options.interactive and options.sudo: M.updateHtml()
   lastFile = M.converter.fiN.GetName()
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
      newEntries = M.converter.fiN.Get('event').GetEntries()
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
               if not options.interactive:  monitorTasks[m].Plot()
            print("run ",lastRun," has finished.")
            quit()  # reinitialize everything with new run number
         if not curPart == lastPart and not options.slowStream:
            lastPart = curPart
            lastFile = options.server+options.path+lastRun+"/"+ lastPart
            M.converter.fiN = ROOT.TFile.Open(lastFile)
         else:
            time.sleep(10) # sleep 10 seconds and check for new events
            print('DAQ inactive for 10sec. Last event = ',M.converter.fiN.Get('event').GetEntries(), curRun,curPart,N0)
            nStart = nLast

