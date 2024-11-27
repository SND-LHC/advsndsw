import ROOT
import os,sys,subprocess,atexit
from array import array
import shipunit as u
import SndlhcMuonReco
import json
from decorators import *
from rootpyPickler import Unpickler
import time
from XRootD import client
import numpy as np
ROOT.gInterpreter.ProcessLine('#include "SiSensor.h"')
import rootUtils as ut

from datetime import datetime
from pathlib import Path

def pyExit():
       "unfortunately need as bypassing an issue related to use xrootd"
       os.system('kill '+str(os.getpid()))
atexit.register(pyExit)


A,B = ROOT.TVector3(),ROOT.TVector3()
freq      =  160.316E6

h={}
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("-r", "--runNumber", dest="runNumber", help="run number", type=int,required=False)
parser.add_argument("-p", "--path", dest="path", help="path to data file",required=False,default=os.environ["EOSSHIP"]+"/eos/experiment/sndlhc/convertedData/physics/2022/")
parser.add_argument("-praw", "--pathRaw", dest="pathRaw", help="path to raw data file",required=False,default="/eos/experiment/sndlhc/raw_data/physics/2022/")
parser.add_argument("-f", "--inputFile", dest="inputFile", help="input file data and MC",default="",required=False)
parser.add_argument("-g", "--geoFile", dest="geoFile", help="geofile", default=os.environ["EOSSHIP"]+"/eos/experiment/sndlhc/convertedData/physics/2022/geofile_sndlhc_TI18_V0_2022.root")
parser.add_argument("-P", "--partition", dest="partition", help="partition of data", type=int,required=False,default=-1)
parser.add_argument("--server", dest="server", help="xrootd server",default=os.environ["EOSSHIP"])
parser.add_argument("-X", dest="extraInfo", help="print extra event info",default=True)

parser.add_argument("-par", "--parFile", dest="parFile", help="parameter file", default=os.environ['ADVSNDSW_ROOT']+"/python/TrackingParams.xml")
parser.add_argument("-hf", "--HoughSpaceFormat", dest="HspaceFormat", help="Hough space representation. Should match the 'Hough_space_format' name in parFile, use quotes", default='linearSlopeIntercept')

options = parser.parse_args()
options.storePic = ''
trans2local = False
runInfo = False
try:
   fg  = ROOT.TFile.Open(options.server+options.p+"RunInfodict.root")
   pkl = Unpickler(fg)
   runInfo = pkl.load('runInfo')
   fg.Close()
except: pass

import SndlhcGeo
geo = SndlhcGeo.GeoInterface(options.geoFile)

detSize = {}
# same sensors used in both the AdvTarget and the AdvMuFilter
detSize[0] =[ROOT.advsnd.sensor_width/ROOT.advsnd.strips, 
             ROOT.advsnd.sensor_width/ROOT.advsnd.strips, 
             geo.snd_geo.AdvTarget.TTZ]# check the channel sizes
detSize[1] = detSize[0]

mc = False


# Initialize FairLogger: set severity and verbosity
logger = ROOT.FairLogger.GetLogger()
logger.SetColoredLog(True)
logger.SetLogVerbosityLevel('low')
logger.SetLogScreenLevel('WARNING')
logger.SetLogToScreen(True)

run      = ROOT.FairRunAna()
ioman = ROOT.FairRootManager.Instance()

if options.inputFile=="":
  f=ROOT.TFile.Open(options.path+'sndsw_raw_'+str(options.runNumber).zfill(6)+'.root')
else:
  f=ROOT.TFile.Open(options.path+options.inputFile)

if f.FindKey('cbmsim'):
        eventTree = f.cbmsim
        runId = 'sim'
        if eventTree.GetBranch('AdvTargetPoint'): mc = True
else:   
        eventTree = f.rawConv
        ioman.SetTreeName('rawConv')

outFile = ROOT.TMemFile('dummy','CREATE')
source = ROOT.FairFileSource(f)
run.SetSource(source)
sink = ROOT.FairRootFileSink(outFile)
run.SetSink(sink)

HT_tasks = {'muon_reco_task_Sf':SndlhcMuonReco.MuonReco(),
            'muon_reco_task_nuInt':SndlhcMuonReco.MuonReco()}
for ht_task in HT_tasks.values():
    run.AddTask(ht_task)

#avoiding some error messages
xrdb = ROOT.FairRuntimeDb.instance()
xrdb.getContainer("FairBaseParSet").setStatic()
xrdb.getContainer("FairGeoParSet").setStatic()

for ht_task in HT_tasks.values():
    ht_task.SetParFile(options.parFile)
    ht_task.SetHoughSpaceFormat(options.HspaceFormat)
    # force the output of reco task to genfit::Track
    # as the display code looks for such output
    ht_task.ForceGenfitTrackFormat()
HT_tasks['muon_reco_task_Sf'].SetTrackingCase('passing_mu_AdvTarget')
HT_tasks['muon_reco_task_nuInt'].SetTrackingCase('nu_interaction_products')

run.Init()
OT = sink.GetOutTree()
eventTree = ioman.GetInTree()
eventTree.GetEvent(0)
#if eventTree.EventHeader.ClassName() == 'SNDLHCEventHeader':
#   geo.modules['Scifi'].InitEvent(eventTree.EventHeader)
#   geo.modules['MuFilter'].InitEvent(eventTree.EventHeader)
## if faireventheader, rely on user to select correct geofile.

nav = ROOT.gGeoManager.GetCurrentNavigator()

# get filling scheme
try:
           runNumber = eventTree.EventHeader.GetRunId()
           fg  = ROOT.TFile.Open(os.environ['EOSSHIP']+options.p+'FSdict.root')
           pkl = Unpickler(fg)
           FSdict = pkl.load('FSdict')
           fg.Close()
           if runNumber in FSdict: fsdict = FSdict[runNumber]
           else:  fsdict = False
except:
           print('continue without knowing filling scheme')
           fsdict = False

startTimeOfRun = {}
def getStartTime(runNumber):
      if runNumber in startTimeOfRun : return startTimeOfRun[runNumber]
      runDir = options.pathRaw+"run_"+str(runNumber).zfill(6)
      jname = "run_timestamps.json"
      dirlist  = str( subprocess.check_output("xrdfs "+options.server+" ls "+runDir,shell=True) ) 
      if not jname in dirlist: return False
      with client.File() as f:
               f.open(options.server+runDir+"/run_timestamps.json")
               status, jsonStr = f.read()
               f.close()
      date = json.loads(jsonStr)
      time_str = date['start_time'].replace('Z','')
      time_obj = time.strptime(time_str, '%Y-%m-%dT%H:%M:%S')
      startTimeOfRun[runNumber] = time.mktime(time_obj)
      return startTimeOfRun[runNumber]


Nlimit = 4
onlyScifi = False

def goodEvent(event):
# can be replaced by any user selection
           # not tuned for AdvSND case
           return True
           stations = {'Scifi':{},'Mufi':{}}
           if event.Digi_AdvTargetHits.GetEntries()>25: return False
           for d in event.Digi_AdvTargetHits:
               stations['Scifi'][d.GetDetectorID()//1000000] = 1
           for d in event.Digi_AdvMuFilterHits:
               plane = d.GetDetectorID()//1000
               stations['Mufi'][plane] = 1
           totalN = len(stations['Mufi'])+len(stations['Scifi'])
           if len(stations['Scifi'])>4 and len(stations['Mufi'])>6: return True
           else: False
           if onlyScifi and len(stations['Scifi'])>Nlimit: return True
           elif not onlyScifi  and totalN >  Nlimit: return True
           else: return False
def bunchXtype():
# check for b1,b2,IP1,IP2
        xing = {'all':True,'B1only':False,'B2noB1':False,'noBeam':False}
        if fsdict:
             T   = eventTree.EventHeader.GetEventTime()
             bunchNumber = int(T%(4*3564)/4+0.5)
             nb1 = (3564 + bunchNumber - fsdict['phaseShift1'])%3564
             nb2 = (3564 + bunchNumber - fsdict['phaseShift1']- fsdict['phaseShift2'])%3564
             b1 = nb1 in fsdict['B1']
             b2 = nb2 in fsdict['B2']
             IP1 = False
             IP2 = False
             if b1:
                IP1 =  fsdict['B1'][nb1]['IP1']
             if b2:
                IP2 =  fsdict['B2'][nb2]['IP2']
             if b2 and not b1:
                xing['B2noB1'] = True
             if b1 and not b2 and not IP1:
                xing['B1only'] = True
             if not b1 and not b2: xing['noBeam'] = True
        return xing

def drawAdvHits(g, colour):
       """Takes TGraph g and draws the graphs markers with the TColor given in list colour."""
       n = g.GetN()
       # Draw highest density points last
       sorted_indices = np.argsort(colour)
       for i_unsorted in range(0, n):
              i = int(sorted_indices[i_unsorted])

              x = g.GetPointX(i)
              y = g.GetPointY(i)

              h["markerCollection"].append(ROOT.TEllipse(x, y, 0.8,0.8))
              h["markerCollection"][-1].SetLineWidth(0)
              h["markerCollection"][-1].SetFillColor(colour[i])
              h["markerCollection"][-1].Draw("SAME")


def drawLegend(max_density, max_QDC, n_legend_points):
       """Draws legend for hit colour"""
       h['simpleDisplay'].cd(1)
       n_legend_points = 5
       padLegTarget = ROOT.TPad("legend","legend",0.4,0.15,0.4+0.2, 0.15+0.2)
       padLegTarget.SetFillStyle(4000)
       padLegTarget.Draw()
       padLegTarget.cd()
       text_target_legend = ROOT.TLatex()
       text_target_legend.SetTextAlign(11)
       text_target_legend.SetTextFont(42)
       text_target_legend.SetTextSize(.15)
       for i in range(n_legend_points) :
              if i < (n_legend_points - 1) :
                     text_target_legend.DrawLatex((i+0.3)*(1./(n_legend_points+2)), 0.2, "{:d}".format(int(i*max_density/(n_legend_points-1))))
                     #text_target_legend.DrawLatex((i+0.3)*(1./(n_legend_points+2)), 0.,  "{:.0f}".format(int(i*max_QDC/(n_legend_points-1))))
              else :
                     text_target_legend.DrawLatex((i+0.3)*(1./(n_legend_points+2)), 0.2, "{:d} hits/cm".format(int(i*max_density/(n_legend_points-1))))
                     #text_target_legend.DrawLatex((i+0.3)*(1./(n_legend_points+2)), 0.,  "{:.0f} QDC units".format(int(i*max_QDC/(n_legend_points-1))))
                     
              h["markerCollection"].append(ROOT.TEllipse((i+0.15)*(1./(n_legend_points+2)), 0.26, 0.05/4, 0.05))
              h["markerCollection"][-1].SetFillColor(ROOT.TColor.GetPalette()[int(float(i*max_density/(n_legend_points-1))/max_density*(len(ROOT.TColor.GetPalette())-1))])
              h["markerCollection"][-1].Draw("SAME")
              
              # in case we want to use that in the future
              #h["markerCollection"].append(ROOT.TBox((i+0.15)*(1./(n_legend_points+2))-0.05/4 , 0.06 - 0.05, (i+0.15)*(1./(n_legend_points+2))+0.05/4, 0.06 + 0.05))
              #h["markerCollection"][-1].SetFillColor(ROOT.TColor.GetPalette()[int(float(i*max_QDC/(n_legend_points-1))/max_QDC*(len(ROOT.TColor.GetPalette())-1))])
              #h["markerCollection"][-1].Draw("SAME")

def getAdvHitDensity(g, x_range=0.5):
       """Takes ROOT TGraph g and returns array with number of hits within x_range cm of each hit."""
       ret = []
       for i in range(g.GetN()):
              x_i = g.GetPointX(i)
              y_i = g.GetPointY(i)
              density = 0
              for j in range(g.GetN()):
                     x_j = g.GetPointX(j)
                     y_j = g.GetPointY(j)
                     if ((x_i - x_j)**2 + (y_i - y_j)**2) <= x_range**2:
                            density += 1
              ret.append(density)
       return ret

#original options of loopEvents for future reference
#def loopEvents(start=0,save=False,goodEvents=False,withTrack=-1,withHoughTrack=-1,nTracks=0,minSipmMult=1,withTiming=False, option=None,Setup='',verbose=0,auto=False):
def loopEvents(start=0,save=False,withHoughTrack=-1,nTracks=0,option=None,Setup='',verbose=0,auto=False,hitColour=None):
 if 'simpleDisplay' not in h: ut.bookCanvas(h,key='simpleDisplay',title='simple event display',nx=1200,ny=1600,cx=1,cy=2)
 h['simpleDisplay'].cd(1)
 zStart = -300. # AdvSND MC
 if Setup == 'H6': zStart = 60.
 if Setup == 'TP': zStart = -50. # old coordinate system with origin in middle of target
 if 'xz' in h: 
    h.pop('xz').Delete()
    h.pop('yz').Delete()
 else:
    h['xmin'],h['xmax'] = -100.,100. # AdvSND MC
    h['ymin'],h['ymax'] = -100.,100. # AdvSND MC
    h['zmin'],h['zmax'] = zStart,zStart+500.
    for d in ['xmin','xmax','ymin','ymax','zmin','zmax']: h['c'+d]=h[d]
 ut.bookHist(h,'xz','; z [cm]; x [cm]',500,h['czmin'],h['czmax'],200,h['cxmin'],h['cxmax'])
 ut.bookHist(h,'yz','; z [cm]; y [cm]',500,h['czmin'],h['czmax'],200,h['cymin'],h['cymax'])

 proj = {1:'xz',2:'yz'}
 h['xz'].SetStats(0)
 h['yz'].SetStats(0)

 N = -1
 Tprev = -1
 A,B = ROOT.TVector3(),ROOT.TVector3()
 ptext={0:'   Y projection',1:'   X projection'}
 text = ROOT.TLatex()
 event = eventTree
 OT = sink.GetOutTree()
 if withHoughTrack==0: OT = eventTree
 if type(start) == type(1):
    s = start
    e = event.GetEntries()
 else:
    s = 0
    e = len(start)
 for N in range(s,e):
    if type(start) == type(1): rc = event.GetEvent(N)
    else: rc = event.GetEvent(start[N])
    #if goodEvents and not goodEvent(event): continue
    nHoughtracks = 0
    OT.Reco_MuonTracks = ROOT.TObjArray(10)
    if withHoughTrack > 0:
       rc = source.GetInTree().GetEvent(N)
       # Delete SndlhcMuonReco kalman tracks container
       for ht_task in HT_tasks.values():
           ht_task.kalman_tracks.Delete()
       if withHoughTrack==1:
            HT_tasks['muon_reco_task_Sf'].Exec(0)
       elif withHoughTrack==2:
            HT_tasks['muon_reco_task_Sf'].Exec(0)
       elif withHoughTrack==4:
            HT_tasks['muon_reco_task_nuInt'].Exec(0)
       # Save the tracks in OT.Reco_MuonTracks object
       for ht_task in HT_tasks.values():
           for trk in ht_task.kalman_tracks:
               OT.Reco_MuonTracks.Add(trk)
       nHoughtracks = OT.Reco_MuonTracks.GetEntries()
       if nHoughtracks>0: print('number of tracks by HT:', nHoughtracks)
       if nHoughtracks<nTracks: continue

    if verbose>0:
       for aTrack in OT.Reco_MuonTracks:
           print(aTrack.__repr__())
           mom    = aTrack.getFittedState().getMom()
           pos      = aTrack.getFittedState().getPos()
           mom.Print()
           pos.Print()
    T,dT = 0,0
    T = event.EventHeader.GetEventTime()
    runId = eventTree.EventHeader.GetRunId()
    if Tprev >0: dT = T-Tprev
    Tprev = T

    digis = []
    if event.FindBranch("Digi_AdvTargetHits"): digis.append(event.Digi_AdvTargetHits)
    if event.FindBranch("Digi_AdvMuFilterHits"): digis.append(event.Digi_AdvMuFilterHits)
    empty = True
    for x in digis:
       if x.GetEntries()>0:
         if empty: print( "event -> %i"%N)
         empty = False
    if empty: continue
    h['hitCollectionX']= {'AdvTarget':[0,ROOT.TGraphErrors()],'AdvMuFilter':[0,ROOT.TGraphErrors()]}
    h['hitCollectionY']= {'AdvTarget':[0,ROOT.TGraphErrors()],'AdvMuFilter':[0,ROOT.TGraphErrors()]}
    if hitColour:
           h['hitColourX'] = {'AdvTarget' : [], 'AdvMuFilter' : []}
           h['hitColourY'] = {'AdvTarget' : [], 'AdvMuFilter' : []}
           h["markerCollection"] = []
    h['firedChannelsX']= {'AdvTarget':[0,0,0],'AdvMuFilter':[0,0,0]}
    h['firedChannelsY']= {'AdvTarget':[0,0,0],'AdvMuFilter':[0,0,0,0]}
    systems = {0:'AdvTarget', 1:'AdvMuFilter'}
    for collection in ['hitCollectionX','hitCollectionY']:
       for c in h[collection]:
          rc=h[collection][c][1].SetName(c)
          rc=h[collection][c][1].Set(0)

    if hitColour:
       h["markerCollection"] = []

    dTs = "%5.2Fns"%(dT/freq*1E9)
    # find detector which triggered
    #minT = firstTimeStamp(event)
    #if minT[0] < 1000000000:
    #    dTs+= "    " + str(minT[1].GetDetectorID())
    for p in proj:
       rc = h[ 'simpleDisplay'].cd(p)
       h[proj[p]].Draw('b')

    drawDetectors()
 
    for D in digis:
      for digi in D:
         detID = digi.GetDetectorID()
         #sipmMult = 1
         if digi.GetName()  == 'AdvMuFilterHit':
            system = 1
            geo.modules['AdvMuFilter'].GetPosition(detID,A,B)
            #sipmMult = len(digi.GetAllSignals(False,False))
            #if sipmMult<minSipmMult and (system==1 or system==2): continue
         else:
            geo.modules['AdvTarget'].GetPosition(detID,A,B)
            system = 0
         curPath = nav.GetPath()
         tmp = curPath.rfind('/')
         nav.cd(curPath[:tmp])
         globA,locA = array('d',[A[0],A[1],A[2]]),array('d',[A[0],A[1],A[2]])
         if trans2local:   nav.MasterToLocal(globA,locA)
         Z = A[2]
         if digi.isVertical():
                   collection = 'hitCollectionX'
                   Y = locA[0]
                   sY = detSize[system][0]
         else:                         
                   collection = 'hitCollectionY'
                   Y = locA[1]
                   sY = detSize[system][1]
         c = h[collection][systems[system]]
         rc = c[1].SetPoint(c[0],Z, Y)
         rc = c[1].SetPointError(c[0],detSize[system][2],sY)
         c[0]+=1 

         #fillNode(curPath)

         if digi.isVertical():  F = 'firedChannelsX'
         else:                     F = 'firedChannelsY'
         '''
         ns = max(1,digi.GetnSides())
         for side in range(ns):
             for m in  range(digi.GetnSiPMs()):
                   qdc = digi.GetSignal(m+side*digi.GetnSiPMs())
                   if qdc < 0 and qdc > -900:  h[F][systems[system]][1]+=1
                   elif not qdc<0:   
                       h[F][systems[system]][0]+=1
                       if len(h[F][systems[system]]) < 2+side: continue
                       h[F][systems[system]][2+side]+=qdc
         '''

    if hitColour == "q" :
       for orientation in ['X', 'Y']:
              max_density = 100              
              max_QDC = 200 # not used for now
              density_AdvTarget = np.clip(0, max_density, getAdvHitDensity(h['hitCollection'+orientation]['AdvTarget'][1]))
              density_AdvMuFilter = np.clip(0, max_density, getAdvHitDensity(h['hitCollection'+orientation]['AdvMuFilter'][1]))
              for det_name in ['AdvTarget', 'AdvMuFilter']:
                if det_name =='AdvTarget': density = density_AdvTarget
                else: density= density_AdvMuFilter
                for i in range(h['hitCollection'+orientation][det_name][1].GetN()) :
                       h['hitColour'+orientation][det_name].append(ROOT.TColor.GetPalette()[int(float(density[i])/max_density*(len(ROOT.TColor.GetPalette())-1))])

       drawLegend(max_density, max_QDC, 5)
    
    k = 1
    moreEventInfo = []
    for collection in ['hitCollectionX','hitCollectionY']:
       h['simpleDisplay'].cd(k)
       drawInfo(h['simpleDisplay'], k, runId, N, T)
       k+=1
       for c in h[collection]:
          F = collection.replace('hitCollection','firedChannels')
          pj = collection.split('ion')[1]
          #if pj =="X" or c=="Scifi":
          #    atext = "%1s %5s %3i  +:%3i -:%3i qdc :%5.1F"%(pj,c,h[collection][c][1].GetN(),h[F][c][0],h[F][c][1],h[F][c][2])
          #else:
          #    atext = "%1s %5s %3i  +:%3i -:%3i qdcL:%5.1F qdcR:%5.1F"%(pj,c,h[collection][c][1].GetN(),h[F][c][0],h[F][c][1],h[F][c][2],h[F][c][3])
          #moreEventInfo.append(atext)
          #print(atext)
          if h[collection][c][1].GetN()<1: continue
          if hitColour not in ["q"] :
              h[collection][c][1].SetMarkerStyle(20)
              h[collection][c][1].SetMarkerSize(0.5)
              if c=='AdvTarget':
                h[collection][c][1].SetMarkerColor(ROOT.kBlue+2)
              else: 
                h[collection][c][1].SetMarkerColor(ROOT.kRed+2)
              rc=h[collection][c][1].Draw('sameP')
              h['display:'+c]=h[collection][c][1]
          elif hitColour == "q" :
              drawAdvHits(h[collection][c][1], h['hitColour'+collection[-1]][c])

    T0 = eventTree.EventHeader.GetEventTime()
    if type(start) == type(1): rc = event.GetEvent(N-1)
    else: rc = event.GetEvent(start[N]-1)
    delTM1 = eventTree.EventHeader.GetEventTime() - T0
    if type(start) == type(1): rc = event.GetEvent(N+1)
    else: rc = event.GetEvent(start[N]+1)
    delTP1 = eventTree.EventHeader.GetEventTime() - T0
    atext = "timing info, prev event: %6i cc  next event: %6i cc"%(delTM1,delTP1)
    #moreEventInfo.append(atext)
    if type(start) == type(1): rc = event.GetEvent(N)
    else: rc = event.GetEvent(start[N])

    k = 1
    for collection in ['hitCollectionX','hitCollectionY']:
       h['simpleDisplay'].cd(k)
       drawInfo(h['simpleDisplay'], k, runId, N, T,moreEventInfo)
       k+=1
            
    h['simpleDisplay'].Update()
    #if withTiming: timingOfEvent()
    addTrack(OT)

    if option == "2tracks": 
          rc = twoTrackEvent(sMin=10,dClMin=7,minDistance=0.5,sepDistance=0.5)
          if not rc: rc = twoTrackEvent(sMin=10,dClMin=7,minDistance=0.5,sepDistance=0.75)
          if not rc: rc = twoTrackEvent(sMin=10,dClMin=7,minDistance=0.5,sepDistance=1.0)
          if not rc: rc = twoTrackEvent(sMin=10,dClMin=7,minDistance=0.5,sepDistance=1.75)
          if not rc: rc = twoTrackEvent(sMin=10,dClMin=7,minDistance=0.5,sepDistance=2.5)
          if not rc: rc = twoTrackEvent(sMin=10,dClMin=7,minDistance=0.5,sepDistance=3.0)

    if verbose>0: dumpChannels()
    if save: h['simpleDisplay'].Print('{:0>2d}-event_{:04d}'.format(runId,N)+'.png')
    if auto:
        h['simpleDisplay'].Print(options.storePic+str(runId)+'-event_'+str(event.EventHeader.GetEventNumber())+'.png')
    if not auto:
       rc = input("hit return for next event or p for print or q for quit: ")
       if rc=='p': 
             h['simpleDisplay'].Print(options.storePic+str(runId)+'-event_'+str(event.EventHeader.GetEventNumber())+'.png')
       if rc=='q': break
 if save: os.system("convert -delay 60 -loop 0 event*.png animated.gif")

def addTrack(OT,scifi=False):
   xax = h['xz'].GetXaxis()
   nTrack = 0
   for   aTrack in OT.Reco_MuonTracks:
      trackColor = ROOT.kRed

      if aTrack.GetUniqueID()==11: trackColor = ROOT.kAzure-2 # HT AdvTarget track
      # HT cross-system track fit
      if aTrack.GetUniqueID()==15: trackColor = ROOT.kOrange+7
      S = aTrack.getFitStatus()
      if not S.isFitConverged():
         print('Kalman Filter track fit did not converge!')
         continue
      for p in [0,1]:
          h['aLine'+str(nTrack*10+p)] = ROOT.TGraph()

      zEx = xax.GetBinCenter(1)
      mom    = aTrack.getFittedState().getMom()
      pos      = aTrack.getFittedState().getPos()
      lam      = (zEx-pos.z())/mom.z()
      Ex        = [pos.x()+lam*mom.x(),pos.y()+lam*mom.y()]
      for p in [0,1]:   h['aLine'+str(nTrack*10+p)].SetPoint(0,zEx,Ex[p])

      for i in range(aTrack.getNumPointsWithMeasurement()):
         state = aTrack.getFittedState(i)
         pos    = state.getPos()
         for p in [0,1]:
             h['aLine'+str(nTrack*10+p)].SetPoint(i+1,pos[2],pos[p])

      zEx = xax.GetBinCenter(xax.GetLast())
      mom    = aTrack.getFittedState().getMom()
      pos      = aTrack.getFittedState().getPos()
      lam      = (zEx-pos.z())/mom.z()
      Ex        = [pos.x()+lam*mom.x(),pos.y()+lam*mom.y()]
      for p in [0,1]:   h['aLine'+str(nTrack*10+p)].SetPoint(i+2,zEx,Ex[p])

      for p in [0,1]:
             tc = h[ 'simpleDisplay'].cd(p+1)
             h['aLine'+str(nTrack*10+p)].SetLineColor(trackColor)
             h['aLine'+str(nTrack*10+p)].SetLineWidth(2)
             h['aLine'+str(nTrack*10+p)].Draw('same')
             tc.Update()
             h[ 'simpleDisplay'].Update()
      nTrack+=1

def twoTrackEvent(sMin=10,dClMin=7,minDistance=1.5,sepDistance=0.5):
        trackTask.clusScifi.Clear()
        trackTask.scifiCluster()
        clusters = trackTask.clusScifi
        sortedClusters={}
        for aCl in clusters:
           so = aCl.GetFirst()//100000
           if not so in sortedClusters: sortedClusters[so]=[]
           sortedClusters[so].append(aCl)
        if len(sortedClusters)<sMin: return False
        M=0
        for x in sortedClusters:
           if len(sortedClusters[x]) == 2:  M+=1
        if M < dClMin: return False
        seeds = {}
        S = [-1,-1]
        for o in range(0,2):
# same procedure for both projections
# take seeds from from first station with 2 clusters
             for s in range(1,6):
                 x = 10*s+o
                 if x in sortedClusters:
                    if len(sortedClusters[x])==2:
                       sortedClusters[x][0].GetPosition(A,B)
                       if o%2==1: pos0 = (A[0]+B[0])/2
                       else: pos0 = (A[1]+B[1])/2
                       sortedClusters[x][1].GetPosition(A,B)
                       if o%2==1: pos1 = (A[0]+B[0])/2
                       else: pos1 = (A[1]+B[1])/2
                       if abs(pos0-pos1) > minDistance:
                         S[o] = s
                         break
             if S[o]<0: break  # no seed found
             seeds[o]={}
             k = -1
             for c in sortedClusters[S[o]*10+o]:
                 k += 1
                 c.GetPosition(A,B)
                 if o%2==1: pos = (A[0]+B[0])/2
                 else: pos = (A[1]+B[1])/2
                 seeds[o][k] = [[c,pos]]
             if k!=1: continue
             if abs(seeds[o][0][0][1] - seeds[o][1][0][1]) < sepDistance: continue
             for s in range(1,6):
               if s==S[o]: continue
               for c in sortedClusters[s*10+o]:
                   c.GetPosition(A,B)
                   if o%2==1: pos = (A[0]+B[0])/2
                   else: pos = (A[1]+B[1])/2
                   for k in range(2):
                        if  abs(seeds[o][k][0][1] - pos) < sepDistance:
                           seeds[o][k].append([c,pos])
        if S[0]<0 or S[1]<0:
            passed = False
        else:
           passed = True
           for o in range(0,2):
              for k in range(2):
                  if len(seeds[o][k])<3:
                      passed = False
                      break
        print(passed)
        if passed:
           tracks = []
           for k in range(2):
             # arbitrarly combine X and Y of combination 0
               n = 0
               hitlist = {}
               for o in range(0,2):
                   for X in seeds[o][k]:
                      hitlist[n] = X[0]
                      n+=1
               theTrack = trackTask.fitTrack(hitlist)
               if not hasattr(theTrack,"getFittedState"):
                    validTrack = False
                    continue
               fitStatus = theTrack.getFitStatus()
               if not fitStatus.isFitConverged():
                    theTrack.Delete()
               else: 
                    tracks.append(theTrack)
           if len(tracks)==2:
                 OT = sink.GetOutTree()
                 OT.Reco_MuonTracks = tracks
                 addTrack(OT,True) 
        return passed

def drawDetectors():
    """ Read the detector geometry and draw its elements on the display. """
    nodes = {'volAdvTarget_1/volTargetWall_1':ROOT.kGray+1}
    geoVer = checkGeoVersion()
    if geoVer!=2:
        print('This is an older version of the geometry and it is likely incompatible with this version of the drawDetectors function of the 2dED.')
    for i in range(geo.snd_geo.AdvTarget.nTT):
        nodes['volAdvTarget_1/volTargetWall_{}'.format(i)]=ROOT.kGray+1
        nodes['volAdvTarget_1/TrackingStation_{}'.format(i)]=ROOT.kBlue

    for i in range(geo.snd_geo.AdvMuFilter.Nplanes):
        nodes['volAdvMuFilter_0/volFeWall_{}'.format(i)] = ROOT.kGreen -6
        nodes['volAdvMuFilter_0/TrackingStation_{}'.format(i)] = ROOT.kRed -2
        if i > geo.snd_geo.AdvMuFilter.Nplanes-2: continue
        for j in range(geo.snd_geo.AdvMuFilter.Nplanes-2):
            if geoVer == 1:
                nodes['volAdvMuFilter_0/volMuonSysDet_{}_{}/volBar_{}'.format(i, i, i*100+j)] = ROOT.kGray
            elif geoVer == 2:
                nodes['volAdvMuFilter_0/volMuonSysDet_{}'.format(i)] = ROOT.kGray
    
    for i in range(2):
        nodes['volAdvMuFilter_0/volVertCoil_{}'.format(i)] = ROOT.kOrange+1
        nodes['volAdvMuFilter_0/volCoil_{}'.format(i)] = ROOT.kOrange+1
    if geoVer != 2:
        nodes['volAdvMuFilter_0/volMagTracker1_10000'] = ROOT.kGray   
        nodes['volAdvMuFilter_0/volMagTracker2_10001'] = ROOT.kGray  
        nodes['volAdvMuFilter_0/volDownMagTracker_10002'] = ROOT.kGray
        nodes['volAdvMuFilter_0/volDownstreamMagnet_0/volDownVertCoil1_0'] = ROOT.kOrange+1
        nodes['volAdvMuFilter_0/volDownstreamMagnet_0/volDownVertCoil2_0'] = ROOT.kOrange+1
        nodes['volAdvMuFilter_0/volDownstreamMagnet_0/volIronYoke_0'] = ROOT.kGreen-2
        nodes['volAdvMuFilter_0/volDownstreamMagnet_0/volDownCoil_20000'] = ROOT.kOrange+1
        nodes['volAdvMuFilter_0/volDownstreamMagnet_0/volIronCore_0'] = ROOT.kGreen+3
        IC_lines = GetIronCoreLines()
    nodes['volAdvMuFilter_0/volMagTracker_10001'] = ROOT.kGray
    nodes['volAdvMuFilter_0/volMagTracker_10002'] = ROOT.kGray

    #passNodes = {'Wall','Coil', 'Block','Yoke'}
    passNodes = {'Wall','Coil'}
    proj = {'X':0, 'Y':1}
    for node_ in nodes:
        node = '/cave_1/Detector_0/'+node_
        for p in proj:
            if node in h and any(passNode in node for passNode in passNodes):
                X = h[node+p]
                c = proj[p]
                h['simpleDisplay'].cd(c+1)
                X.Draw('f&&same')
                X.Draw('same')
            else:
                if not nav.CheckPath(node): continue
                nav.cd(node)
                N = nav.GetCurrentNode()
                # get rid of the Coil element for the vertical projection
                # FIXME change the naming of the Coil and VertCoil in the geo file
                if (N.GetVolume().GetName().find('Coil')>0
                   and N.GetVolume().GetName().find('VertCoil')<0
                   and p=='X'): continue                
                S = N.GetVolume().GetShape()
                scale = 1
                if N.GetVolume().GetName().find('FeWall')>0:
                   scale = 1 #2
                dx,dy,dz = S.GetDX()/scale,S.GetDY()/scale,S.GetDZ()
                ox,oy,oz = S.GetOrigin()[0],S.GetOrigin()[1],S.GetOrigin()[2]
                P = {}
                M = {}                
                if p=='X':
                    P['LeftBottom'] = array('d',[-dx+ox,oy,-dz+oz])
                    P['LeftTop'] = array('d',[dx+ox,oy,-dz+oz])
                    P['RightBottom'] = array('d',[-dx+ox,oy,dz+oz])
                    P['RightTop'] = array('d',[dx+ox,oy,dz+oz])
                elif p=='Y':
                    P['LeftBottom'] = array('d',[ox,-dy+oy,-dz+oz])
                    P['LeftTop'] = array('d',[ox,dy+oy,-dz+oz])
                    P['RightBottom'] = array('d',[ox,-dy+oy,dz+oz])
                    P['RightTop'] = array('d',[ox,dy+oy,dz+oz])
                else: continue
                for C in P:
                    M[C] = array('d',[0,0,0])
                    nav.LocalToMaster(P[C],M[C])
                h[node+p] = ROOT.TPolyLine()
                X = h[node+p]
                c = proj[p]
                if node == '/cave_1/Detector_0/volAdvMuFilter_0/volDownstreamMagnet_0/volIronCore_0':
                    X.SetPoint(0,M['LeftBottom'][2],IC_lines[p]['LeftBottom'][c])
                    X.SetPoint(1,M['LeftTop'][2],IC_lines[p]['LeftTop'][c])
                    X.SetPoint(2,M['RightTop'][2],IC_lines[p]['RightTop'][c])
                    X.SetPoint(3,M['RightBottom'][2],IC_lines[p]['RightBottom'][c])
                    X.SetPoint(4,M['LeftBottom'][2],IC_lines[p]['LeftBottom'][c])
                elif node == '/cave_1/Detector_0/volAdvMuFilter_0/volDownstreamMagnet_0/volDownCoil_20000' and p == 'Y':
                    X.SetPoint(0,M['LeftBottom'][2],IC_lines[p]['LeftBottom'][c]-7.)
                    X.SetPoint(1,M['LeftTop'][2],IC_lines[p]['LeftTop'][c]+7.)
                    X.SetPoint(2,M['RightTop'][2],IC_lines[p]['RightTop'][c]+7.)
                    X.SetPoint(3,M['RightBottom'][2],IC_lines[p]['RightBottom'][c]-7.)
                    X.SetPoint(4,M['LeftBottom'][2],IC_lines[p]['LeftBottom'][c]-7.)
                else:
                    X.SetPoint(0,M['LeftBottom'][2],M['LeftBottom'][c])
                    X.SetPoint(1,M['LeftTop'][2],M['LeftTop'][c])
                    X.SetPoint(2,M['RightTop'][2],M['RightTop'][c])
                    X.SetPoint(3,M['RightBottom'][2],M['RightBottom'][c])
                    X.SetPoint(4,M['LeftBottom'][2],M['LeftBottom'][c])
                X.SetLineColor(nodes[node_])
                X.SetLineWidth(1)
                h['simpleDisplay'].cd(c+1)
                if any(passNode in node for passNode in passNodes):
                   X.SetFillColorAlpha(nodes[node_], 0.5)
                   X.Draw('f&&same')
                X.Draw('same')

def checkGeoVersion():
    """ Check the geometry version of the AdvSND  """
    vol_list = ROOT.gGeoManager.GetListOfVolumes()
    pylist = list()
    for i in range(vol_list.GetEntries()):
        vol = vol_list.At(i)
        pylist.append(vol.GetName())
    if 'volDownstreamMagnet' in pylist: return 1
    else: return 2


def zoom(xmin=None,xmax=None,ymin=None,ymax=None,zmin=None,zmax=None):
# zoom() will reset to default setting
  for d in ['xmin','xmax','ymin','ymax','zmin','zmax']:
     if eval(d): h['c'+d]=eval(d)
     else: h['c'+d]=h[d]
  h['xz'].GetXaxis().SetRangeUser(h['czmin'],h['czmax'])
  h['yz'].GetXaxis().SetRangeUser(h['czmin'],h['czmax'])
  h['xz'].GetYaxis().SetRangeUser(h['cxmin'],h['cxmax'])
  h['yz'].GetYaxis().SetRangeUser(h['cymin'],h['cymax'])
  tc = h['simpleDisplay'].cd(1)
  tc.Update()
  tc = h['simpleDisplay'].cd(2)
  tc.Update()
  h['simpleDisplay'].Update()

def dumpVeto():
    muHits = {10:[],11:[]}
    for aHit in eventTree.Digi_AdvMuFilterHits:
         if not aHit.isValid(): continue
         s = aHit.GetDetectorID()//10000
         if s>1: continue
         p = (aHit.GetDetectorID()//1000)%10
         bar = (aHit.GetDetectorID()%1000)%60
         plane = s*10+p
         muHits[plane].append(aHit)
    for plane in [10,11]:
        for aHit in muHits[plane]:
          S =aHit.GetAllSignals(False,False)
          txt = ""
          for x in S:
              if x[1]>0: txt+=str(x[1])+" "
          print(plane, (aHit.GetDetectorID()%1000)%60, txt)

# decode MuFilter detID
def MuFilter_PlaneBars(detID):
         s = detID//10000
         l  = (detID%10000)//1000  # plane number
         bar = (detID%1000)
         if s>2:
             l=2*l
             if bar>59:
                  bar=bar-60
                  if l<6: l+=1
         return {'station':s,'plane':l,'bar':bar}

def checkOtherTriggers(event,deadTime = 100,debug=False):
      T0 = event.EventHeader.GetEventTime()
      N = event.EventHeader.GetEventNumber()
      Nprev = 1
      rc = event.GetEvent(N-Nprev)
      dt = T0 - event.EventHeader.GetEventTime()
      otherFastTrigger = False
      otherAdvTrigger = False
      tightNoiseFilter = False
      while dt < deadTime:
         otherFastTrigger = False
         for x in event.EventHeader.GetFastNoiseFilters():
             if debug: print('fast:', x.first, x.second )
             if x.second and not x.first == 'Veto_Total': otherFastTrigger = True
         otherAdvTrigger = False
         for x in event.EventHeader.GetAdvNoiseFilters():
             if debug: print('adv:', x.first, x.second )
             if x.second and not x.first == 'VETO_Planes': otherAdvTrigger = True
         if debug: print('pre event ',Nprev,dt,otherFastTrigger,otherAdvTrigger)
         if otherFastTrigger and otherAdvTrigger:
             rc = event.GetEvent(N)
             return otherFastTrigger, otherAdvTrigger, tightNoiseFilter, Nprev, dt
         Nprev+=1
         rc = event.GetEvent(N-Nprev)
         dt = T0 - event.EventHeader.GetEventTime()
      Nprev = 1
      rc = event.GetEvent(N-Nprev)
      dt = T0 - event.EventHeader.GetEventTime()
      while dt < deadTime:
         hits = {1:0,0:0}
         for aHit in event.Digi_AdvMuFilterHits:
            Minfo = MuFilter_PlaneBars(aHit.GetDetectorID())
            s,l,bar = Minfo['station'],Minfo['plane'],Minfo['bar']
            if s>1: continue
            allChannels = aHit.GetAllSignals(False,False)
            hits[l]+=len(allChannels)
         noiseFilter0 = (hits[0]+hits[1])>4.5
         noiseFilter1 = hits[0]>0 and hits[1]>0
         if debug: print('veto hits:',hits)
         if noiseFilter0 and noiseFilter1: 
            tightNoiseFilter = True
            rc = event.GetEvent(N)
            return otherFastTrigger, otherAdvTrigger, tightNoiseFilter, Nprev-1, dt
         Nprev+=1
         rc = event.GetEvent(N-Nprev)
         dt = T0 - event.EventHeader.GetEventTime()
      if Nprev>1: 
            rc = event.GetEvent(N-Nprev+1)
            dt = T0 - event.EventHeader.GetEventTime()
      rc = event.GetEvent(N)
      return otherFastTrigger, otherAdvTrigger, tightNoiseFilter, Nprev-1, dt
      
def cleanTracks():
    OT = sink.GetOutTree()
    listOfDetIDs = {}
    n = 0
    for aTrack in OT.Reco_MuonTracks:
        listOfDetIDs[n] = []
        for i in range(aTrack.getNumPointsWithMeasurement()):
           M =  aTrack.getPointWithMeasurement(i)
           R =  M.getRawMeasurement()
           listOfDetIDs[n].append(R.getDetId())
           if R.getDetId()>0: listOfDetIDs[n].append(R.getDetId()-1)
           listOfDetIDs[n].append(R.getDetId()+1)
        n+=1
    uniqueTracks = []
    for n1 in range( len(listOfDetIDs) ):
       unique = True
       for n2 in range( len(listOfDetIDs) ):
             if n1==n2: continue
             I = set(listOfDetIDs[n1]).intersection(listOfDetIDs[n2])
             if len(I)>0:  unique = False
       if unique: uniqueTracks.append(n1)
    if len(uniqueTracks)>1: 
         for n1 in range( len(listOfDetIDs) ): print(listOfDetIDs[n1])
    return uniqueTracks

def timingOfEvent(makeCluster=False,debug=False):
   firstScifi_z = 300*u.cm
   TDC2ns = 1E9/160.316E6
   ut.bookHist(h,'evTimeDS','cor time of hits;[ns]',70,-5.,30)
   ut.bookHist(h,'evTimeScifi','cor time of hits blue DS red Scifi;[ns]',70,-5.,30)
   ut.bookCanvas(h,'tevTime','cor time of hits',1024,768,1,1)
   h['evTimeScifi'].SetLineColor(ROOT.kRed)
   h['evTimeDS'].SetLineColor(ROOT.kBlue)
   h['evTimeScifi'].SetStats(0)
   h['evTimeDS'].SetStats(0)
   h['evTimeScifi'].SetLineWidth(2)
   h['evTimeDS'].SetLineWidth(2)
   if makeCluster: trackTask.scifiCluster()
   meanXY = {}
   for siCl in trackTask.clusScifi:
       detID = siCl.GetFirst()
       s = detID//1000000
       isVertical = detID%1000000//100000
       siCl.GetPosition(A,B)
       z=(A[2]+B[2])/2.
       pos = (A[1]+B[1])/2.
       L = abs(A[0]-B[0])/2.
       if isVertical:
          pos = (A[0]+B[0])/2.
          L = abs(A[1]-B[1])/2.
       corTime = geo.modules['AdvTarget'].GetCorrectedTime(detID, siCl.GetTime(), 0) - (z-firstScifi_z)/u.speedOfLight
       h['evTimeScifi'].Fill(corTime)
       if debug: print(detID,corTime,pos)
   for aHit in eventTree.Digi_AdvMuFilterHits:
       detID = aHit.GetDetectorID()
       if not detID//10000==3: continue
       if aHit.isVertical(): nmax = 1
       else: nmax=2
       geo.modules['AdvMuFilter'].GetPosition(detID,A,B)
       z=(A[2]+B[2])/2.
       pos = (A[1]+B[1])/2.
       L = abs(A[0]-B[0])/2.
       if isVertical: 
          pos = (A[0]+B[0])/2.
          L = abs(A[1]-B[1])/2.
       for i in range(nmax):
            corTime = geo.modules['AdvMuFilter'].GetCorrectedTime(detID, i, aHit.GetTime(i)*TDC2ns, 0)- (z-firstScifi_z)/u.speedOfLight
            h['evTimeDS'].Fill(corTime)
            if debug: print(detID,i,corTime,pos)
   tc=h['tevTime'].cd()
   h['evTimeScifi'].Draw()
   h['evTimeDS'].Draw('same')
   tc.Update()
def mufiNoise():
  for s in range(1,4): 
    ut.bookHist(h,'mult'+str(s),'hit mult for system '+str(s),100,-0.5,99.5)
    ut.bookHist(h,'multb'+str(s),'hit mult per bar for system '+str(s),20,-0.5,19.5)
    ut.bookHist(h,'res'+str(s),'residual system '+str(s),20,-10.,10.)
  OT = sink.GetOutTree()
  N=0
  for event in eventTree:
       N+=1
       if N%1000==0: print(N)
       OT.Reco_MuonTracks.Delete()
       rc = trackTask.ExecuteTask("Scifi")
       for aTrack in OT.Reco_MuonTracks:
           mom    = aTrack.getFittedState().getMom()
           pos      = aTrack.getFittedState().getPos()
           if not aTrack.getFitStatus().isFitConverged(): continue
           mult = {1:0,2:0,3:0}
           for aHit in eventTree.Digi_AdvMuFilterHits:
              if not aHit.isValid(): continue
              s = aHit.GetDetectorID()//10000
              S = aHit.GetAllSignals(False,False)
              rc = h['multb'+str(s)].Fill(len(S))
              mult[s]+=len(S)
              if s==2 or s==1:
                 geo.modules['AdvMuFilter'].GetPosition(aHit.GetDetectorID(),A,B)
                 y = (A[1]+B[1])/2.
                 zEx = (A[2]+B[2])/2.
                 lam      = (zEx-pos.z())/mom.z()
                 Ey        = pos.y()+lam*mom.y()
                 rc = h['res'+str(s)].Fill(Ey-y)
           for s in mult: rc = h['mult'+str(s)].Fill(mult[s])
  ut.bookCanvas(h,'noise','',1200,1200,2,3)
  for s in range(1,4):
   tc = h['noise'].cd(s*2-1)
   tc.SetLogy(1)
   h['mult'+str(s)].Draw()
   h['noise'].cd(s*2)
   h['multb'+str(s)].Draw()
  ut.bookCanvas(h,'res','',600,1200,1,3)
  for s in range(1,4):
   tc = h['res'].cd(s)
   h['res'+str(s)].Draw()

def firstTimeStamp(event):
        tmin = [1E9,'']
        digis = [event.Digi_AdvMuFilterHits,event.Digi_AdvTargetHits]
        for digi in event.Digi_AdvTargetHits:
               dt = digi.GetTime()
               if dt<tmin[0]:
                    tmin[0]=dt
                    tmin[1]=digi
        for digi in event.Digi_AdvMuFilterHits:
           for t in digi.GetAllTimes():    # will not give time if QDC<0!
               dt = t.second
               if dt<tmin[0]:
                    tmin[0]=dt
                    tmin[1]=digi
        return tmin

def dumpChannels(D='Digi_AdvMuFilterHits'):
     X = eval("eventTree."+D)
     text = {}
     for aHit in X:
         side = 'L'
         txt = "%8i"%(aHit.GetDetectorID())
         for k in range(aHit.GetnSiPMs()*aHit.GetnSides()):
              qdc = aHit.GetSignal(k)
              if qdc < -900: continue
              i = k
              if not k<aHit.GetnSiPMs():
                   i = k-aHit.GetnSiPMs()
                   if side == 'L':
                        txt += " | "
                        side = 'R'
              txt+= "  %2i:%4.1F  "%(i,qdc)
         text[aHit.GetDetectorID()] = txt
     keys = list(text.keys())
     keys.sort()
     for k in keys: print(text[k])

def fillNode(node):
   xNodes = {'UpstreamBar', 'VetoBar', 'hor'}
   proj = {'X':0,'Y':1}
   color = ROOT.kBlack
   thick = 5
   for p in proj:
      if node+p in h:
         X = h[node+p]
         if 'Veto' in node:
            color = ROOT.kRed+1
         if 'Downstream' in node:
            thick = 5
         c = proj[p]
         h[ 'simpleDisplay'].cd(c+1)
         X.SetFillColor(color)
         X.SetLineColor(color)
         X.SetLineWidth(thick)
         X.Draw('f&&same')
         X.Draw('same')   

def drawInfo(pad, k, run, event, timestamp,moreEventInfo=[]):
   drawLogo = True
   drawText = True
   if drawLogo:
      padLogo = ROOT.TPad("logo","logo",0.1,0.1,0.2,0.3)
      padLogo.SetFillStyle(4000)
      padLogo.SetFillColorAlpha(0, 0)
      padLogo.Draw()
      logo = ROOT.TImage.Open('$SNDSW_ROOT/shipLHC/Large__SND_Logo_black_cut.png')
      logo.SetConstRatio(True)
      logo.DrawText(0, 0, 'SND', 98)
      padLogo.cd()
      #logo.Draw()
      pad.cd(k)

   if drawText:
    if k==1 or len(moreEventInfo)<5:
      runNumber = eventTree.EventHeader.GetRunId()
      if eventTree.GetBranch('MCTrack'):
        timestamp_start = False
      else:
        timestamp_start = getStartTime(runNumber)
        if  timestamp_start:
           TDC2ns = 6.23768   #conversion factor from 160MHz clock to ns
           timestamp_s = timestamp * TDC2ns * 1E-9
           timestamp_event = int(timestamp_start + timestamp_s)
           time_event = datetime.fromtimestamp(timestamp_event)
      padText = ROOT.TPad("info","info",0.19,0.1,0.6,0.3)
      padText.SetFillStyle(4000)
      padText.Draw()
      padText.cd()
      textInfo = ROOT.TLatex()
      textInfo.SetTextAlign(11)
      textInfo.SetTextFont(42)
      textInfo.SetTextSize(.15)
      textInfo.DrawLatex(0, 0.6, 'Proposed AdvSND Experiment, CERN')
      if hasattr(eventTree.EventHeader,'GetEventNumber'): N = eventTree.EventHeader.GetEventNumber()
      else: N = event
      textInfo.DrawLatex(0, 0.4, 'Run / Event: '+str(run)+' / '+str(N))
      if timestamp_start:
           textInfo.DrawLatex(0, 0.2, 'Time (GMT): {}'.format(time_event))
      else:
           textInfo.DrawLatex(0, 0.2, 'MC simulation')
      pad.cd(k)
    elif options.extraInfo:
      padText = ROOT.TPad("info","info",0.29,0.12,0.9,0.35)
      padText.SetFillStyle(4000)
      padText.Draw()
      padText.cd()
      textInfo = ROOT.TLatex()
      textInfo.SetTextAlign(11)
      textInfo.SetTextFont(42)
      textInfo.SetTextSize(.1)
      textInfo.SetTextColor(ROOT.kMagenta+2)
      dely = 0.12
      ystart = 0.85
      for i in range(7):
        textInfo.DrawLatex(0.4, 0.9-dely*i, moreEventInfo[i])
      pad.cd(k)

