import ROOT
from argparse import ArgumentParser
import os

import SndlhcMuonReco

parser = ArgumentParser()
parser.add_argument("-f", "--inputFile", dest="inputFile", help="single input file", required=True)
parser.add_argument("-g", "--geoFile", dest="geoFile", help="geofile", required=True)
parser.add_argument("-o", "--withOutput", dest="withOutput", help="persistent output", action='store_true',default=False)
parser.add_argument("-s", "--saveTo", dest="outPath", help="output storage path", type=str,default="",required=False)
parser.add_argument("-par", "--parFile", dest="parFile", help="parameter file", required=False, default="${SNDSW_ROOT}/python/TrackingParams.xml")
parser.add_argument("-c", "--case", dest="trackingCase", help="type of tracks to build. Should match the 'tracking_case' name in parFile, use quotes", required=True)
parser.add_argument("-n", "--nEvents", dest="nEvents",  type=int, help="number of events to process", default=1100000)
parser.add_argument("-i", "--firstEvent",dest="firstEvent",  help="First event of input file to use", required=False,  default=0, type=int)

options = parser.parse_args()

import SndlhcGeo
geo = SndlhcGeo.GeoInterface(options.geoFile)

lsOfGlobals = ROOT.gROOT.GetListOfGlobals()
lsOfGlobals.Add(geo.modules['Scifi'])
lsOfGlobals.Add(geo.modules['MuFilter'])

x = options.inputFile
filename = x[x.rfind('/')+1:]
if x.rfind('/run_')>0:
   runN = x[x.rfind('/run_')+4:x.rfind('/run_')+11]
else: runN = ""
outFileName = options.outPath+filename.replace('.root',runN+'_muonReco.root')

fullPath = options.inputFile
if options.inputFile.find('/eos')==0:
     fullPath = os.environ['EOSSHIP']+options.inputFile
F = ROOT.TFile.Open(fullPath)

if options.withOutput:
  outFile = ROOT.TFile(outFileName, 'RECREATE')
else:
  outFile = ROOT.TMemFile(outFileName,'CREATE')

treename = None
for test_treename in ["cbmsim", "rawConv"] :
     if hasattr(F, test_treename) :
          treename = test_treename
          break

if treename is None :
     raise RuntimeError("File {0} contains no object with a valid SND@LHC TTree name".format(fullPath))

fairRootManager = ROOT.FairRootManager.Instance()
fairRootManager.SetTreeName(treename)

run = ROOT.FairRunAna()
print("Initialized FairRunAna")

source = ROOT.FairFileSource(F)
run.SetSource(source)

sink = ROOT.FairRootFileSink(outFile)
run.SetSink(sink)

muon_reco_task = SndlhcMuonReco.MuonReco()
run.AddTask(muon_reco_task)
w = ROOT.TStopwatch()

# Set the parameter file - must be called before Init()
muon_reco_task.SetParFile(options.parFile)
muon_reco_task.SetTrackingCase(options.trackingCase)

run.Init()

# Set number of events to process
nEvents = min( options.nEvents, source.GetEntries())
print('Number of events for task', nEvents)

run.Run(options.firstEvent, options.firstEvent + nEvents)
w.Stop()
print('Real Time:', w.RealTime(), ' CPU time: ', w.CpuTime())
print("Done running")

