#!/usr/bin/env python
import ROOT
import shipRoot_conf
from ShipGeoConfig import ConfigRegistry
from argparse import ArgumentParser

mcEngine     = "TGeant4"

parser = ArgumentParser()

parser.add_argument("-c",   dest="config",   help="configuration file", required=True)
parser.add_argument("-g",   dest="geofile",   help="geo file output name", required=True)
options = parser.parse_args()

shipRoot_conf.configure(0)     # load basic libraries, prepare atexit for python

snd_geo = ConfigRegistry.loadpy(options.config)

# -----Create simulation run----------------------------------------
run = ROOT.FairRunSim()
run.SetName(mcEngine)  # Transport engine
run.SetSink(ROOT.FairRootFileSink(ROOT.TMemFile('output', 'recreate')))  # output file
run.SetUserConfig("g4Config.C") # user configuration file default g4Config.C 
rtdb = run.GetRuntimeDb() 

# -----Create geometry----------------------------------------------
import shipLHC_conf as sndDet_conf
modules = sndDet_conf.configure(run,snd_geo)

# -----Create PrimaryGenerator--------------------------------------
primGen = ROOT.FairPrimaryGenerator()

# -----Particle Gun-----------------------
myPgun = ROOT.FairBoxGenerator(13,1)
primGen.AddGenerator(myPgun)

run.SetGenerator(primGen)
# -----Initialize simulation run------------------------------------
run.Init()

# -----Start run----------------------------------------------------
# run.Run(1)
# ------------------------------------------------------------------------
run.CreateGeometryFile(options.geofile)
# save detector parameters dictionary in geofile
import saveBasicParameters
saveBasicParameters.execute(options.geofile,snd_geo)


def checkOverlaps(removeScifi=False):
 sGeo = ROOT.gGeoManager
 sGeo.SetNmeshPoints(10000)
 sGeo.CheckOverlaps(0.1)  # 1 micron takes 5minutes
 sGeo.PrintOverlaps()
# check subsystems in more detail
 for x in sGeo.GetTopNode().GetNodes(): 
   x.CheckOverlaps(0.0001)
   sGeo.PrintOverlaps()

checkOverlaps()
