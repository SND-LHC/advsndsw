from array import array
import ROOT,os
from argparse import ArgumentParser
import atexit
def pyExit():
       os.system('kill '+str(os.getpid()))
atexit.register(pyExit)


parser = ArgumentParser()
parser.add_argument("-g", "--geometry", dest="geofile", help="input geometry file", required=True)
parser.add_argument("-sz", dest="startz", help="start z", type=float,required=True)
parser.add_argument("-ez", dest="endz", help="end z", type=float,default=True)
parser.add_argument("-sx", dest="startx", help="start x", type=float,default=-30)
parser.add_argument("-ex", dest="endx", help="end x", type=float,default=-30)
parser.add_argument("-sy", dest="starty", help="start y", type=float,default=30)
parser.add_argument("-ey", dest="endy", help="end y", type=float,default=30)

options = parser.parse_args()

Geniegen = ROOT.GenieGenerator()
import SndlhcGeo
geo = SndlhcGeo.GeoInterface(options.geofile)

# get dimensions by running getGeoInformation on the geofile 
# python $FAIRSHIP/macro/getGeoInformation.py -g geofile_full.conical.Genie-TGeant4.root

# python -i $FAIRSHIP/macro/run_simScript.py --Genie -f /eos/experiment/ship/data/GenieEvents/genie-nu_mu.root
start=array('d',[options.startx,options.starty,options.startz])
end=array('d',[options.endx,options.endy,options.endz])
mparam=array('d',[0,0,0,0,0,0,0,0,0,0,0,0])
Geniegen.MeanMaterialBudget(start, end, mparam)
print(mparam[8], " equivalent interaction length fraction")
print(mparam[1], " equivalent rad length fraction")

