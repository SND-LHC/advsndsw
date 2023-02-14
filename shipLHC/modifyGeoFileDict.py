#!/usr/bin/env python
import ROOT,os
from rootpyPickler import Pickler
from rootpyPickler import Unpickler
import shipunit as u
import saveBasicParameters

path = os.environ['EOSSHIP']+"/eos/experiment/sndlhc/convertedData/physics/2022/"
supportedGeoFiles = ["geofile_sndlhc_TI18_V1_06July2022.root","geofile_sndlhc_TI18_V2_12July2022.root","geofile_sndlhc_TI18_V3_08August2022.root",
                     "geofile_sndlhc_TI18_V4_10August2022.root","geofile_sndlhc_TI18_V5_14August2022.root","geofile_sndlhc_TI18_V6_08October2022.root",
                     "geofile_sndlhc_TI18_V7_22November2022.root"]
            
def modifyDicts():
   for g in supportedGeoFiles:
         os.system('xrdcp -f '+os.environ['EOSSHIP']+'/eos/experiment/sndlhc/convertedData/commissioning/TI18/'+g+' '+g)
         fg  = ROOT.TFile(g)
         pkl = Unpickler(fg)
         sGeo = pkl.load('ShipGeo')
         fg.Close()
         
         setattr(sGeo.MuFilter,'DsPropSpeed',14.9*u.cm/u.nanosecond)
         sGeo.MuFilter['DsPropSpeed'] = 14.9*u.cm/u.nanosecond
       #time delay corrections first order, only for DS at the moment
         setattr(sGeo.MuFilter,'DSTcorslope',0.09)
         sGeo.MuFilter['DSTcorslope'] = 0.09
         constants = [-5.59,-5.60,-5.74,-5.36,-4.84,-5.04,-5.97,-6.05,-6.25,-6.39,-5.39,-5.68,-5.43,-5.48,-5.72,-5.87,-4.95,-5.03,-5.05,-5.19]
         for i in range(len(constants)): 
            setattr(sGeo.MuFilter,'DSTcorC'+str(i),constants[i])
            sGeo.MuFilter['DSTcorC'+str(i)] = constants[i]
         print('save',g)
         saveBasicParameters.execute(g,sGeo)
