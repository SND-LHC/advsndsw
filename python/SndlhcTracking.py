import ROOT
from array import array
import shipunit as u
A,B = ROOT.TVector3(),ROOT.TVector3()

ROOT.gInterpreter.Declare("""
#include <KalmanFitterInfo.h>
#include <Track.h>
#include <MeasuredStateOnPlane.h>
#include <stddef.h>     

const genfit::MeasuredStateOnPlane& getFittedState(genfit::Track* theTrack, int nM){
      try{
        return theTrack->getFittedState(nM);
      }
      catch(genfit::Exception& e){
        std::cerr<<"Exception "<< e.what() <<std::endl;
        const genfit::MeasuredStateOnPlane* state(NULL);
        return *state;
      }
}
""")

class Tracking(ROOT.FairTask):
 " Tracking "
 def Init(self,online=False):
   geoMat = ROOT.genfit.TGeoMaterialInterface()
   bfield = ROOT.genfit.ConstField(0,0,0)   # constant field of zero
   fM = ROOT.genfit.FieldManager.getInstance()
   fM.init(bfield)
   ROOT.genfit.MaterialEffects.getInstance().init(geoMat)
   ROOT.genfit.MaterialEffects.getInstance().setNoEffects()
   lsOfGlobals  = ROOT.gROOT.GetListOfGlobals()
   self.scifiDet = lsOfGlobals.FindObject('Scifi')
   self.mufiDet = lsOfGlobals.FindObject('MuFilter')
   
   # internal storage of clusters
   self.clusScifi   = ROOT.TObjArray(100)
   self.DetID2Key = {}
   self.clusMufi  = ROOT.TObjArray(100)
   
   self.fitter = ROOT.genfit.KalmanFitter()
   self.fitter.setMaxIterations(50)
   #internal storage of fitted tracks
   # output is in genfit::track or sndRecoTrack format
   # default is genfit::Track
   if hasattr(self, "genfitTrack"): pass
   else: self.genfitTrack = True
   
   if self.genfitTrack:
        self.fittedTracks = ROOT.TObjArray(10)
   else:
        self.fittedTracks = ROOT.TClonesArray("sndRecoTrack", 10)
   self.sigmaScifi_spatial = 2*150.*u.um
   self.sigmaMufiUS_spatial = 2.*u.cm
   self.sigmaMufiDS_spatial = 0.3*u.cm
   self.scifi_vsignal = 15.*u.cm/u.ns
   self.Debug = False
   self.ioman = ROOT.FairRootManager.Instance()
   self.sink = self.ioman.GetSink()
   # online mode:         raw data in, converted data in output
   # offline read only:   converted data in, no output
   # offline read/write:  converted data in, converted data out

   # for scifi tracking
   self.nPlanes = 3
   self.nClusters = 5
   self.sigma=150*u.um
   self.maxRes=50
   self.maskPlane=-1
   # for DS tracking
   self.DSnPlanes = 2
   self.DSnHits = 2
   self.nDSPlanesVert  = self.mufiDet.GetConfParI("MuFilter/NDownstreamPlanes")
   self.nDSPlanesHor = self.nDSPlanesVert-1

   if online:
      self.event = self.sink.GetOutTree()
   else: 
      self.event = self.ioman.GetInChain()     # should contain all digis, but not necessarily the tracks and scifi clusters

   self.systemAndPlanes  = {1:2,2:5,3:7}
   return 0
 
 def SetTrackClassType(self,genfitTrack):
     self.genfitTrack = genfitTrack
 
 def FinishEvent(self):
  pass

 def ExecuteTask(self,option='ScifiDS'):
    self.trackCandidates = {}
    if not option.find('DS')<0:
           self.clusMufi.Delete()
           self.dsCluster()
           self.trackCandidates['DS'] = self.DStrack()
    if not option.find('Scifi')<0:
           self.clusScifi.Delete()
           self.scifiCluster()
           self.trackCandidates['Scifi'] = self.Scifi_track()
    i_muon = -1
    for x in self.trackCandidates:
      for aTrack in self.trackCandidates[x]:
           rc = self.fitTrack(aTrack)
           if type(rc)==type(1):
                print('trackfit failed',rc,aTrack)
           else:
                i_muon += 1
                if x=='DS':   rc.SetUniqueID(3)
                if x=='Scifi': rc.SetUniqueID(1)
                # add the tracks
                if self.genfitTrack:
                    self.fittedTracks.Add(rc)
                else:
                    # Load items into snd track class object
                    #if not rc.getFitStatus().isFitConverged(): continue
                    this_track = ROOT.sndRecoTrack(rc)
                    pointTimes = []
                    if x=='DS':
                       for pnt in rc.getPointsWithMeasurement():
                           hitID = pnt.getRawMeasurement().getHitId()
                           aCl = self.clusMufi[hitID]
                           pointTimes.append(aCl.GetTime())
                    if x=='Scifi':
                       for pnt in rc.getPointsWithMeasurement():
                           hitID = pnt.getRawMeasurement().getHitId()
                           aCl = self.clusScifi[hitID]
                           pointTimes.append(aCl.GetTime())
                    this_track.setRawMeasTimes(pointTimes)
                    this_track.setTrackType(rc.GetUniqueID())
                    # Store the track in sndRecoTrack format
                    self.fittedTracks[i_muon] = this_track
            

 def DStrack(self):
    event = self.event
    trackCandidates = []
    clusters = self.clusMufi

    stations = {}
    s = 3
    for plane in range(self.systemAndPlanes[s]): 
          stations[s*10+plane] = {}
    k=-1
    for aCl in clusters:
         k+=1
         detID = aCl.GetFirst()
         p = (detID//1000)%10
         bar = detID%1000
         plane = s*10+p
         if s==3:
           if bar<60 or p==3: plane = s*10+2*p
           else:  plane = s*10+2*p+1
           stations[plane][k] = aCl
    success = True
    pXWithHits = 0
    pYWithHits = 0
    clustPer_pX = {}
    clustPer_pY = {}
    mask_pX = []
    mask_pY = []
    for p in range(30,30+self.systemAndPlanes[s]):
         if p%2==0 and p<36: clustPer_pY[p]=len(stations[p])
         else: clustPer_pX[p]=len(stations[p])
         if len(stations[p])>self.DSnHits or len(stations[p])<1: continue
         if p%2==0 and p<36: pYWithHits+=1
         else: pXWithHits+=1
    if pXWithHits<self.DSnPlanes or pYWithHits<self.DSnPlanes: success = False
    if success:
 # build trackCandidate
      # define planes to mask
      clustPer_pX = dict(sorted(clustPer_pX.items(), key=lambda item: item[1], reverse = True))
      clustPer_pY = dict(sorted(clustPer_pY.items(), key=lambda item: item[1], reverse = True))
      # count planes with clusters
      nX = self.nDSPlanesVert - list(clustPer_pX.values()).count(0)
      nY = self.nDSPlanesHor - list(clustPer_pY.values()).count(0)
      # mask busiest planes until there are at least DSnPlanes planes with clusters left
      for ii in range(nX-self.DSnPlanes):
        if list(clustPer_pX.values())[ii] >=self.DSnHits:
           mask_pX.append(list(clustPer_pX.keys())[ii])
      for ii in range(nY-self.DSnPlanes):
        if list(clustPer_pY.values())[ii] >=self.DSnHits:
           mask_pY.append(list(clustPer_pY.keys())[ii])

      hitlist = {}
      for p in stations:
         if p in mask_pX or p in mask_pY: continue
         for k in stations[p]:
             hitlist[k] = stations[p][k]
      trackCandidates.append(hitlist)
    return trackCandidates

 def Scifi_track(self):
# check for low occupancy and enough hits in Scifi
        event = self.event
        trackCandidates = []
        clusters = self.clusScifi
        stations = {}
        projClusters = {0:{},1:{}}
        for s in range(1,6):
           for o in range(2):
              stations[s*10+o] = []
        k=0      
        for cl in clusters:
            detID = cl.GetFirst()
            s  = detID//1000000
            o = (detID//100000)%10
            if self.maskPlane != s:
                stations[s*10+o].append(detID)
                projClusters[o][detID] = [cl,k]
                k+=1
        nclusters = 0
        check = {}
        ignore = []
        for o in range(2):
            check[o]={}
            for s in range(1,6):
                if len(stations[s*10+o]) > self.nClusters: 
                  ignore.append(s*10+o)
                elif len(stations[s*10+o]) > 0 : 
                  check[o][s] = 1
                  nclusters+=len(stations[s*10+o])
        
        if len(check[0])<self.nPlanes or len(check[1])<self.nPlanes: return trackCandidates
# build trackCandidate
# PR logic, fit straight line in x/y projection, remove outliers. Ignore tilt.
        hitlist = {}
        sortedClusters = {}
        check[0]=0
        check[1]=0
        for o in range(2):
           sortedClusters[o]=sorted(projClusters[o])
           g = ROOT.TGraph()
           n = 0
           mean = {}
           points = {}
           for detID in sortedClusters[o]:
               s  = detID//1000000
               if  (s*10+o) in ignore: continue
               if not s in mean: 
                  mean[s]=[0,0,0]
               projClusters[o][detID][0].GetPosition(A,B)
               z = (A[2]+B[2])/2.
               if o==0: y = (A[1]+B[1])/2.
               else: y = (A[0]+B[0])/2.
               points[detID] = [z,y]
               mean[s][0]+=z
               mean[s][1]+=y
               mean[s][2]+=1
           for s in mean:
              for n in range(2):
                 mean[s][n]=mean[s][n]/mean[s][2]
              g.AddPoint(mean[s][0],mean[s][1])
           rc = g.Fit('pol1','SQ')
           fun = g.GetFunction('pol1')
           for detID in points:
               z = points[detID][0]
               y = points[detID][1]
               res = abs(y-fun.Eval(z))/self.sigma
               if res < self.maxRes:
                 k = projClusters[o][detID][1]
                 hitlist[k] = projClusters[o][detID][0]
                 check[o]+=1
        if check[0]>2 and check[1]>2:
              trackCandidates.append(hitlist)
        return trackCandidates

 def scifiCluster(self):
       clusters = []
       self.DetID2Key.clear()
       hitDict = {}
       for k in range(self.event.Digi_ScifiHits.GetEntries()):
            d = self.event.Digi_ScifiHits[k]
            if not d.isValid(): continue 
            hitDict[d.GetDetectorID()] = k
       hitList = list(hitDict.keys())
       if len(hitList)>0:
              hitList.sort()
              tmp = [ hitList[0] ]
              cprev = hitList[0]
              ncl = 0
              last = len(hitList)-1
              hitvector = ROOT.std.vector("sndScifiHit*")()
              for i in range(len(hitList)):
                   if i==0 and len(hitList)>1: continue
                   c=hitList[i]
                   neighbour = False
                   if (c-cprev)==1:    # does not account for neighbours across sipms
                        neighbour = True
                        tmp.append(c)
                   if not neighbour  or c==hitList[last]:
                        first = tmp[0]
                        N = len(tmp)
                        hitvector.clear()
                        for aHit in tmp: hitvector.push_back( self.event.Digi_ScifiHits[hitDict[aHit]])
                        aCluster = ROOT.sndCluster(first,N,hitvector,self.scifiDet,False)
                        clusters.append(aCluster)
                        if c!=hitList[last]:
                             ncl+=1
                             tmp = [c]
                        elif not neighbour :   # save last channel
                            hitvector.clear()
                            hitvector.push_back( self.event.Digi_ScifiHits[hitDict[c]])
                            aCluster = ROOT.sndCluster(c,1,hitvector,self.scifiDet,False)
                            clusters.append(aCluster)
                   cprev = c
       self.clusScifi.Delete()            
       for c in clusters:  
            self.clusScifi.Add(c)
# map clusters to hit keys
            for nHit in range(self.event.Digi_ScifiHits.GetEntries()):
              if self.event.Digi_ScifiHits[nHit].GetDetectorID()==c.GetFirst():
                 self.DetID2Key[c.GetFirst()] = nHit

 def dsCluster(self):
       clusters = []
       hitDict = {}
       for k in range(self.event.Digi_MuFilterHits.GetEntries()):
            d = self.event.Digi_MuFilterHits[k]
            if (d.GetDetectorID()//10000)<3 or (not d.isValid()): continue
            hitDict[d.GetDetectorID()] = k
       hitList = list(hitDict.keys())
       if len(hitList)>0:
              hitList.sort()
              tmp = [ hitList[0] ]
              cprev = hitList[0]
              ncl = 0
              last = len(hitList)-1
              hitvector = ROOT.std.vector("MuFilterHit*")()
              for i in range(len(hitList)):
                   if i==0 and len(hitList)>1: continue
                   c=hitList[i]
                   neighbour = False
                   if (c-cprev)==1 or (c-cprev)==2:    # allow for one missing channel
                        neighbour = True
                        tmp.append(c)
                   if not neighbour  or c==hitList[last] or c%1000==59:
                        first = tmp[0]
                        N = len(tmp)
                        hitvector.clear()
                        for aHit in tmp: hitvector.push_back( self.event.Digi_MuFilterHits[hitDict[aHit]])
                        aCluster = ROOT.sndCluster(first,N,hitvector,self.mufiDet,False)
                        clusters.append(aCluster)
                        if c!=hitList[last]:
                             ncl+=1
                             tmp = [c]
                        elif not neighbour :   # save last channel
                            hitvector.clear()
                            hitvector.push_back( self.event.Digi_MuFilterHits[hitDict[c]])
                            aCluster = ROOT.sndCluster(c,1,hitvector,self.mufiDet,False)
                            clusters.append(aCluster)
                   cprev = c
       self.clusMufi.Delete()
       for c in clusters:  self.clusMufi.Add(c)

 def patternReco(self):
# very simple for the moment, take all scifi clusters
    trackCandidates = []
    hitlist = {}
    ScifiStations = {}
    for k in range(len(self.event.Cluster_Scifi)):
           hitlist[k] = self.event.Cluster_Scifi[k]
           ScifiStations[hitlist[k].GetFirst()//1000000] = 1
# take fired muonFilter bars if more than 2 SiPMs have fired
    nMin = 1
    MuFiPlanes = {}
    for k in range(self.event.Digi_MuFilterHits.GetEntries()):
         aHit = self.event.Digi_MuFilterHits[k]
         if not aHit.isValid(): continue
         detID = aHit.GetDetectorID()
         sy    = detID//10000
         l       = (detID%10000)//1000  # plane number
         bar = (detID%1000)
         nSiPMs = aHit.GetnSiPMs()
         nSides  = aHit.GetnSides()
         nFired = 0
         for i in range(nSides*nSiPMs):
              if aHit.GetSignal(i) > 0: nFired+=1
         if nMin > nFired: continue
         hitlist[k*1000] = self.event.Digi_MuFilterHits[k]
         MuFiPlanes[sy*100+l] = 1
    if (len(ScifiStations) == 5 or len(MuFiPlanes)>4) and len(hitlist)<20:
           trackCandidates.append(hitlist)
    return trackCandidates

 def fitTrack(self,hitlist):
# hitlist:  clusterID: [A,B] endpoints of scifiCluster
    hitPosLists={}
    trID = 0

    posM    = ROOT.TVector3(0, 0, 0.)
    momM = ROOT.TVector3(0,0,100.)  # default track with high momentum

# approximate covariance
    covM = ROOT.TMatrixDSym(6)
    res = self.sigmaScifi_spatial
    for  i in range(3):   covM[i][i] = res*res
    for  i in range(3,6): covM[i][i] = ROOT.TMath.Power(res / (4.*2.) / ROOT.TMath.Sqrt(3), 2)
    rep = ROOT.genfit.RKTrackRep(13)

# start state
    state = ROOT.genfit.MeasuredStateOnPlane(rep)
    rep.setPosMomCov(state, posM, momM, covM)

# create track
    seedState = ROOT.TVectorD(6)
    seedCov   = ROOT.TMatrixDSym(6)
    rep.get6DStateCov(state, seedState, seedCov)
    theTrack = ROOT.genfit.Track(rep, seedState, seedCov)

# make measurements sorted in z
    unSortedList = {}
    tmpList = {}
    A,B = ROOT.TVector3(),ROOT.TVector3()
    for k in hitlist:
        aCl = hitlist[k]
        if hasattr(aCl,"GetFirst"):
            detID = aCl.GetFirst()
            aCl.GetPosition(A,B)
            detSys  = 1
            if detID>50000: detSys=3
        else:
            detID = aCl.GetDetectorID()
            if detID//10000 < 2: detSys  = 2
            else: detSys  = 3
            self.mufiDet.GetPosition(detID,A,B)
        distance = 0
        tmp = array('d',[A[0],A[1],A[2],B[0],B[1],B[2],distance])
        unSortedList[A[2]] = [ROOT.TVectorD(7,tmp),detID,k,detSys]
    sorted_z=list(unSortedList.keys())
    sorted_z.sort()
    for z in sorted_z:
        tp = ROOT.genfit.TrackPoint() # note how the point is told which track it belongs to
        hitCov = ROOT.TMatrixDSym(7)
        detSys = unSortedList[z][3]
        if detSys==3:      
              res = self.sigmaMufiDS_spatial
              maxDis = 1.0
        elif detSys==2:  
              res = self.sigmaMufiUS_spatial
              maxDis = 5.0
        else:         
              res = self.sigmaScifi_spatial
              maxDis = 0.1
        hitCov[6][6] = res*res
        measurement = ROOT.genfit.WireMeasurement(unSortedList[z][0],hitCov,1,6,tp) # the measurement is told which trackpoint it belongs to
        measurement.setMaxDistance(maxDis)
        measurement.setDetId(unSortedList[z][1])
        measurement.setHitId(unSortedList[z][2])
        tp.addRawMeasurement(measurement) # package measurement in the TrackPoint                                          
        theTrack.insertPoint(tp)  # add point to Track
    if not theTrack.checkConsistency():
        print("track not consistent")
        theTrack.Delete()
        return -2
# do the fit
    self.fitter.processTrack(theTrack) # processTrackWithRep(theTrack,rep,True)
    fitStatus   = theTrack.getFitStatus()
    if self.Debug: print("Fit result: converged chi2 Ndf",fitStatus.isFitConverged(),fitStatus.getChi2(),fitStatus.getNdf())
    if not fitStatus.isFitConverged() and 0>1:
        theTrack.Delete()
        return -1
    if self.Debug: 
        chi2 = fitStatus.getChi2()/(fitStatus.getNdf()+1E-15)
        fittedState = theTrack.getFittedState()
        P = fittedState.getMomMag()
        print("track fitted Ndf #Meas P",fitStatus.getNdf(), theTrack.getNumPointsWithMeasurement(),P)
        for p in theTrack.getPointsWithMeasurement():
            rawM = p.getRawMeasurement()
            info = p.getFitterInfo()
            if not info: continue
            detID = rawM.getDetId()
            print(detID,"weights",info.getWeights()[0],info.getWeights()[1],fitStatus.getNdf())
    return theTrack

 def trackDir(self,theTrack):
      if theTrack.GetUniqueID()>1: return False # for the moment, only the scifi is time calibrated
      fitStatus   = theTrack.getFitStatus()
      if not fitStatus.isFitConverged() : return [100,-100]
      state = ROOT.getFittedState(theTrack,0)
      pos = state.getPos()
# start with first measurement
      M = theTrack.getPointWithMeasurement(0)
      W      = M.getRawMeasurement()
      detID = W.getDetId()
      aHit   = self.event.Digi_ScifiHits[ self.DetID2Key[detID] ]
      self.scifiDet.GetSiPMPosition(detID,A,B)
      X = B-pos
      L0 = X.Mag()/self.scifi_vsignal
      # need to correct for signal propagation along fibre
      clkey  = W.getHitId()
      aCl = self.clusScifi[clkey]
      T0track = aCl.GetTime() - L0
      TZero    = aCl.GetTime()
      self.Tline = ROOT.TGraph()
      for nM in range(theTrack.getNumPointsWithMeasurement()):
            state = ROOT.getFittedState(theTrack,nM)
            if not state: continue
            posM   = state.getPos()
            M = theTrack.getPointWithMeasurement(nM)
            W = M.getRawMeasurement()
            detID = W.getDetId()
            clkey = W.getHitId()
            aCl = self.clusScifi[clkey]
            aHit = self.event.Digi_ScifiHits[ self.DetID2Key[detID] ]
            self.scifiDet.GetSiPMPosition(detID,A,B)
            if aHit.isVertical(): X = B-posM
            else: X = A-posM
            L = X.Mag()/self.scifi_vsignal
         # need to correct for signal propagation along fibre
            corTime = self.scifiDet.GetCorrectedTime(detID, aCl.GetTime(), 0)
            trajLength = (posM-pos).Mag()
            if posM[2] -pos[2] < 0: trajLength = -trajLength
            dT = corTime - L - T0track - trajLength/u.speedOfLight
            self.Tline.AddPoint(trajLength,dT)
      rc = self.Tline.Fit('pol1','SQ')
      fitResult =  rc.Get()
      slope = fitResult.Parameter(1)
      return [slope,slope/(fitResult.ParError(1)+1E-13)]

