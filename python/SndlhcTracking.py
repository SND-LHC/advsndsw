import ROOT
from array import array
import shipunit as u
A,B = ROOT.TVector3(),ROOT.TVector3()

class Tracking(ROOT.FairTask):
 " Tracking "
 def Init(self):
   geoMat =  ROOT.genfit.TGeoMaterialInterface()
   bfield     = ROOT.genfit.ConstField(0,0,0)   # constant field of zero
   fM = ROOT.genfit.FieldManager.getInstance()
   fM.init(bfield)
   ROOT.genfit.MaterialEffects.getInstance().init(geoMat)
   ROOT.genfit.MaterialEffects.getInstance().setNoEffects()
   lsOfGlobals  = ROOT.gROOT.GetListOfGlobals()
   self.scifiDet = lsOfGlobals.FindObject('Scifi')
   self.mufiDet = lsOfGlobals.FindObject('MuFilter')

   self.fitter = ROOT.genfit.KalmanFitter()
   self.fitter.setMaxIterations(50)
   self.sigmaScifi_spatial = 2*150.*u.um
   self.sigmaMufiUS_spatial = 2.*u.cm
   self.sigmaMufiDS_spatial = 0.3*u.cm
   self.Debug = False
   self.ioman = ROOT.FairRootManager.Instance()
   self.sink = self.ioman.GetSink()
   # online mode:         raw data in, converted data in output, Reco_MuonTracks defined in ConvRawData
   # offline read only:   converted data in, no output
   # offline read/write:  converted data in, converted data out
   offlineRO   = False
   online        = False
   offlineRW   = False
   if not self.sink.GetOutTree():     offlineRO = True
   if not self.ioman.GetInTree() and not self.ioman.GetInChain():     online = True
   if not online and not offlineRO: offlineRW = True

   self.makeScifiClusters = False
   self.clusMufi   = ROOT.TObjArray(100);

   if offlineRO or offlineRW:
         remakeClusters = False
         self.event = self.ioman.GetInChain()     # should contain all digis, but not necessarily the tracks and scifi clusters
         if self.event.FindBranch("Cluster_Scifi"):
             if not self.event.GetBranchStatus("Cluster_Scifi*"): remakeClusters = True
         if self.sink.GetOutTree():  # somebody else in charge
            self.kalman_tracks = ROOT.TObjArray(10)
            if not self.event.FindBranch("Cluster_Scifi") or remakeClusters:
                self.clusScifi   = ROOT.TObjArray(100);
                self.makeScifiClusters = True
         else:
            self.kalman_tracks = ROOT.TObjArray(10)
            self.ioman.Register("Reco_MuonTracks", "", self.kalman_tracks, ROOT.kTRUE)  # user asks for tracking, independent if tracks exist in inputfile.
            if not self.event.FindBranch("Cluster_Scifi") or remakeClusters:   # no scifi clusters on input file, create them
               self.makeScifiClusters = True
               self.clusScifi   = ROOT.TObjArray(100);
               self.ioman.Register("Cluster_Scifi","",self.clusScifi,ROOT.kTRUE)

   if online:
         self.event = self.sink.GetOutTree()
         self.kalman_tracks = self.sink.GetOutTree().Reco_MuonTracks

   self.systemAndPlanes  = {1:2,2:5,3:7}
   return 0

 def FinishEvent(self):
  pass

 def ExecuteTask(self,option=''):
    if self.makeScifiClusters:
          self.clusScifi.Clear()
          self.scifiCluster()
          self.event.Cluster_Scifi = self.sink.GetOutTree().Cluster_Scifi
    self.trackCandidates = {}
    if option=='DS':
           self.clusMufi.Clear()
           self.trackCandidates['DS'] = self.DStrack()
    elif option=='Scifi':
           self.trackCandidates['Scifi'] = self.Scifi_track()
    elif option=='ScifiDS':
           self.trackCandidates['DS'] = self.DStrack()
           self.trackCandidates['Scifi'] = self.Scifi_track()
    else:
           self.trackCandidates['comb'] = self.patternReco()
    for x in self.trackCandidates:
      for aTrack in self.trackCandidates[x]:
           rc = self.fitTrack(aTrack)
           if type(rc)==type(1):
                print('trackfit failed',rc,aTrack)
           else:
                if x=='DS':   rc.SetUniqueID(3)
                if x=='Scifi': rc.SetUniqueID(1)
                self.kalman_tracks.Add(rc)

 def DStrack(self,nPlanes = 2, nHits = 2):
    event = self.event
    trackCandidates = []
    if not event.GetBranch("Cluster_Mufi"):
            self.dtCluster()
            clusters = self.clusMufi
    else:
            clusters = self.event.Cluster_Mufi

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
    for p in range(30,30+self.systemAndPlanes[s]):
         if len(stations[p])>nHits or len(stations[p])<1: continue
         if p%2==0 and p<36: pYWithHits+=1
         else: pXWithHits+=1
    if pXWithHits<nPlanes or pYWithHits<nPlanes: success = False
    if success:
 # build trackCandidate
      hitlist = {}
      for p in stations:
         for k in stations[p]:
             hitlist[k] = stations[p][k]
      trackCandidates.append(hitlist)
    return trackCandidates

 def Scifi_track(self,nPlanes = 3, nClusters = 20,sigma=150*u.um,maxRes=50):
# check for low occupancy and enough hits in Scifi
        event = self.event
        trackCandidates = []
        if not event.GetBranch("Cluster_Scifi"):
            self.scifiCluster()
            clusters = self.clusScifi
        else:
            clusters = self.event.Cluster_Scifi
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
            stations[s*10+o].append(detID)
            projClusters[o][detID] = [cl,k]
            k+=1
        nclusters = 0
        check = {}
        for o in range(2):
            check[o]={}
            for s in range(1,6):
                if len(stations[s*10+o]) > 0: check[o][s]=1
                nclusters+=len(stations[s*10+o])
        if len(check[0])<nPlanes or len(check[1])<nPlanes or nclusters > nClusters: return trackCandidates
# build trackCandidate
# PR logic, fit straight line in x/y projection, remove outliers. Ignore tilt.
        hitlist = {}
        sortedClusters = {}
        masked = {}
        check[0]=0
        check[1]=0
        for o in range(2): 
           sortedClusters[o]=sorted(projClusters[o])
           g = ROOT.TGraph()
           n = 0
           for detID in sortedClusters[o]:
               projClusters[o][detID][0].GetPosition(A,B)
               z = (A[2]+B[2])/2.
               if o==0: y = (A[1]+B[1])/2.
               else: y = (A[0]+B[0])/2.
               g.SetPoint(n,z,y)
               n+=1
           rc = g.Fit('pol1','SQ')
           fun = g.GetFunction('pol1')
           masked[o] = []
           for i in range(n):
               z = g.GetPointX(i)
               res = abs(g.GetPointY(i)-fun.Eval(z))/sigma
               if res < maxRes:
                 detID = sortedClusters[o][i]
                 k = projClusters[o][detID][1]
                 hitlist[k] = projClusters[o][detID][0]
                 check[o]+=1
        if check[0]>2 and check[1]>2:
              trackCandidates.append(hitlist)
        return trackCandidates

 def scifiCluster(self):
       if self.sink.GetOutTree().GetBranch("Cluster_Scifi"):  self.sink.GetOutTree().Cluster_Scifi.Delete()
       clusters = []
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
       for c in clusters:  self.clusScifi.Add(c)

 def dtCluster(self):
       if self.sink.GetOutTree().GetBranch("Cluster_DT"):  self.sink.GetOutTree().Cluster_DT.Delete()
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
