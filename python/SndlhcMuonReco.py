import ROOT
import numpy as np
import scipy.ndimage
import warnings
from array import array
import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt
import shipunit as unit

def hit_finder(slope, intercept, box_centers, box_ds, tol = 0.) :
    """ Finds hits intersected by Hough line """
    
    # First check if track at center of box is within box limits
    d = np.abs(box_centers[0,:,1] - (box_centers[0,:,0]*slope + intercept))
    center_in_box = d < (box_ds[0,:,1]+tol)/2.

    # Now check if, assuming line is not in box at box center, the slope is large enough for line to clip the box at corner
    clips_corner = np.abs(slope) > np.abs((d - (box_ds[0,:,1]+tol)/2.)/(box_ds[0,:,0]+tol)/2.)
    
    # If either of these is true, line goes through hit:
    hit_mask = np.logical_or(center_in_box, clips_corner)

    # Return indices
    return np.where(hit_mask)[0]

class hough() :
    """ Hough transform implementation """

    def __init__(self, n_yH, yH_range, n_xH, xH_range, z_offset, Hformat, det_Zlen, squaretheta = False, smooth = True) :

        self.n_yH = n_yH
        self.n_xH = n_xH

        self.yH_range = yH_range
        self.xH_range = xH_range

        self.z_offset = z_offset
        self.HoughSpace_format = Hformat
        
        self.det_Zlen = det_Zlen

        self.smooth = smooth

        self.yH_bins = np.linspace(self.yH_range[0], self.yH_range[1], n_yH)
        if not squaretheta :
            self.xH_bins = np.linspace(self.xH_range[0], self.xH_range[1], n_xH)
        else :
            self.xH_bins = np.linspace(np.sign(self.xH_range[0])*(self.xH_range[0]**0.5), np.sign(self.xH_range[1])*(self.xH_range[1]**0.5), n_xH)
            self.xH_bins = np.sign(self.xH_bins)*np.square(self.xH_bins)
        
        self.cos_thetas = np.cos(self.xH_bins)
        self.sin_thetas = np.sin(self.xH_bins)
        
        self.xH_i = np.array(list(range(n_xH)))

    def fit(self, hit_collection, draw = False, weights = None) :

        self.accumulator = np.zeros((self.n_yH, self.n_xH))
        for i_hit, hit in enumerate(hit_collection) :
            shifted_hitZ = hit[0] - self.z_offset
            if self.HoughSpace_format == 'normal':
                 hit_yH = shifted_hitZ*self.cos_thetas + hit[1]*self.sin_thetas
            elif self.HoughSpace_format == 'linearSlopeIntercept':
                 hit_yH = hit[1] - shifted_hitZ*self.xH_bins
            elif self.HoughSpace_format== 'linearIntercepts':
                 hit_yH = (self.det_Zlen*hit[1] - shifted_hitZ*self.xH_bins)/(self.det_Zlen - shifted_hitZ)
            out_of_range = np.logical_and(hit_yH > self.yH_range[0], hit_yH < self.yH_range[1]) 
            hit_yH_i = np.floor((hit_yH[out_of_range] - self.yH_range[0])/(self.yH_range[1] - self.yH_range[0])*self.n_yH).astype(np.int)

            if weights is not None :
                self.accumulator[hit_yH_i, self.xH_i[out_of_range]] += weights[i_hit]
            else :
                self.accumulator[hit_yH_i, self.xH_i[out_of_range]] += 1

        # Smooth accumulator
        if self.smooth :
            self.accumulator = scipy.ndimage.gaussian_filter(self.accumulator, 3)

        # This might be useful for debugging, but leave out for now.
        if draw :
            pass
            '''
            plt.figure()
            plt.imshow(self.accumulator, origin = "lower", extent = [self.xH_range[0], self.xH_range[-1], self.yH_range[0], self.yH_range[-1]], aspect = "auto")
            if self.HoughSpace_format == 'normal':
               plt.xlabel(r"$\theta$ [rad]")
               plt.ylabel("r [cm]")
            elif self.HoughSpace_format == 'linearSlopeIntercept':
               plt.xlabel("slope")
               plt.ylabel("intercept @ 1st plane [cm]")
            elif self.HoughSpace_format == 'linearIntercepts':
               plt.xlabel("intercept @ last plane [cm]")
               plt.ylabel("intercept @ 1st plane [cm]")
            plt.tight_layout()
            plt.show()
            '''
        
        if self.smooth == True:
          # In case of multiple occurrences of the maximum values, argmax returns
          # the indices corresponding to the first occurrence(along 1st axis).
          # With smootihg it is very unlikely to have 2 bins with same Nentries
          i_max = np.unravel_index(self.accumulator.argmax(), self.accumulator.shape)
        else:
          # In case there are more than 1 bins with the maximal Nentries, choose the one having most entries
          # after smoothing in an sub-area of the found max. This approach is choosen since
          # smoothing the whole accumulator array is time consuming for large Nbins
          maxima = np.argwhere(self.accumulator == np.amax(self.accumulator))
          if len(maxima) == 1:
            i_max = maxima[0]
          else:
            sigma = 1
            truncate = 4
            maxima_smooth = []
            while not len(maxima_smooth) == 1 and len(maxima_smooth) <= len(maxima) and sigma < 4:
                 smooth_max = []
                 maxima_smooth = []
                 at  = 2*int(sigma*truncate+0.5)+1
                 up  = at + 1
                 low = at - 1
                 for item in maxima:
                     if item[0] > low and item[1] > low:
                        subset = scipy.ndimage.gaussian_filter(self.accumulator[item[0]-at:item[0]+up,item[1]-at:item[1]+up], sigma)
                     elif item[0] < at and item[1] > low:
                        subset = scipy.ndimage.gaussian_filter(self.accumulator[0:item[0]+up,item[1]-at:item[1]+up], sigma)
                     elif item[0] > low and item[1] < at:
                        subset = scipy.ndimage.gaussian_filter(self.accumulator[item[0]-at:item[0]+up,0:item[1]+up], sigma)
                     else:
                        subset = scipy.ndimage.gaussian_filter(self.accumulator[0:item[0]+up,0:item[1]+up], sigma)
                     smooth_max.append(np.amax(subset))
                 smooth_max = np.asarray(smooth_max)
                 many = np.argwhere(smooth_max == np.amax(smooth_max))
                 for i in many: 
                   for items in maxima[i]:
                     i_max = items
                   maxima_smooth.append(maxima[i])
                 if len(maxima_smooth) > 1:
                    sigma += 1
            # In case there are still more than 1 bins with the maximal Nentries, choose the one closest to middle in xH axis
            # i.e. the bin corresponding to the smallest abs(track_slope), thus enhancing track reconstruction for IP1 origin
            if len(maxima_smooth) > 1:
               if not self.HoughSpace_format == 'linearIntercepts':
                  i_x = min([x[1] for x in maxima], key=lambda b: abs(b-self.n_xH/2.))
                  for im in maxima:
                      if im[1] == i_x : i_max = im

        found_yH = self.yH_bins[i_max[0]]
        found_xH = self.xH_bins[i_max[1]]
        
        if self.HoughSpace_format == 'normal':
           slope = -1./np.tan(found_xH)
           interceptShift = found_yH/np.sin(found_xH)
           intercept = (np.tan(found_xH)*interceptShift + self.z_offset)/np.tan(found_xH)
        elif self.HoughSpace_format == 'linearSlopeIntercept':
           slope = found_xH
           intercept = found_yH - slope*self.z_offset
        elif self.HoughSpace_format == 'linearIntercepts':
           slope = (found_xH - found_yH)/self.det_Zlen
           intercept = found_yH - slope*self.z_offset
        
        return (slope, intercept)

    def fit_randomize(self, hit_collection, hit_d, n_random, draw = False, weights = None) :
        success = True
        if not len(hit_collection) :
            return (-1, -1, [[],[]], [], False)

        # Randomize hits
        if (n_random > 0) :
            random_hit_collection = []
            for i_random in range(n_random) :
                random_hits = np.random.uniform(size = hit_collection.shape) - 0.5
                random_hits *= hit_d
                random_hits += hit_collection
                random_hit_collection.append(random_hits)

            random_hit_collection = np.concatenate(random_hit_collection)
            if weights is not None :
                weights = np.tile(weights, n_random)

            fit = self.fit(random_hit_collection, draw, weights)
        else :
            fit = self.fit(hit_collection, draw, weights)

        return fit

def numPlanesHit(systems, detector_ids) :
    scifi_stations = []
    mufi_ds_planes = []
    mufi_us_planes = []

    scifi_stations.append( detector_ids[systems == 0]//1000000 )
    mufi_ds_planes.append( (detector_ids[systems == 3]%10000)//1000 )
    mufi_us_planes.append( (detector_ids[systems == 2]%10000)//1000 )

    return len(np.unique(scifi_stations)) + len(np.unique(mufi_ds_planes)) + len(np.unique(mufi_us_planes))
    
class MuonReco(ROOT.FairTask) :
    " Muon reconstruction "

    def Init(self) :
        
        print("Initializing muon reconstruction task!")
        self.lsOfGlobals  = ROOT.gROOT.GetListOfGlobals()
        self.scifiDet = self.lsOfGlobals.FindObject('Scifi')
        self.mufiDet = self.lsOfGlobals.FindObject('MuFilter')
        self.ioman = ROOT.FairRootManager.Instance()

        # Pass input data through to output.
        # self.Passthrough()
        
        # MC or data - needed for hit timing unit
        if self.ioman.GetInTree().GetName() == 'cbmsim': self.isMC = True
        else: self.isMC = False

        # Fetch digi hit collections from online if exist
        sink = self.ioman.GetSink()
        eventTree = None
        if sink:   eventTree = sink.GetOutTree()
        if eventTree:
            self.MuFilterHits = eventTree.Digi_MuFilterHits
            self.ScifiHits       = eventTree.Digi_ScifiHits
        else:
        # Try the FairRoot way 
            self.MuFilterHits = self.ioman.GetObject("Digi_MuFilterHits")
            self.ScifiHits = self.ioman.GetObject("Digi_ScifiHits")

        # If that doesn't work, try using standard ROOT
            if self.MuFilterHits == None :
               warnings.warn("Digi_MuFilterHits not in branch list")
               self.MuFilterHits = self.ioman.GetInTree().Digi_MuFilterHits
            if self.ScifiHits == None :
               warnings.warn("Digi_ScifiHits not in branch list")
               self.ScifiHits = self.ioman.GetInTree().Digi_ScifiHits
        
        if self.MuFilterHits == None :
            raise RuntimeException("Digi_MuFilterHits not found in input file.")
        if self.ScifiHits == None :
            raise RuntimeException("Digi_ScifiHits not found in input file.")
        
        # Initialize event counters in case scaling of events is required
        self.scale = 1
        self.current_event = -1
        self.events_run = 0
        
        # Initialize hough transform - reading parameter xml file
        tree = ET.parse(self.par_file)
        root = tree.getroot()
        
        # Output track in genfit::Track or sndRecoTrack format
        # Check if genfit::Track format is already forced
        if hasattr(self, "genfitTrack"): pass
        else: self.genfitTrack = int(root[0].text)
                
        track_case_exists = False
        for case in root.findall('tracking_case'):
            if case.get('name') == self.tracking_case:
               track_case_exists = True
               # Use SciFi hits or clusters
               self.Scifi_meas = int(case.find('use_Scifi_clust').text)
               # Maximum absolute value of reconstructed angle (+/- 1 rad is the maximum angle to form a triplet in the SciFi)
               max_angle = float(case.find('max_angle').text)
               
               # Hough space representation
               Hspace_format_exists = False 
               for rep in case.findall('Hough_space_format'):
                   if rep.get('name') == self.Hough_space_format:
                      Hspace_format_exists = True
                      # Number of bins per Hough accumulator axes and range
                      ''' xH and yH are the abscissa and ordinate of the Hough parameter space
                          xz and yz represent horizontal and vertical projections 
                          in the SNDLHC physics coord. system '''
                      n_accumulator_yH = int(rep.find('N_yH_bins').text)
                      yH_min_xz = float(rep.find('yH_min_xz').text)
                      yH_max_xz = float(rep.find('yH_max_xz').text)
                      yH_min_yz = float(rep.find('yH_min_yz').text)
                      yH_max_yz = float(rep.find('yH_max_yz').text)
                      n_accumulator_xH = int(rep.find('N_xH_bins').text)
                      xH_min_xz = float(rep.find('xH_min_xz').text)
                      xH_max_xz = float(rep.find('xH_max_xz').text)
                      xH_min_yz = float(rep.find('xH_min_yz').text)
                      xH_max_yz = float(rep.find('xH_max_yz').text)
                   else: continue
               if not Hspace_format_exists:
                  raise RuntimeException("Unknown Hough space format, check naming in parameter xml file.") 
               # Number of random throws per hit
               self.n_random = int(case.find('n_random').text)
               # MuFilter weight. Muon filter hits are thrown more times than scifi
               self.muon_weight = int(case.find('mufi_weight').text)
               # Minimum number of planes hit in each of the downstream muon filter (if muon filter hits used) or scifi (if muon filter hits not used) views to try to reconstruct a muon
               self.min_planes_hit = int(case.find('min_planes_hit').text)

               # Maximum number of muons to find. To avoid spending too much time on events with lots of downstream activity.
               self.max_reco_muons = int(case.find('max_reco_muons').text)

               # How far away from Hough line hits will be assigned to the muon, for Kalman tracking
               self.tolerance = float(case.find('tolerance').text)

               # Which hits to use for track fitting. By default use both scifi and muon filter.
               self.hits_to_fit = case.find('hits_to_fit').text.strip()
               # Which hits to use for triplet condition. By default use only downstream muon system hits.
               self.hits_for_triplet = case.find('hits_for_hough').text.strip()
            else: continue
        if not track_case_exists:
           raise RuntimeException("Unknown tracking case, check naming in parameter xml file.")
                
        # Get sensor dimensions from geometry
        self.MuFilter_ds_dx = self.mufiDet.GetConfParF("MuFilter/DownstreamBarY") # Assume y dimensions in vertical bars are the same as x dimensions in horizontal bars.
        self.MuFilter_ds_dy = self.mufiDet.GetConfParF("MuFilter/DownstreamBarY") # Assume y dimensions in vertical bars are the same as x dimensions in horizontal bars.
        self.MuFilter_ds_dz = self.mufiDet.GetConfParF("MuFilter/DownstreamBarZ")

        self.MuFilter_us_dx = self.mufiDet.GetConfParF("MuFilter/UpstreamBarX")
        self.MuFilter_us_dy = self.mufiDet.GetConfParF("MuFilter/UpstreamBarY")
        self.MuFilter_us_dz = self.mufiDet.GetConfParF("MuFilter/UpstreamBarZ")

        self.Scifi_dx = self.scifiDet.GetConfParF("Scifi/channel_width")
        self.Scifi_dy = self.scifiDet.GetConfParF("Scifi/channel_width")
        self.Scifi_dz = self.scifiDet.GetConfParF("Scifi/epoxymat_z") # From Scifi.cxx This is the variable used to define the z dimension of SiPM channels, so seems like the right dimension to use.
        
        self.Scifi_nPlanes = self.scifiDet.GetConfParI("Scifi/nscifi")
        
        # get the distance between 1st and last detector planes to be used in the track fit.
        # a z_offset is used to shift detector hits so to have smaller Hough parameter space
        # Using geometers measurements! For safety, add a 5-cm-buffer in detector lengths and a 2.5-cm one to z_offset.
        # This is done to account for possible det. position shifts/mismatches going from geom. measurements and sndsw physics CS.
        if self.hits_for_triplet.find('sf') >= 0 and self.hits_for_triplet.find('ds') >= 0:
           det_Zlen = (self.mufiDet.GetConfParF("MuFilter/Muon9Dy") - self.scifiDet.GetConfParF("Scifi/Ypos0"))*unit.cm + 5.0*unit.cm
           z_offset = self.scifiDet.GetConfParF("Scifi/Ypos0")*unit.cm - 2.5*unit.cm        
        elif self.hits_for_triplet == 'sf':
           det_Zlen = (self.scifiDet.GetConfParF("Scifi/Ypos4") - self.scifiDet.GetConfParF("Scifi/Ypos0"))*unit.cm + 5.0*unit.cm
           z_offset = self.scifiDet.GetConfParF("Scifi/Ypos0")*unit.cm - 2.5*unit.cm
        elif self.hits_for_triplet == 'ds':
           det_Zlen = (self.mufiDet.GetConfParF("MuFilter/Muon9Dy") - self.mufiDet.GetConfParF("MuFilter/Muon6Dy"))*unit.cm + 5.0*unit.cm
           z_offset = self.mufiDet.GetConfParF("MuFilter/Muon6Dy")*unit.cm - 2.5*unit.cm
        # this use case is not tested with an z offset yet
        if self.tracking_case.find('nu_') >= 0: z_offset = 0*unit.cm 
        #other use cases come here if ever added        
        
        # Initialize Hough transforms for both views:
        if self.Hough_space_format == 'normal':
            # rho-theta representation - must exclude theta range for tracks passing < 3 det. planes
            self.h_ZX = hough(n_accumulator_yH, [yH_min_xz, yH_max_xz], n_accumulator_xH, [-max_angle+np.pi/2., max_angle+np.pi/2.], z_offset, self.Hough_space_format, det_Zlen)
            self.h_ZY = hough(n_accumulator_yH, [yH_min_yz, yH_max_yz], n_accumulator_xH, [-max_angle+np.pi/2., max_angle+np.pi/2.], z_offset, self.Hough_space_format, det_Zlen)
        else:
            self.h_ZX = hough(n_accumulator_yH, [yH_min_xz, yH_max_xz], n_accumulator_xH, [xH_min_xz, xH_max_xz], z_offset, self.Hough_space_format, det_Zlen)
            self.h_ZY = hough(n_accumulator_yH, [yH_min_yz, yH_max_yz], n_accumulator_xH, [xH_min_yz, xH_max_yz], z_offset, self.Hough_space_format, det_Zlen)

        # If only Scifi hits used, no need for accumulator smoothing.
        if self.hits_to_fit == "sf" :
            self.h_ZX.smooth = False
            self.h_ZY.smooth = False
            self.track_type = 11
        elif self.hits_to_fit == "ds": self.track_type = 13
        else : self.track_type = 15
        
        # To keep temporary detector information
        self.a = ROOT.TVector3()
        self.b = ROOT.TVector3()

        # check if track container exists
        sink = self.ioman.GetSink()
        self.event = sink.GetOutTree()
        if self.event:
            self.kalman_tracks = sink.GetOutTree().Reco_MuonTracks 
        else:
        # Now initialize output in genfit::track or sndRecoTrack format
           if self.genfitTrack:
              self.kalman_tracks = ROOT.TObjArray(self.max_reco_muons)
              self.ioman.Register("Reco_MuonTracks", self.ioman.GetFolderName(), self.kalman_tracks, ROOT.kTRUE)
           else:
              self.kalman_tracks = ROOT.TClonesArray("sndRecoTrack", self.max_reco_muons)
              self.ioman.Register("Reco_MuonTracks", "", self.kalman_tracks, ROOT.kTRUE)

        # internal storage of clusters
        if self.Scifi_meas: self.clusScifi = ROOT.TObjArray(100)
        
        # Kalman filter stuff

        geoMat = ROOT.genfit.TGeoMaterialInterface()
        bfield = ROOT.genfit.ConstField(0, 0, 0)
        fM = ROOT.genfit.FieldManager.getInstance()
        fM.init(bfield)
        ROOT.genfit.MaterialEffects.getInstance().init(geoMat)
        ROOT.genfit.MaterialEffects.getInstance().setNoEffects()
        
        self.kalman_fitter = ROOT.genfit.KalmanFitter()
        self.kalman_fitter.setMaxIterations(50)
        self.kalman_sigmaScifi_spatial = self.Scifi_dx / 12**0.5
        self.kalman_sigmaMufiUS_spatial = self.MuFilter_us_dy / 12**0.5
        self.kalman_sigmaMufiDS_spatial = self.MuFilter_ds_dy/ 12**0.5
        
        # Init() MUST return int
        return 0
    
    def SetScaleFactor(self, scale):
        self.scale = scale
        
    def SetParFile(self, file_name):
        self.par_file = file_name
    
    def SetTrackingCase(self, case):
        self.tracking_case = case

    def SetHoughSpaceFormat(self, Hspace_format):
        self.Hough_space_format = Hspace_format
        
    def ForceGenfitTrackFormat(self):
        self.genfitTrack = 1
    
    def Passthrough(self) :
        T = self.ioman.GetInTree()
        
        for x in T.GetListOfBranches():
             obj_name = x.GetName()
             if self.ioman.GetObject(obj_name) == None :
                 continue
             self.ioman.Register(obj_name, self.ioman.GetFolderName(), self.ioman.GetObject(obj_name), ROOT.kTRUE) 

    def Exec(self, opt) :
        self.kalman_tracks.Clear('C')
                
        if self.scale>1:
           self.current_event += 1
           if self.current_event == 0: 
              self.event_to_process = int(ROOT.gRandom.Rndm()*self.scale)
           if not self.current_event == self.event_to_process: 
              if self.current_event == self.scale - 1:
                 self.current_event = -1
              return
           else:
                if self.current_event == self.scale - 1:
                   self.current_event = -1

        self.events_run += 1
        hit_collection = {"pos" : [[], [], []], 
                          "d" : [[], [], []], 
                          "vert" : [], 
                          "index" : [],
                          "system" : [],
                          "detectorID" : [],
                          "B" : [[], [], []],
                          "time": [],
                          "mask": []}

        if ("us" in self.hits_to_fit) or ("ds" in self.hits_to_fit) or ("ve" in self.hits_to_fit) :
            # Loop through muon filter hits
            for i_hit, muFilterHit in enumerate(self.MuFilterHits) :
                # Don't use veto for fitting
                if muFilterHit.GetSystem() == 1 :
                    if "ve" not in self.hits_to_fit :
                        continue
                elif muFilterHit.GetSystem() == 2 :
                    if "us" not in self.hits_to_fit :
                        continue
                elif muFilterHit.GetSystem() == 3 :
                    if "ds" not in self.hits_to_fit :
                        continue
                else :
                    print("WARNING! Unknown MuFilter system!!")
            
                self.mufiDet.GetPosition(muFilterHit.GetDetectorID(), self.a, self.b)
            
                hit_collection["pos"][0].append(self.a.X())
                hit_collection["pos"][1].append(self.a.Y())
                hit_collection["pos"][2].append(self.a.Z())
            
                hit_collection["B"][0].append(self.b.X())
                hit_collection["B"][1].append(self.b.Y())
                hit_collection["B"][2].append(self.b.Z())
            
                hit_collection["vert"].append(muFilterHit.isVertical())
                hit_collection["system"].append(muFilterHit.GetSystem())
            
                hit_collection["d"][0].append(self.MuFilter_ds_dx)
                hit_collection["d"][2].append(self.MuFilter_ds_dz)
            
                hit_collection["index"].append(i_hit)
                
                hit_collection["detectorID"].append(muFilterHit.GetDetectorID())
                hit_collection["mask"].append(False)
            
                # Downstream
                if muFilterHit.GetSystem() == 3 :
                    hit_collection["d"][1].append(self.MuFilter_ds_dx)
                    if muFilterHit.isVertical(): 
                        if self.isMC:
                           hit_collection["time"].append(muFilterHit.GetAllTimes()[0]) #already in ns
                        else: hit_collection["time"].append(muFilterHit.GetAllTimes()[0]*6.25) #tdc2ns
                    else:
                        if self.isMC:
                           hit_collection["time"].append(-1) # FIXME - don't know what is best here
                        else: hit_collection["time"].append(muFilterHit.GetImpactT()) #already in ns
                # Upstream
                else :
                    hit_collection["d"][1].append(self.MuFilter_us_dy)
                    if self.isMC:
                           hit_collection["time"].append(-1) # FIXME - don't know what is best here
                    else: hit_collection["time"].append(muFilterHit.GetImpactT()) #already in ns
        
        if "sf" in self.hits_to_fit :
            if self.Scifi_meas:
               # Make scifi clusters
               self.clusScifi.Clear()
               self.scifiCluster()
               
               # Loop through scifi clusters
               for i_clust, scifiCl in enumerate(self.clusScifi) :
                   scifiCl.GetPosition(self.a,self.b)
                
                   hit_collection["pos"][0].append(self.a.X())
                   hit_collection["pos"][1].append(self.a.Y())
                   hit_collection["pos"][2].append(self.a.Z())
                
                   hit_collection["B"][0].append(self.b.X())
                   hit_collection["B"][1].append(self.b.Y())
                   hit_collection["B"][2].append(self.b.Z())
                   
                   # take the cluster size as the active area size
                   hit_collection["d"][0].append(scifiCl.GetN()*self.Scifi_dx)
                   hit_collection["d"][1].append(scifiCl.GetN()*self.Scifi_dy)
                   hit_collection["d"][2].append(self.Scifi_dz)
                
                   if int(scifiCl.GetFirst()/100000)%10==1:
                      hit_collection["vert"].append(True)
                   else: hit_collection["vert"].append(False)
                   hit_collection["index"].append(i_clust)
                
                   hit_collection["system"].append(0)
                   hit_collection["detectorID"].append(scifiCl.GetFirst())
                   hit_collection["mask"].append(False)
                   
                   if self.isMC : hit_collection["time"].append(scifiCl.GetTime()/6.25) # for MC, hit time is in ns. Then for MC Scifi cluster time one has to divide by tdc2ns
                   else: hit_collection["time"].append(scifiCl.GetTime()) # already in ns

            else:
                 if self.hits_for_triplet == 'sf' and self.hits_to_fit == 'sf':
                   # Loop through scifi hits and count hits per projection and plane
                   N_plane_ZY = {1:0, 2:0, 3:0, 4:0, 5:0}
                   N_plane_ZX = {1:0, 2:0, 3:0, 4:0, 5:0}
                   for scifiHit in self.ScifiHits:
                      if not scifiHit.isValid(): continue
                      if scifiHit.isVertical(): 
                         N_plane_ZX[scifiHit.GetStation()] += 1
                      else:
                         N_plane_ZY[scifiHit.GetStation()] += 1
                   mask_plane_ZY = []
                   mask_plane_ZX = []
                   # sorting
                   N_plane_ZY = dict(sorted(N_plane_ZY.items(), key=lambda item: item[1], reverse = True))
                   N_plane_ZX = dict(sorted(N_plane_ZX.items(), key=lambda item: item[1], reverse = True))
                   # count planes with hits
                   n_zx = self.Scifi_nPlanes - list(N_plane_ZX.values()).count(0)
                   n_zy = self.Scifi_nPlanes - list(N_plane_ZY.values()).count(0)
                   # always use the 2 least busy planes for Hough transform(HT).
                   # Addition to the HT pool of other more busy planes depends on number of hits per plane and hit ratio wrt the previous less-busy plane
                   # If a plane is to be masked, all other planes having higher hit occupancy are masked too.
                   for ii in range(n_zx-self.min_planes_hit, -1, -1):
                      if list(N_plane_ZX.values())[ii]/list(N_plane_ZX.values())[ii+1] > 3 and list(N_plane_ZX.values())[ii] > 10:
                         for iii in range(ii, -1, -1):
                           mask_plane_ZX.append(list(N_plane_ZX.keys())[iii])
                         break
                   for ii in range(n_zy-self.min_planes_hit, -1, -1):
                      if list(N_plane_ZY.values())[ii]/list(N_plane_ZY.values())[ii+1] > 3 and list(N_plane_ZY.values())[ii] > 10:
                         for iii in range(ii, -1, -1):
                           mask_plane_ZY.append(list(N_plane_ZY.keys())[iii]) 
                         break
                 
                 # Loop through scifi hits
                 for i_hit, scifiHit in enumerate(self.ScifiHits) :
                     if not scifiHit.isValid(): continue 
                     self.scifiDet.GetSiPMPosition(scifiHit.GetDetectorID(), self.a, self.b)
                     hit_collection["pos"][0].append(self.a.X())
                     hit_collection["pos"][1].append(self.a.Y())
                     hit_collection["pos"][2].append(self.a.Z())
            
                     hit_collection["B"][0].append(self.b.X())
                     hit_collection["B"][1].append(self.b.Y())
                     hit_collection["B"][2].append(self.b.Z())
            
                     hit_collection["d"][0].append(self.Scifi_dx)
                     hit_collection["d"][1].append(self.Scifi_dy)
                     hit_collection["d"][2].append(self.Scifi_dz)
                
                     hit_collection["vert"].append(scifiHit.isVertical())
                     hit_collection["index"].append(i_hit)
                
                     hit_collection["system"].append(0)
            
                     hit_collection["detectorID"].append(scifiHit.GetDetectorID())
                     
                     if self.hits_for_triplet == 'sf' and self.hits_to_fit == 'sf':
                       if (scifiHit.isVertical()==0 and scifiHit.GetStation() in mask_plane_ZY) or (scifiHit.isVertical() and scifiHit.GetStation() in mask_plane_ZX):
                          hit_collection["mask"].append(True)
                       else: hit_collection["mask"].append(False)
                     else:
                          hit_collection["mask"].append(False)

                     if self.isMC : hit_collection["time"].append(scifiHit.GetTime()) # already in ns
                     else:
                          hit_collection["time"].append(scifiHit.GetTime()*6.25) #tdc2ns
    
        # Make the hit collection numpy arrays.
        for key, item in hit_collection.items() :
            if key == 'vert' :
                this_dtype = np.bool
            elif key == 'mask' :
                this_dtype = np.bool
            elif key == "index" or key == "system" or key == "detectorID" :
                this_dtype = np.int32
            else :
                this_dtype = np.float
            hit_collection[key] = np.array(item, dtype = this_dtype)

        # Useful for later
        triplet_condition_system = []
        if "sf" in self.hits_for_triplet :
            triplet_condition_system.append(0)
        if "ve" in self.hits_for_triplet :
            triplet_condition_system.append(1)
        if "us" in self.hits_for_triplet :
            triplet_condition_system.append(2)
        if "ds" in self.hits_for_triplet :
            triplet_condition_system.append(3)

        # Reconstruct muons until there are not enough hits in downstream muon filter
        for i_muon in range(self.max_reco_muons) :
            triplet_hits_horizontal = np.logical_and( ~hit_collection["vert"],
                                                      np.isin(hit_collection["system"], triplet_condition_system) )
            triplet_hits_vertical = np.logical_and( hit_collection["vert"],
                                                    np.isin(hit_collection["system"], triplet_condition_system) )

            n_planes_ZY = numPlanesHit(hit_collection["system"][triplet_hits_horizontal],
                                       hit_collection["detectorID"][triplet_hits_horizontal])
            if n_planes_ZY < self.min_planes_hit :
                break

            n_planes_ZX = numPlanesHit(hit_collection["system"][triplet_hits_vertical],
                                       hit_collection["detectorID"][triplet_hits_vertical])
            if n_planes_ZX < self.min_planes_hit :
                break

            # Get hits in hough transform format
            muon_hits_horizontal = np.logical_and( np.logical_and( ~hit_collection["vert"], ~hit_collection["mask"]),
                                                   np.isin(hit_collection["system"], [1, 2, 3]))
            muon_hits_vertical = np.logical_and( np.logical_and( hit_collection["vert"], ~hit_collection["mask"]),
                                                 np.isin(hit_collection["system"], [1, 2, 3]))
            scifi_hits_horizontal = np.logical_and( np.logical_and( ~hit_collection["vert"], ~hit_collection["mask"]),
                                                    np.isin(hit_collection["system"], [0]))
            scifi_hits_vertical = np.logical_and( np.logical_and( hit_collection["vert"], ~hit_collection["mask"]),
                                                  np.isin(hit_collection["system"], [0]))


            ZY = np.dstack([np.concatenate([np.tile(hit_collection["pos"][2][muon_hits_horizontal], self.muon_weight),
                                            hit_collection["pos"][2][scifi_hits_horizontal]]),
                            np.concatenate([np.tile(hit_collection["pos"][1][muon_hits_horizontal], self.muon_weight),
                                            hit_collection["pos"][1][scifi_hits_horizontal]])])[0]

            d_ZY = np.dstack([np.concatenate([np.tile(hit_collection["d"][2][muon_hits_horizontal], self.muon_weight),
                                              hit_collection["d"][2][scifi_hits_horizontal]]),
                              np.concatenate([np.tile(hit_collection["d"][1][muon_hits_horizontal], self.muon_weight),
                                              hit_collection["d"][1][scifi_hits_horizontal]])])[0]

            ZX = np.dstack([np.concatenate([np.tile(hit_collection["pos"][2][muon_hits_vertical], self.muon_weight),
                                            hit_collection["pos"][2][scifi_hits_vertical]]),
                            np.concatenate([np.tile(hit_collection["pos"][0][muon_hits_vertical], self.muon_weight),
                                            hit_collection["pos"][0][scifi_hits_vertical]])])[0]

            d_ZX = np.dstack([np.concatenate([np.tile(hit_collection["d"][2][muon_hits_vertical], self.muon_weight),
                                              hit_collection["d"][2][scifi_hits_vertical]]),
                              np.concatenate([np.tile(hit_collection["d"][0][muon_hits_vertical], self.muon_weight),
                                              hit_collection["d"][0][scifi_hits_vertical]])])[0]

            ZY_hough = self.h_ZY.fit_randomize(ZY, d_ZY, self.n_random)
            ZX_hough = self.h_ZX.fit_randomize(ZX, d_ZX, self.n_random)

            tol = self.tolerance
            # Special treatment for events with low hit occupancy - increase tolerance
            # For Scifi-only tracks
            if len(hit_collection["detectorID"]) < 31 and self.hits_for_triplet == 'sf' and self.hits_to_fit == 'sf' :
               # as there are masked Scifi planes, make sure to use hit counts before the masking
               # make requirement on hits per plane consistent with plane mask one
               if max(N_plane_ZX.values()) < 5 and max(N_plane_ZY.values()) < 5:
                  tol = 5*self.tolerance
            # for DS-only tracks
            if len(hit_collection["detectorID"]) < 22 and self.hits_for_triplet == 'ds' and self.hits_to_fit == 'ds' :
               # Loop through hits and count hits per projection and plane
               N_plane_ZY = {0:0, 1:0, 2:0, 3:0}
               N_plane_ZX = {0:0, 1:0, 2:0, 3:0}
               for item in range(len(hit_collection["detectorID"])):
                   if hit_collection["vert"][item]: N_plane_ZX[(hit_collection["detectorID"][item]%10000)//1000] += 1
                   else: N_plane_ZY[(hit_collection["detectorID"][item]%10000)//1000] += 1
               if max(N_plane_ZX.values()) < 4 and max(N_plane_ZY.values()) < 4:
                  tol = 3*self.tolerance

            # Check if track intersects minimum number of hits in each plane.
            track_hits_for_triplet_ZY = hit_finder(ZY_hough[0], ZY_hough[1], 
                                                   np.dstack([hit_collection["pos"][2][triplet_hits_horizontal],
                                                              hit_collection["pos"][1][triplet_hits_horizontal]]),
                                                   np.dstack([hit_collection["d"][2][triplet_hits_horizontal],
                                                              hit_collection["d"][1][triplet_hits_horizontal]]), tol)

            track_hits_for_triplet_ZX = hit_finder(ZX_hough[0], ZX_hough[1], 
                                                   np.dstack([hit_collection["pos"][2][triplet_hits_vertical],
                                                              hit_collection["pos"][0][triplet_hits_vertical]]),
                                                   np.dstack([hit_collection["d"][2][triplet_hits_vertical],
                                                              hit_collection["d"][0][triplet_hits_vertical]]), tol)
                                                   
            n_planes_hit_ZY = numPlanesHit(hit_collection["system"][triplet_hits_horizontal][track_hits_for_triplet_ZY],
                                           hit_collection["detectorID"][triplet_hits_horizontal][track_hits_for_triplet_ZY])
            if n_planes_hit_ZY < self.min_planes_hit :
                break
            n_planes_hit_ZX = numPlanesHit(hit_collection["system"][triplet_hits_vertical][track_hits_for_triplet_ZX],
                                           hit_collection["detectorID"][triplet_hits_vertical][track_hits_for_triplet_ZX])
            if n_planes_hit_ZX < self.min_planes_hit :
                break

#                print("Found {0} downstream ZX planes associated to muon track".format(n_planes_ds_ZX))
#                print("Found {0} downstream ZY planes associated to muon track".format(n_planes_ds_ZY))
            
            # This time with all the hits, not just triplet condition.
            track_hits_ZY = hit_finder(ZY_hough[0], ZY_hough[1], 
                                       np.dstack([hit_collection["pos"][2][~hit_collection["vert"]], 
                                                  hit_collection["pos"][1][~hit_collection["vert"]]]), 
                                       np.dstack([hit_collection["d"][2][~hit_collection["vert"]],
                                                  hit_collection["d"][1][~hit_collection["vert"]]]), tol)

            track_hits_ZX = hit_finder(ZX_hough[0], ZX_hough[1], 
                                       np.dstack([hit_collection["pos"][2][hit_collection["vert"]], 
                                                  hit_collection["pos"][0][hit_collection["vert"]]]), 
                                       np.dstack([hit_collection["d"][2][hit_collection["vert"]], 
                                                  hit_collection["d"][0][hit_collection["vert"]]]), tol)
            # Onto Kalman fitter (based on SndlhcTracking.py)
            posM    = ROOT.TVector3(0, 0, 0.)
            momM = ROOT.TVector3(0,0,100.)  # default track with high momentum

            # approximate covariance
            covM = ROOT.TMatrixDSym(6)
            res = self.kalman_sigmaScifi_spatial
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
            
            # Sort measurements in Z
            hit_z = np.concatenate([hit_collection["pos"][2][hit_collection["vert"]][track_hits_ZX],
                                    hit_collection["pos"][2][~hit_collection["vert"]][track_hits_ZY]])

            hit_A0 = np.concatenate([hit_collection["pos"][0][hit_collection["vert"]][track_hits_ZX],
                                     hit_collection["pos"][0][~hit_collection["vert"]][track_hits_ZY]])

            hit_A1 = np.concatenate([hit_collection["pos"][1][hit_collection["vert"]][track_hits_ZX],
                                     hit_collection["pos"][1][~hit_collection["vert"]][track_hits_ZY]])
            
            hit_B0 = np.concatenate([hit_collection["B"][0][hit_collection["vert"]][track_hits_ZX],
                                     hit_collection["B"][0][~hit_collection["vert"]][track_hits_ZY]])

            hit_B1 = np.concatenate([hit_collection["B"][1][hit_collection["vert"]][track_hits_ZX],
                                     hit_collection["B"][1][~hit_collection["vert"]][track_hits_ZY]])
            
            hit_B2 = np.concatenate([hit_collection["B"][2][hit_collection["vert"]][track_hits_ZX],
                                     hit_collection["B"][2][~hit_collection["vert"]][track_hits_ZY]])
            
            hit_detid = np.concatenate([hit_collection["detectorID"][hit_collection["vert"]][track_hits_ZX],
                                        hit_collection["detectorID"][~hit_collection["vert"]][track_hits_ZY]])

            kalman_spatial_sigma = np.concatenate([hit_collection["d"][0][hit_collection["vert"]][track_hits_ZX] / 12**0.5,
                                                   hit_collection["d"][1][~hit_collection["vert"]][track_hits_ZY] / 12**0.5])

            # Maximum distance. Use (d_xy/2**2 + d_z/2**2)**0.5
            kalman_max_dis = np.concatenate([((hit_collection["d"][0][hit_collection["vert"]][track_hits_ZX]/2.)**2 +
                                              (hit_collection["d"][2][hit_collection["vert"]][track_hits_ZX]/2.)**2)**0.5,
                                             ((hit_collection["d"][1][~hit_collection["vert"]][track_hits_ZY]/2.)**2 +
                                              (hit_collection["d"][2][~hit_collection["vert"]][track_hits_ZY]/2.)**2)**0.5])
            
            hitID = 0 # Does it matter? We don't have a global hit ID.

            hit_time = np.concatenate([hit_collection["time"][hit_collection["vert"]][track_hits_ZX],
                                      hit_collection["time"][~hit_collection["vert"]][track_hits_ZY]])

            for i_z_sorted in hit_z.argsort() :
                tp = ROOT.genfit.TrackPoint()
                hitCov = ROOT.TMatrixDSym(7)
                hitCov[6][6] = kalman_spatial_sigma[i_z_sorted]**2
                
                measurement = ROOT.genfit.WireMeasurement(ROOT.TVectorD(7, array('d', [hit_A0[i_z_sorted],
                                                                                       hit_A1[i_z_sorted],
                                                                                       hit_z[i_z_sorted],
                                                                                       hit_B0[i_z_sorted],
                                                                                       hit_B1[i_z_sorted],
                                                                                       hit_B2[i_z_sorted],
                                                                                       0.])),
                                                          hitCov,
                                                          1, # detid?
                                                          6, # hitid?
                                                          tp)
                
                measurement.setMaxDistance(kalman_max_dis[i_z_sorted])
                measurement.setDetId(int(hit_detid[i_z_sorted]))
                measurement.setHitId(int(hitID))
                hitID += 1
                tp.addRawMeasurement(measurement)
                theTrack.insertPoint(tp)

            if not theTrack.checkConsistency():
                theTrack.Delete()
                raise RuntimeException("Kalman fitter track consistency check failed.")

            # do the fit
            self.kalman_fitter.processTrack(theTrack) # processTrackWithRep(theTrack,rep,True)

            fitStatus = theTrack.getFitStatus()
            if not fitStatus.isFitConverged() and 0>1:
                theTrack.Delete()
                raise RuntimeException("Kalman fit did not converge.")
            
            # Now save the track!
            if self.genfitTrack: self.kalman_tracks.Add(theTrack)
            else :            
                # Load items into snd track class object
                this_track = ROOT.sndRecoTrack(theTrack)
                pointTimes = []
                for i_z_sorted in hit_z.argsort() :
                    pointTimes.append(hit_time[i_z_sorted])
                this_track.setRawMeasTimes(pointTimes)
                this_track.setTrackType(self.track_type)
                # Save the track in sndRecoTrack format
                self.kalman_tracks[i_muon] = this_track
            

            # Remove track hits and try to find an additional track
            # Find array index to be removed
            index_to_remove_ZX = np.where(np.in1d(hit_collection["detectorID"], hit_collection["detectorID"][hit_collection["vert"]][track_hits_ZX]))[0]
            index_to_remove_ZY = np.where(np.in1d(hit_collection["detectorID"], hit_collection["detectorID"][~hit_collection["vert"]][track_hits_ZY]))[0]

            index_to_remove = np.concatenate([index_to_remove_ZX, index_to_remove_ZY])
            
            # Remove dictionary entries 
            for key in hit_collection.keys() :
                if len(hit_collection[key].shape) == 1 :
                    hit_collection[key] = np.delete(hit_collection[key], index_to_remove)
                elif len(hit_collection[key].shape) == 2 :
                    hit_collection[key] = np.delete(hit_collection[key], index_to_remove, axis = 1)
                else :
                    raise Exception("Wrong number of dimensions found when deleting hits in iterative muon identification algorithm.")

    def FinishTask(self) :
        print("Processed" ,self.events_run)
        if not self.genfitTrack : self.kalman_tracks.Delete()
        else : pass

    # this is a copy on SndlhcTracking function with small adjustments in event object names to make it work here
    # FIXME Should find a way to use this function straight from the SndlhcTracking!
    def scifiCluster(self):
       clusters = []
       hitDict = {}
       for k in range(self.ScifiHits.GetEntries()):
            d = self.ScifiHits[k]
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
                        for aHit in tmp: hitvector.push_back( self.ScifiHits[hitDict[aHit]])
                        aCluster = ROOT.sndCluster(first,N,hitvector,self.scifiDet,False)
                        clusters.append(aCluster)
                        if c!=hitList[last]:
                             ncl+=1
                             tmp = [c]
                        elif not neighbour :   # save last channel
                            hitvector.clear()
                            hitvector.push_back(self.ScifiHits[hitDict[c]])
                            aCluster = ROOT.sndCluster(c,1,hitvector,self.scifiDet,False)
                            clusters.append(aCluster)
                   cprev = c
       self.clusScifi.Delete()
       for c in clusters:  
            self.clusScifi.Add(c)
