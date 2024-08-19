"""Perform Hough Transform pattern recognition and a track fit."""

import ROOT
import numpy as np
import scipy.ndimage
from array import array
import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt
import shipunit as unit

ROOT.gInterpreter.ProcessLine('#include "SiSensor.h"')


def hit_finder(slope, intercept, box_centers, box_ds, tol=0.0):
    """Find hits intersected by Hough line."""
    # First check if track at center of box is within box limits
    d = np.abs(box_centers[0, :, 1] - (box_centers[0, :, 0] * slope + intercept))
    center_in_box = d < (box_ds[0, :, 1] + tol) / 2.0

    # Now check if, assuming line is not in box at box center, the slope is large enough for line to clip the box at corner
    clips_corner = np.abs(slope) > np.abs(
        (d - (box_ds[0, :, 1] + tol) / 2.0) / (box_ds[0, :, 0] + tol) / 2.0
    )

    # If either of these is true, line goes through hit:
    hit_mask = np.logical_or(center_in_box, clips_corner)

    # Return indices
    return np.where(hit_mask)[0]


class hough:
    """Hough transform implementation."""

    def __init__(
        self,
        n_yH,
        yH_range,
        n_xH,
        xH_range,
        z_offset,
        Hformat,
        space_scale,
        det_Zlen,
        squaretheta=False,
        smooth=True,
    ):
        """Initialize."""
        self.n_yH = n_yH
        self.n_xH = n_xH

        self.yH_range = yH_range
        self.xH_range = xH_range

        self.z_offset = z_offset
        self.HoughSpace_format = Hformat
        self.space_scale = space_scale

        self.det_Zlen = det_Zlen

        self.smooth = smooth

        self.yH_bins = np.linspace(self.yH_range[0], self.yH_range[1], n_yH)
        if not squaretheta:
            self.xH_bins = np.linspace(self.xH_range[0], self.xH_range[1], n_xH)
        else:
            self.xH_bins = np.linspace(
                np.sign(self.xH_range[0]) * (self.xH_range[0] ** 0.5),
                np.sign(self.xH_range[1]) * (self.xH_range[1] ** 0.5),
                n_xH,
            )
            self.xH_bins = np.sign(self.xH_bins) * np.square(self.xH_bins)

        self.cos_thetas = np.cos(self.xH_bins)
        self.sin_thetas = np.sin(self.xH_bins)

        self.xH_i = np.array(list(range(n_xH)))

        # A back-up Hough space designed to have more/less bins wrt the default one above.
        # It is useful when fitting some low-E muon tracks, which are curved due to mult. scattering.
        self.n_yH_scaled = int(n_yH * space_scale)
        self.n_xH_scaled = int(n_xH * space_scale)
        self.yH_bins_scaled = np.linspace(
            self.yH_range[0], self.yH_range[1], self.n_yH_scaled
        )
        if not squaretheta:
            self.xH_bins_scaled = np.linspace(
                self.xH_range[0], self.xH_range[1], self.n_xH_scaled
            )
        else:
            self.xH_bins_scaled = np.linspace(
                np.sign(self.xH_range[0]) * (self.xH_range[0] ** 0.5),
                np.sign(self.xH_range[1]) * (self.xH_range[1] ** 0.5),
                self.n_xH_scaled,
            )
            self.xH_bins_scaled = np.sign(self.xH_bins_scaled) * np.square(
                self.xH_bins_scaled
            )

        self.cos_thetas_scaled = np.cos(self.xH_bins_scaled)
        self.sin_thetas_scaled = np.sin(self.xH_bins_scaled)

        self.xH_i_scaled = np.array(list(range(self.n_xH_scaled)))

    def fit(self, hit_collection, is_scaled, draw, weights=None):
        """Find Hough Transform space maxima."""
        if not is_scaled:
            n_xH = self.n_xH
            n_yH = self.n_yH
            cos_thetas = self.cos_thetas
            sin_thetas = self.sin_thetas
            xH_bins = self.xH_bins
            yH_bins = self.yH_bins
            xH_i = self.xH_i
            res = self.res
        else:
            n_xH = self.n_xH_scaled
            n_yH = self.n_yH_scaled
            cos_thetas = self.cos_thetas_scaled
            sin_thetas = self.sin_thetas_scaled
            xH_bins = self.xH_bins_scaled
            yH_bins = self.yH_bins_scaled
            xH_i = self.xH_i_scaled
            res = self.res * self.space_scale

        self.accumulator = np.zeros((n_yH, n_xH))
        for i_hit, hit in enumerate(hit_collection):
            shifted_hitZ = hit[0] - self.z_offset
            if self.HoughSpace_format == "normal":
                hit_yH = shifted_hitZ * cos_thetas + hit[1] * sin_thetas
            elif self.HoughSpace_format == "linearSlopeIntercept":
                hit_yH = hit[1] - shifted_hitZ * xH_bins
            elif self.HoughSpace_format == "linearIntercepts":
                hit_yH = (self.det_Zlen * hit[1] - shifted_hitZ * xH_bins) / (
                    self.det_Zlen - shifted_hitZ
                )
            out_of_range = np.logical_and(
                hit_yH > self.yH_range[0], hit_yH < self.yH_range[1]
            )
            hit_yH_i = np.floor(
                (hit_yH[out_of_range] - self.yH_range[0])
                / (self.yH_range[1] - self.yH_range[0])
                * n_yH
            ).astype(np.int_)

            if weights is not None:
                self.accumulator[hit_yH_i, xH_i[out_of_range]] += weights[i_hit]
            else:
                self.accumulator[hit_yH_i, xH_i[out_of_range]] += 1

        # Smooth accumulator
        if self.smooth_full:
            self.accumulator = scipy.ndimage.gaussian_filter(
                self.accumulator, self.sigma, truncate=self.truncate
            )

        # This might be useful for debugging, but leave out for now.
        if draw:
            plt.figure()
            plt.imshow(
                self.accumulator,
                origin="lower",
                extent=[
                    self.xH_range[0],
                    self.xH_range[-1],
                    self.yH_range[0],
                    self.yH_range[-1],
                ],
                aspect="auto",
            )
            if self.HoughSpace_format == "normal":
                plt.xlabel(r"$\theta$ [rad]")
                plt.ylabel("r [cm]")
            elif self.HoughSpace_format == "linearSlopeIntercept":
                plt.xlabel("slope")
                plt.ylabel("intercept @ 1st plane [cm]")
            elif self.HoughSpace_format == "linearIntercepts":
                plt.xlabel("intercept @ last plane [cm]")
                plt.ylabel("intercept @ 1st plane [cm]")
            plt.tight_layout()
            plt.show()

        if self.smooth_full:
            # In case of multiple occurrences of the maximum values, argmax returns
            # the indices corresponding to the first occurrence(along 1st axis).
            i_max = np.unravel_index(self.accumulator.argmax(), self.accumulator.shape)
        else:
            # In case there are more than 1 bins with the maximal Nentries, check if the n-th quantile of
            # found peaks along 'slope' axis(yH axis) enloses <n-th quantile> portion of all maxima 'slope' bins
            # within a specified range. It is advisable that that range is consistent with
            # the detector angular resolution.
            maxima = np.argwhere(self.accumulator == np.amax(self.accumulator))
            if len(maxima) == 1:
                i_max = maxima[0]
            elif len(maxima) == 0 or (
                len(maxima) > 1 and self.HoughSpace_format == "linearIntercepts"
            ):
                # if no reasonable way to select btw maxima, force track-build failure
                # to be decided what to do in 'linearIntercepts' case and multiple maxima
                return (-999, -999)
            elif (
                len(maxima) > 1
                and abs(min([x[1] for x in maxima]) - max([x[1] for x in maxima])) < res
            ):
                i_max = maxima[0]
            else:
                # FIXME: the next two lines can be a single line command for sure
                maxima_slopesAxis_list = list([x[1] for x in maxima])
                maxima_slopesAxis_list = np.asarray(maxima_slopesAxis_list)
                quantile = np.quantile(maxima_slopesAxis_list, self.n_quantile)
                Nwithin = 0
                # FIXME: a more elegant 'hidden' loop is maybe possible here too
                for item in maxima:
                    if abs(item[1] - quantile) < res:
                        Nwithin += 1
                if Nwithin / len(maxima) > self.n_quantile:
                    i_x = min([x[1] for x in maxima], key=lambda b: abs(b - quantile))
                    for im in maxima:
                        if im[1] == i_x:
                            i_max = im
                else:
                    # if no reasonable way to select btw maxima, force track-build failure
                    return (-999, -999)

        found_yH = yH_bins[int(i_max[0])]
        found_xH = xH_bins[int(i_max[1])]

        if self.HoughSpace_format == "normal":
            slope = -1.0 / np.tan(found_xH)
            interceptShift = found_yH / np.sin(found_xH)
            intercept = (np.tan(found_xH) * interceptShift + self.z_offset) / np.tan(
                found_xH
            )
        elif self.HoughSpace_format == "linearSlopeIntercept":
            slope = found_xH
            intercept = found_yH - slope * self.z_offset
        elif self.HoughSpace_format == "linearIntercepts":
            slope = (found_xH - found_yH) / self.det_Zlen
            intercept = found_yH - slope * self.z_offset

        return (slope, intercept)

    def fit_randomize(
        self, hit_collection, hit_d, n_random, is_scaled, draw, weights=None
    ):
        """Do n_random throws per measurement."""
        if not len(hit_collection):
            return (-1, -1, [[], []], [], False)

        # Randomize hits
        if n_random > 0:
            random_hit_collection = []
            for i_random in range(n_random):
                random_hits = np.random.uniform(size=hit_collection.shape) - 0.5
                random_hits *= hit_d
                random_hits += hit_collection
                random_hit_collection.append(random_hits)

            random_hit_collection = np.concatenate(random_hit_collection)
            if weights is not None:
                weights = np.tile(weights, n_random)

            fit = self.fit(random_hit_collection, is_scaled, draw, weights)
        else:
            fit = self.fit(hit_collection, is_scaled, draw, weights)

        return fit


def numPlanesHit(detector_ids):
    """Return how many planes have a hit."""
    adv_target_stations = []

    adv_target_stations.append(detector_ids >> 17)

    return len(np.unique(adv_target_stations))


class MuonReco(ROOT.FairTask):
    """Muon reconstruction."""

    def Init(self):
        """Initialize the task."""
        self.logger = ROOT.FairLogger.GetLogger()
        if self.logger.IsLogNeeded(ROOT.fair.Severity.info):
            print("Initializing muon reconstruction task!")

        self.lsOfGlobals = ROOT.gROOT.GetListOfGlobals()
        self.AdvTargetDet = self.lsOfGlobals.FindObject("AdvTarget")
        self.ioman = ROOT.FairRootManager.Instance()

        # MC or data - needed for hit timing unit
        if self.ioman.GetInTree().GetName() == "rawConv":
            self.isMC = False
        else:
            self.isMC = True

        # Fetch digi hit collections from online if exist
        sink = self.ioman.GetSink()
        eventTree = None
        if sink:
            eventTree = sink.GetOutTree()
        if eventTree:
            self.AdvTargetHits = eventTree.Digi_AdvTargetHits
            self.EventHeader = eventTree.EventHeader
        else:
            # Try the FairRoot way
            self.AdvTargetHits = self.ioman.GetObject("Digi_AdvTargetHits")
            self.EventHeader = self.ioman.GetObject("EventHeader")

            # If that doesn't work, try using standard ROOT
            if self.AdvTargetHits is None:
                if self.logger.IsLogNeeded(ROOT.fair.Severity.info):
                    print("Digi_AdvTargetHits not in branch list")
                self.AdvTargetHits = self.ioman.GetInTree().Digi_AdvTargetHits
            if self.EventHeader is None:
                if self.logger.IsLogNeeded(ROOT.fair.Severity.info):
                    print("EventHeader not in branch list")
                self.EventHeader = self.ioman.GetInTree().EventHeader

        if self.AdvTargetHits is None:
            raise RuntimeError("Digi_AdvTargetHits not found in input file.")
        if self.EventHeader is None:
            raise RuntimeError("EventHeader not found in input file.")

        # Initialize event counters in case scaling of events is required
        self.scale = 1
        self.events_run = 0

        # Initialize hough transform - reading parameter xml file
        tree = ET.parse(self.par_file)
        root = tree.getroot()

        # Output track in genfit::Track or sndRecoTrack format
        # Check if genfit::Track format is already forced
        if hasattr(self, "genfitTrack"):
            pass
        else:
            self.genfitTrack = int(root[0].text)

        self.draw = int(root[1].text)

        track_case_exists = False
        for case in root.findall("tracking_case"):
            if case.get("name") == self.tracking_case:
                track_case_exists = True
                # Use AdvTarget hits or clusters
                self.AdvTarget_meas = int(case.find("use_Scifi_clust").text)
                # Maximum absolute value of reconstructed angle (+/- 1 rad is the maximum angle to form a triplet in the SciFi)
                max_angle = float(case.find("max_angle").text)

                # Hough space representation
                Hspace_format_exists = False
                for rep in case.findall("Hough_space_format"):
                    if rep.get("name") == self.Hough_space_format:
                        Hspace_format_exists = True
                        # Number of bins per Hough accumulator axes and range
                        """ xH and yH are the abscissa and ordinate of the Hough parameter space
                          xz and yz represent horizontal and vertical projections
                          in the SNDLHC physics coord. system """
                        n_accumulator_yH = int(rep.find("N_yH_bins").text)
                        yH_min_xz = float(rep.find("yH_min_xz").text)
                        yH_max_xz = float(rep.find("yH_max_xz").text)
                        yH_min_yz = float(rep.find("yH_min_yz").text)
                        yH_max_yz = float(rep.find("yH_max_yz").text)
                        n_accumulator_xH = int(rep.find("N_xH_bins").text)
                        xH_min_xz = float(rep.find("xH_min_xz").text)
                        xH_max_xz = float(rep.find("xH_max_xz").text)
                        xH_min_yz = float(rep.find("xH_min_yz").text)
                        xH_max_yz = float(rep.find("xH_max_yz").text)

                    else:
                        continue
                if not Hspace_format_exists:
                    raise RuntimeError(
                        "Unknown Hough space format, check naming in parameter xml file."
                    )

                # A scale factor for a back-up Hough space having more/less bins than the default one
                # It is useful when fitting some low-E muon tracks, which are curved due to mult. scattering.
                self.HT_space_scale = (
                    1
                    if case.find("HT_space_scale") is None
                    else float(case.find("HT_space_scale").text)
                )

                # Number of random throws per hit
                self.n_random = int(case.find("n_random").text)
                # Minimum number of planes hit in each of the downstream muon filter (if muon filter hits used) or scifi (if muon filter hits not used) views to try to reconstruct a muon
                self.min_planes_hit = int(case.find("min_planes_hit").text)

                # Maximum number of muons to find. To avoid spending too much time on events with lots of downstream activity.
                self.max_reco_muons = int(case.find("max_reco_muons").text)

                # How far away from Hough line hits will be assigned to the muon, for Kalman tracking
                self.tolerance = float(case.find("tolerance").text)

                # Detector plane masking. If flag is active, a plane will be masked if its N_hits > Nhits_per_plane.
                # In any case, plane masking will only be applied if solely Scifi hits are used in HT as it is
                # a measure against having many maxima in HT space.
                self.mask_plane = int(case.find("mask_plane").text)
                self.Nhits_per_plane = int(case.find("Nhits_per_plane").text)

                # Enable Gaussian smoothing over the full accumulator space.
                self.smooth_full = int(case.find("smooth_full").text)
                # Gaussian smoothing parameters. The kernel size is determined as 2*int(truncate*sigma+0.5)+1
                self.sigma = int(case.find("sigma").text)
                self.truncate = int(case.find("truncate").text)
                # Helpers to pick up one of many HT space maxima
                self.n_quantile = float(case.find("n_quantile").text)
                self.res = int(case.find("res").text)

            else:
                continue
        if not track_case_exists:
            raise RuntimeError(
                "Unknown tracking case, check naming in parameter xml file."
            )

        # Get sensor dimensions from geometry
        self.AdvTarget_dx = ROOT.advsnd.sensor_width / ROOT.advsnd.strips
        self.AdvTarget_dy = self.AdvTarget_dx
        self.AdvTarget_dz = self.AdvTargetDet.GetConfParI("AdvTarget/TTZ")

        # Get number of readout channels
        self.AdvTarget_nPlanes = self.AdvTargetDet.GetConfParI("AdvTarget/nTT")
        self.max_n_hits_plane = 3
        self.max_n_AdvTarget_hits = self.max_n_hits_plane * 2 * self.AdvTarget_nPlanes

        # get the distance between 1st and last detector planes to be used in the track fit.
        # a z_offset is used to shift detector hits so to have smaller Hough parameter space
        # Using geometers measurements! For safety, add a 5-cm-buffer in detector lengths and a 2.5-cm one to z_offset.
        # This is done to account for possible det. position shifts/mismatches going from geom. measurements and sndsw physics CS.
        # FIXME get these form the geometry
        det_Zlen = 200 * unit.cm
        z_offset = -125 * unit.cm

        # Initialize Hough transforms for both views:
        if self.Hough_space_format == "normal":
            # rho-theta representation - must exclude theta range for tracks passing < 3 det. planes
            self.h_ZX = hough(
                n_accumulator_yH,
                [yH_min_xz, yH_max_xz],
                n_accumulator_xH,
                [-max_angle + np.pi / 2.0, max_angle + np.pi / 2.0],
                z_offset,
                self.Hough_space_format,
                self.HT_space_scale,
                det_Zlen,
            )
            self.h_ZY = hough(
                n_accumulator_yH,
                [yH_min_yz, yH_max_yz],
                n_accumulator_xH,
                [-max_angle + np.pi / 2.0, max_angle + np.pi / 2.0],
                z_offset,
                self.Hough_space_format,
                self.HT_space_scale,
                det_Zlen,
            )
        else:
            self.h_ZX = hough(
                n_accumulator_yH,
                [yH_min_xz, yH_max_xz],
                n_accumulator_xH,
                [xH_min_xz, xH_max_xz],
                z_offset,
                self.Hough_space_format,
                self.HT_space_scale,
                det_Zlen,
            )
            self.h_ZY = hough(
                n_accumulator_yH,
                [yH_min_yz, yH_max_yz],
                n_accumulator_xH,
                [xH_min_yz, xH_max_yz],
                z_offset,
                self.Hough_space_format,
                self.HT_space_scale,
                det_Zlen,
            )

        self.h_ZX.smooth_full = self.smooth_full
        self.h_ZY.smooth_full = self.smooth_full
        self.h_ZX.sigma = self.sigma
        self.h_ZX.truncate = self.truncate
        self.h_ZY.sigma = self.sigma
        self.h_ZY.truncate = self.truncate

        self.h_ZX.n_quantile = self.n_quantile
        self.h_ZX.res = self.res
        self.h_ZY.n_quantile = self.n_quantile
        self.h_ZY.res = self.res

        self.track_type = 11

        # To keep temporary detector information
        self.a = ROOT.TVector3()
        self.b = ROOT.TVector3()

        # check if track container exists
        if self.ioman.GetObject("Reco_MuonTracks") != None:
            self.kalman_tracks = self.ioman.GetObject("Reco_MuonTracks")
            if self.logger.IsLogNeeded(ROOT.fair.Severity.info):
                print("Branch activated by another task!")
        else:
            # Now initialize output in genfit::track or sndRecoTrack format
            if self.genfitTrack:
                self.kalman_tracks = ROOT.TObjArray(50)
                if hasattr(self, "standalone") and self.standalone:
                    self.ioman.Register(
                        "Reco_MuonTracks",
                        self.ioman.GetFolderName(),
                        self.kalman_tracks,
                        ROOT.kTRUE,
                    )
            else:
                self.kalman_tracks = ROOT.TClonesArray("sndRecoTrack", 50)
                if hasattr(self, "standalone") and self.standalone:
                    self.ioman.Register(
                        "Reco_MuonTracks", "", self.kalman_tracks, ROOT.kTRUE
                    )

        # initialize detector class with EventHeader(runN), if SNDLHCEventHeader detected
        # only needed if using HT tracking manager, i.e. standalone
        '''
        if (
            self.EventHeader.IsA().GetName() == "SNDLHCEventHeader"
            and hasattr(self, "standalone")
            and self.standalone
            and not self.isMC
        ):
            self.AdvTargetDet.InitEvent(self.EventHeader)
        '''

        # internal storage of clusters
        if self.AdvTarget_meas:
            self.clusScifi = ROOT.TObjArray(100)

        # Kalman filter stuff

        geoMat = ROOT.genfit.TGeoMaterialInterface()
        bfield = ROOT.genfit.ConstField(0, 0, 0)
        fM = ROOT.genfit.FieldManager.getInstance()
        fM.init(bfield)
        ROOT.genfit.MaterialEffects.getInstance().init(geoMat)
        ROOT.genfit.MaterialEffects.getInstance().setNoEffects()

        # self.kalman_fitter = ROOT.genfit.KalmanFitter()
        # self.kalman_fitter.setMaxIterations(50)
        # KalmanFitter(unsigned int maxIterations = 4, double deltaPval = 1e-3, double blowUpFactor = 1e3, bool squareRootFormalism = false)
        # DAF (bool useRefKalman=true, double deltaPval=1e-3, double deltaWeight=1e-3)
        self.kalman_fitter = ROOT.genfit.DAF(True, 1e-3, 1e-3)
        self.kalman_sigmaAdvTarget_spatial = self.AdvTarget_dx / 12**0.5

        # Init() MUST return int
        return 0

    def SetScaleFactor(self, scale):
        """Set the scale factor for the Hough space."""
        self.scale = scale

    def SetParFile(self, file_name):
        """Set the parameter file."""
        self.par_file = file_name

    def SetTrackingCase(self, case):
        """Set the tracking case."""
        self.tracking_case = case

    def SetHoughSpaceFormat(self, Hspace_format):
        """Set the format of the Hough space."""
        self.Hough_space_format = Hspace_format

    def ForceGenfitTrackFormat(self):
        """Set the track format."""
        self.genfitTrack = 1

    # flag showing the task is run seperately from other tracking tasks
    def SetStandalone(self):
        """See if the task is run alone."""
        self.standalone = 1

    def Exec(self, opt):
        """Core part of teh task."""
        self.kalman_tracks.Clear("C")

        # Set scaling in case task is run seperately from other tracking tasks
        if self.scale > 1 and self.standalone:
            if ROOT.gRandom.Rndm() > 1.0 / self.scale:
                return

        self.events_run += 1
        hit_collection = {
            "pos": [[], [], []],
            "d": [[], [], []],
            "vert": [],
            "index": [],
            "detectorID": [],
            "B": [[], [], []],
            "time": [],
            "mask": [],
        }

        if 1:
            # if self.AdvTarget_meas:
            # no AdvTarget clustering available for now
            # baraly any adaptation from SND to AdvSND done for this case
            if 0:
                # Make AdvTarget clusters
                self.clusScifi.Clear()
                self.scifiCluster()

                # Loop through scifi clusters
                for i_clust, scifiCl in enumerate(self.clusScifi):
                    scifiCl.GetPosition(self.a, self.b)

                    hit_collection["pos"][0].append(self.a.X())
                    hit_collection["pos"][1].append(self.a.Y())
                    hit_collection["pos"][2].append(self.a.Z())

                    hit_collection["B"][0].append(self.b.X())
                    hit_collection["B"][1].append(self.b.Y())
                    hit_collection["B"][2].append(self.b.Z())

                    # take the cluster size as the active area size
                    hit_collection["d"][0].append(scifiCl.GetN() * self.AdvTarget_dx)
                    hit_collection["d"][1].append(scifiCl.GetN() * self.AdvTarget_dy)
                    hit_collection["d"][2].append(self.AdvTarget_dz)

                    if int(scifiCl.GetFirst() / 100000) % 10 == 1:
                        hit_collection["vert"].append(True)
                    else:
                        hit_collection["vert"].append(False)
                    hit_collection["index"].append(i_clust)

                    hit_collection["detectorID"].append(scifiCl.GetFirst())
                    hit_collection["mask"].append(False)

                    times = []
                    if self.isMC:
                        times.append(
                            scifiCl.GetTime() / 6.25
                        )  # for MC, hit time is in ns. Then for MC Scifi cluster time one has to divide by tdc2ns
                    else:
                        times.append(scifiCl.GetTime())  # already in ns
                    hit_collection["time"].append(times)

            else:
                # Loop through AdvTarget hits and count hits per projection and plane
                N_plane_ZY = {}
                N_plane_ZX = {}
                for i in range(self.AdvTarget_nPlanes):
                    N_plane_ZY[i] = 0
                    N_plane_ZX[i] = 0
                for adv_target_hit in self.AdvTargetHits:
                    # if not adv_target_hit.isValid(): continue
                    if adv_target_hit.isVertical():
                        N_plane_ZX[adv_target_hit.GetStation()] += 1
                    else:
                        N_plane_ZY[adv_target_hit.GetStation()] += 1
                if self.mask_plane:
                    mask_plane_ZY = []
                    mask_plane_ZX = []
                    # sorting
                    N_plane_ZY = dict(
                        sorted(
                            N_plane_ZY.items(), key=lambda item: item[1], reverse=True
                        )
                    )
                    N_plane_ZX = dict(
                        sorted(
                            N_plane_ZX.items(), key=lambda item: item[1], reverse=True
                        )
                    )
                    # count planes with hits
                    n_zx = self.AdvTarget_nPlanes - list(N_plane_ZX.values()).count(0)
                    n_zy = self.AdvTarget_nPlanes - list(N_plane_ZY.values()).count(0)
                    # check with min number of hit planes
                    if n_zx < self.min_planes_hit or n_zy < self.min_planes_hit:
                        return
                    # mask busiest planes until there are at least 3 planes with hits left
                    for ii in range(n_zx - self.min_planes_hit):
                        if list(N_plane_ZX.values())[ii] > self.Nhits_per_plane:
                            mask_plane_ZX.append(list(N_plane_ZX.keys())[ii])
                    for ii in range(n_zy - self.min_planes_hit):
                        if list(N_plane_ZY.values())[ii] > self.Nhits_per_plane:
                            mask_plane_ZY.append(list(N_plane_ZY.keys())[ii])

                # Loop through scifi hits
                for i_hit, adv_target_hit in enumerate(self.AdvTargetHits):
                    # if not adv_target_hit.isValid(): continue
                    self.AdvTargetDet.GetPosition(
                        adv_target_hit.GetDetectorID(), self.a, self.b
                    )
                    hit_collection["pos"][0].append(self.a.X())
                    hit_collection["pos"][1].append(self.a.Y())
                    hit_collection["pos"][2].append(self.a.Z())

                    hit_collection["B"][0].append(self.b.X())
                    hit_collection["B"][1].append(self.b.Y())
                    hit_collection["B"][2].append(self.b.Z())

                    hit_collection["d"][0].append(self.AdvTarget_dx)
                    hit_collection["d"][1].append(self.AdvTarget_dy)
                    hit_collection["d"][2].append(self.AdvTarget_dz)

                    hit_collection["vert"].append(adv_target_hit.isVertical())
                    hit_collection["index"].append(i_hit)

                    hit_collection["detectorID"].append(adv_target_hit.GetDetectorID())

                    if self.mask_plane:
                        if (
                            adv_target_hit.isVertical() == 0
                            and adv_target_hit.GetStation() in mask_plane_ZY
                        ) or (
                            adv_target_hit.isVertical()
                            and adv_target_hit.GetStation() in mask_plane_ZX
                        ):
                            hit_collection["mask"].append(True)
                        else:
                            hit_collection["mask"].append(False)
                    else:
                        hit_collection["mask"].append(False)

                    times = []

                    if self.isMC:
                        times.append(adv_target_hit.GetTime())  # already in ns
                    else:
                        times.append(adv_target_hit.GetTime() * 6.25)  # tdc2ns
                    hit_collection["time"].append(times)

        # If no hits, return
        if len(hit_collection["pos"][0]) == 0:
            return

        # Make the hit collection numpy arrays.
        for key, item in hit_collection.items():
            if key == "vert":
                this_dtype = np.bool_
            elif key == "mask":
                this_dtype = np.bool_
            elif key == "index" or key == "detectorID":
                this_dtype = np.int32
            elif key != "time":
                this_dtype = np.float32
            if key == "time":
                length = max(map(len, item))
                hit_collection[key] = np.stack(
                    np.array([xi + [None] * (length - len(xi)) for xi in item]), axis=1
                )
            else:
                hit_collection[key] = np.array(item, dtype=this_dtype)

        # Reconstruct muons until there are not enough hits in downstream muon filter
        for i_muon in range(self.max_reco_muons):
            n_planes_ZY = numPlanesHit(
                hit_collection["detectorID"][~hit_collection["vert"]]
            )
            if n_planes_ZY < self.min_planes_hit:
                break

            n_planes_ZX = numPlanesHit(
                hit_collection["detectorID"][hit_collection["vert"]]
            )
            if n_planes_ZX < self.min_planes_hit:
                break

            # Get hits in hough transform format
            adv_target_hits_horizontal = np.logical_and(
                ~hit_collection["vert"], ~hit_collection["mask"]
            )
            adv_target_hits_vertical = np.logical_and(
                hit_collection["vert"], ~hit_collection["mask"]
            )

            ZY = np.dstack(
                [
                    hit_collection["pos"][2][adv_target_hits_horizontal],
                    hit_collection["pos"][1][adv_target_hits_horizontal],
                ]
            )[0]

            d_ZY = np.dstack(
                [
                    hit_collection["d"][2][adv_target_hits_horizontal],
                    hit_collection["d"][1][adv_target_hits_horizontal],
                ]
            )[0]

            ZX = np.dstack(
                [
                    hit_collection["pos"][2][adv_target_hits_vertical],
                    hit_collection["pos"][0][adv_target_hits_vertical],
                ]
            )[0]

            d_ZX = np.dstack(
                [
                    hit_collection["d"][2][adv_target_hits_vertical],
                    hit_collection["d"][0][adv_target_hits_vertical],
                ]
            )[0]

            is_scaled = False
            ZY_hough = self.h_ZY.fit_randomize(
                ZY, d_ZY, self.n_random, is_scaled, self.draw
            )
            ZX_hough = self.h_ZX.fit_randomize(
                ZX, d_ZX, self.n_random, is_scaled, self.draw
            )

            tol = self.tolerance
            # Special treatment for events with low hit occupancy - increase tolerance
            # For Scifi-only tracks
            if len(hit_collection["detectorID"]) <= self.max_n_AdvTarget_hits:
                # as there are masked Scifi planes, make sure to use hit counts before the masking
                if (
                    max(N_plane_ZX.values()) <= self.max_n_hits_plane
                    and max(N_plane_ZY.values()) <= self.max_n_hits_plane
                ):
                    tol = 5 * self.tolerance

            # Check if track intersects minimum number of hits in each plane.
            track_hits_for_triplet_ZY = hit_finder(
                ZY_hough[0],
                ZY_hough[1],
                np.dstack(
                    [
                        hit_collection["pos"][2][~hit_collection["vert"]],
                        hit_collection["pos"][1][~hit_collection["vert"]],
                    ]
                ),
                np.dstack(
                    [
                        hit_collection["d"][2][~hit_collection["vert"]],
                        hit_collection["d"][1][~hit_collection["vert"]],
                    ]
                ),
                tol,
            )

            track_hits_for_triplet_ZX = hit_finder(
                ZX_hough[0],
                ZX_hough[1],
                np.dstack(
                    [
                        hit_collection["pos"][2][hit_collection["vert"]],
                        hit_collection["pos"][0][hit_collection["vert"]],
                    ]
                ),
                np.dstack(
                    [
                        hit_collection["d"][2][hit_collection["vert"]],
                        hit_collection["d"][0][hit_collection["vert"]],
                    ]
                ),
                tol,
            )

            n_planes_hit_ZY = numPlanesHit(
                hit_collection["detectorID"][~hit_collection["vert"]][
                    track_hits_for_triplet_ZY
                ]
            )

            n_planes_hit_ZX = numPlanesHit(
                hit_collection["detectorID"][hit_collection["vert"]][
                    track_hits_for_triplet_ZX
                ]
            )

            if (
                n_planes_hit_ZY < self.min_planes_hit
                or n_planes_hit_ZX < self.min_planes_hit
            ):
                break

            #                print("Found {0} downstream ZX planes associated to muon track".format(n_planes_ds_ZX))
            #                print("Found {0} downstream ZY planes associated to muon track".format(n_planes_ds_ZY))

            # This time with all the hits, not just triplet condition.
            track_hits_ZY = hit_finder(
                ZY_hough[0],
                ZY_hough[1],
                np.dstack(
                    [
                        hit_collection["pos"][2][~hit_collection["vert"]],
                        hit_collection["pos"][1][~hit_collection["vert"]],
                    ]
                ),
                np.dstack(
                    [
                        hit_collection["d"][2][~hit_collection["vert"]],
                        hit_collection["d"][1][~hit_collection["vert"]],
                    ]
                ),
                tol,
            )

            track_hits_ZX = hit_finder(
                ZX_hough[0],
                ZX_hough[1],
                np.dstack(
                    [
                        hit_collection["pos"][2][hit_collection["vert"]],
                        hit_collection["pos"][0][hit_collection["vert"]],
                    ]
                ),
                np.dstack(
                    [
                        hit_collection["d"][2][hit_collection["vert"]],
                        hit_collection["d"][0][hit_collection["vert"]],
                    ]
                ),
                tol,
            )

            # Onto Kalman fitter (based on SndlhcTracking.py)
            posM = ROOT.TVector3(0, 0, 0.0)
            momM = ROOT.TVector3(0, 0, 100.0)  # default track with high momentum

            # approximate covariance
            covM = ROOT.TMatrixDSym(6)
            res = self.kalman_sigmaAdvTarget_spatial
            for i in range(3):
                covM[i][i] = res * res
            for i in range(3, 6):
                covM[i][i] = ROOT.TMath.Power(res / (4.0 * 2.0) / ROOT.TMath.Sqrt(3), 2)
            rep = ROOT.genfit.RKTrackRep(13)

            # start state
            state = ROOT.genfit.MeasuredStateOnPlane(rep)
            rep.setPosMomCov(state, posM, momM, covM)

            # create track
            seedState = ROOT.TVectorD(6)
            seedCov = ROOT.TMatrixDSym(6)
            rep.get6DStateCov(state, seedState, seedCov)

            theTrack = ROOT.genfit.Track(rep, seedState, seedCov)

            # Sort measurements in Z
            hit_z = np.concatenate(
                [
                    hit_collection["pos"][2][hit_collection["vert"]][track_hits_ZX],
                    hit_collection["pos"][2][~hit_collection["vert"]][track_hits_ZY],
                ]
            )

            hit_A0 = np.concatenate(
                [
                    hit_collection["pos"][0][hit_collection["vert"]][track_hits_ZX],
                    hit_collection["pos"][0][~hit_collection["vert"]][track_hits_ZY],
                ]
            )

            hit_A1 = np.concatenate(
                [
                    hit_collection["pos"][1][hit_collection["vert"]][track_hits_ZX],
                    hit_collection["pos"][1][~hit_collection["vert"]][track_hits_ZY],
                ]
            )

            hit_B0 = np.concatenate(
                [
                    hit_collection["B"][0][hit_collection["vert"]][track_hits_ZX],
                    hit_collection["B"][0][~hit_collection["vert"]][track_hits_ZY],
                ]
            )

            hit_B1 = np.concatenate(
                [
                    hit_collection["B"][1][hit_collection["vert"]][track_hits_ZX],
                    hit_collection["B"][1][~hit_collection["vert"]][track_hits_ZY],
                ]
            )

            hit_B2 = np.concatenate(
                [
                    hit_collection["B"][2][hit_collection["vert"]][track_hits_ZX],
                    hit_collection["B"][2][~hit_collection["vert"]][track_hits_ZY],
                ]
            )

            hit_detid = np.concatenate(
                [
                    hit_collection["detectorID"][hit_collection["vert"]][track_hits_ZX],
                    hit_collection["detectorID"][~hit_collection["vert"]][
                        track_hits_ZY
                    ],
                ]
            )

            kalman_spatial_sigma = np.concatenate(
                [
                    hit_collection["d"][0][hit_collection["vert"]][track_hits_ZX]
                    / 12**0.5,
                    hit_collection["d"][1][~hit_collection["vert"]][track_hits_ZY]
                    / 12**0.5,
                ]
            )

            # Maximum distance. Use (d_xy/2**2 + d_z/2**2)**0.5
            kalman_max_dis = np.concatenate(
                [
                    (
                        (
                            hit_collection["d"][0][hit_collection["vert"]][
                                track_hits_ZX
                            ]
                            / 2.0
                        )
                        ** 2
                        + (
                            hit_collection["d"][2][hit_collection["vert"]][
                                track_hits_ZX
                            ]
                            / 2.0
                        )
                        ** 2
                    )
                    ** 0.5,
                    (
                        (
                            hit_collection["d"][1][~hit_collection["vert"]][
                                track_hits_ZY
                            ]
                            / 2.0
                        )
                        ** 2
                        + (
                            hit_collection["d"][2][~hit_collection["vert"]][
                                track_hits_ZY
                            ]
                            / 2.0
                        )
                        ** 2
                    )
                    ** 0.5,
                ]
            )

            hitID = 0  # Does it matter? We don't have a global hit ID.

            hit_time = {}
            for ch in range(hit_collection["time"].shape[0]):
                hit_time[ch] = np.concatenate(
                    [
                        hit_collection["time"][ch][hit_collection["vert"]][
                            track_hits_ZX
                        ],
                        hit_collection["time"][ch][~hit_collection["vert"]][
                            track_hits_ZY
                        ],
                    ]
                )

            for i_z_sorted in hit_z.argsort():
                tp = ROOT.genfit.TrackPoint()
                hitCov = ROOT.TMatrixDSym(7)
                hitCov[6][6] = kalman_spatial_sigma[i_z_sorted] ** 2

                measurement = ROOT.genfit.WireMeasurement(
                    ROOT.TVectorD(
                        7,
                        array(
                            "d",
                            [
                                hit_A0[i_z_sorted],
                                hit_A1[i_z_sorted],
                                hit_z[i_z_sorted],
                                hit_B0[i_z_sorted],
                                hit_B1[i_z_sorted],
                                hit_B2[i_z_sorted],
                                0.0,
                            ],
                        ),
                    ),
                    hitCov,
                    1,  # detid?
                    6,  # hitid?
                    tp,
                )

                measurement.setMaxDistance(kalman_max_dis[i_z_sorted])
                measurement.setDetId(int(hit_detid[i_z_sorted]))
                measurement.setHitId(int(hitID))
                hitID += 1
                tp.addRawMeasurement(measurement)
                theTrack.insertPoint(tp)

            if not theTrack.checkConsistency():
                theTrack.Delete()
                raise RuntimeError("Kalman fitter track consistency check failed.")

            # do the fit
            self.kalman_fitter.processTrack(
                theTrack
            )  # processTrackWithRep(theTrack,rep,True)

            fitStatus = theTrack.getFitStatus()
            if not fitStatus.isFitConverged() and 0 > 1:
                theTrack.Delete()
                raise RuntimeError("Kalman fit did not converge.")

            # Now save the track if fit converged!
            theTrack.SetUniqueID(self.track_type)
            if fitStatus.isFitConverged():
                if self.genfitTrack:
                    self.kalman_tracks.Add(theTrack)
                else:
                    # Load items into snd track class object
                    this_track = ROOT.sndRecoTrack(theTrack)
                    pointTimes = ROOT.std.vector(ROOT.std.vector("float"))()
                    for n, i_z_sorted in enumerate(hit_z.argsort()):
                        t_per_hit = []
                        for ch in range(len(hit_time)):
                            if hit_time[ch][i_z_sorted] is not None:
                                t_per_hit.append(hit_time[ch][i_z_sorted])
                        pointTimes.push_back(t_per_hit)
                    this_track.setRawMeasTimes(pointTimes)
                    this_track.setTrackType(self.track_type)
                    # Save the track in sndRecoTrack format
                    self.kalman_tracks[i_muon] = this_track
                    # Delete the Kalman track object
                    theTrack.Delete()

            # Remove track hits and try to find an additional track
            # Find array index to be removed
            index_to_remove_ZX = np.where(
                np.in1d(
                    hit_collection["detectorID"],
                    hit_collection["detectorID"][hit_collection["vert"]][track_hits_ZX],
                )
            )[0]
            index_to_remove_ZY = np.where(
                np.in1d(
                    hit_collection["detectorID"],
                    hit_collection["detectorID"][~hit_collection["vert"]][
                        track_hits_ZY
                    ],
                )
            )[0]

            index_to_remove = np.concatenate([index_to_remove_ZX, index_to_remove_ZY])

            # Remove dictionary entries
            for key in hit_collection.keys():
                if len(hit_collection[key].shape) == 1:
                    hit_collection[key] = np.delete(
                        hit_collection[key], index_to_remove
                    )
                elif len(hit_collection[key].shape) == 2:
                    hit_collection[key] = np.delete(
                        hit_collection[key], index_to_remove, axis=1
                    )
                else:
                    raise Exception(
                        "Wrong number of dimensions found when deleting hits in iterative muon identification algorithm."
                    )

    def FinishTask(self):
        """End of the task."""
        print("Processed", self.events_run)
        if not self.genfitTrack:
            self.kalman_tracks.Delete()
        else:
            pass

    # this is a copy on SndlhcTracking function with small adjustments in event object names to make it work here
    # FIXME Should find a way to use this function straight from the SndlhcTracking!
    def scifiCluster(self):
        """Make clusters."""
        clusters = []
        hitDict = {}
        for k in range(self.AdvTargetHits.GetEntries()):
            d = self.AdvTargetHits[k]
            if not d.isValid():
                continue
            hitDict[d.GetDetectorID()] = k
        hitList = list(hitDict.keys())
        if len(hitList) > 0:
            hitList.sort()
            tmp = [hitList[0]]
            cprev = hitList[0]
            ncl = 0
            last = len(hitList) - 1
            hitvector = ROOT.std.vector("sndadv_target_hit*")()
            for i in range(len(hitList)):
                if i == 0 and len(hitList) > 1:
                    continue
                c = hitList[i]
                neighbour = False
                if (c - cprev) == 1:  # does not account for neighbours across sipms
                    neighbour = True
                    tmp.append(c)
                if not neighbour or c == hitList[last]:
                    first = tmp[0]
                    N = len(tmp)
                    hitvector.clear()
                    for aHit in tmp:
                        hitvector.push_back(self.AdvTargetHits[hitDict[aHit]])
                    aCluster = ROOT.sndCluster(
                        first, N, hitvector, self.AdvTargetDet, False
                    )
                    clusters.append(aCluster)
                    if c != hitList[last]:
                        ncl += 1
                        tmp = [c]
                    elif not neighbour:  # save last channel
                        hitvector.clear()
                        hitvector.push_back(self.AdvTargetHits[hitDict[c]])
                        aCluster = ROOT.sndCluster(
                            c, 1, hitvector, self.AdvTargetDet, False
                        )
                        clusters.append(aCluster)
                cprev = c
        self.clusScifi.Delete()
        for c in clusters:
            self.clusScifi.Add(c)
