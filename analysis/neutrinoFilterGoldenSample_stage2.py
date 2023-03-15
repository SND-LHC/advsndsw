import ROOT
import numpy as np

import SndlhcGeo
import shipunit as u

import argparse
parser = argparse.ArgumentParser()

parser.add_argument("-f", dest = "inputFile", required = True)
parser.add_argument("-o", dest = "outputFile", required = True)
parser.add_argument("-t", dest = "trackFile", required = True)
parser.add_argument("-g", dest = "geoFile", required = True)

args = parser.parse_args()

snd_geo = SndlhcGeo.GeoInterface(args.geoFile)
scifiDet = ROOT.gROOT.GetListOfGlobals().FindObject('Scifi')
muFilterDet = ROOT.gROOT.GetListOfGlobals().FindObject('MuFilter')

# Hard-coded?!
TDC2ns = 1E9/160.316E6
dsV = muFilterDet.GetConfParF("MuFilter/DsPropSpeed")
dsL = muFilterDet.GetConfParF("MuFilter/UpstreamBarX")

# Set up TTrees
isMC = False
treeName = "rawConv"

ch = ROOT.TChain(treeName)
ch.Add(args.inputFile)

if ch.GetEntries() == 0 :
    treeName = "cbmsim"
    isMC = True
    TDC2ns = 1. # In the simulation hit times are always in ns
    del ch
    ch = ROOT.TChain(treeName)
    ch.Add(args.inputFile)

if ch.GetEntries() == 0 :
    print("Chain is empty. Exitting")
    exit(-1)

ch_tracks = ROOT.TChain(treeName)
ch_tracks.Add(args.trackFile)
ch.AddFriend(ch_tracks)

dts = []

# Set up cuts
cuts = []

################################################################################
# Event has one reconstructed DS track
################################################################################
def eventHasOneTrack(event) :
    if event.Reco_MuonTracks.GetEntries() >= 1 :
        return True
    else :
        return False
cuts.append(["Event has one reconstructed DS track", eventHasOneTrack])
################################################################################
# Event direction cut
################################################################################
def direction(event) :

    t_scifi = []
    t_muon = []

    for hit in event.Digi_ScifiHits :
        if not hit.isValid() : 
            continue
        t_scifi.append(hit.GetTime())

    for hit in event.Digi_MuFilterHits :
        if hit.GetSystem() != 3 :
            continue
        if not hit.isValid() :
            continue
        if len(hit.GetAllTimes()) == 0 :
            continue

        mufi_hit_times = []
        
        for hit_times in hit.GetAllTimes() :
            mufi_hit_times.append(hit_times.second)
        t_muon.append(np.average(mufi_hit_times))

    if len(t_scifi) == 0 :
        return 999
    if len(t_muon) == 0 :
        return 999

    return (np.sort(t_scifi)[0] - np.sort(t_muon)[-1]) < 0

#delta_t_cut = -9
#def direction(event) :
#    if not isMC :
#        scifiDet.InitEvent(event.EventHeader)
#        muFilterDet.InitEvent(event.EventHeader)
#    
#    a = ROOT.TVector3()
#    b = ROOT.TVector3()
#
#    t_scifi = []
#    for hit in event.Digi_ScifiHits :
#        if not hit.isValid() : 
#            continue
#        scifiDet.GetSiPMPosition(hit.GetDetectorID(), a, b)
#        z = (a.Z() + b.Z())/2.
#        t = hit.GetTime()*TDC2ns
#        if not isMC :
#            t = scifiDet.GetCorrectedTime(hit.GetDetectorID(), t, 0)
#        t_scifi.append(t - z/u.speedOfLight)
#
#    t_muFilter = []
#    for hit in event.Digi_MuFilterHits :
#        if hit.GetSystem() != 3 :
#            continue
#        if hit.GetPlane() != 2 :
#            continue
#        if hit.isVertical() : 
#            continue
#        if not hit.isValid() :
#            continue
#        
#        muFilterDet.GetPosition(hit.GetDetectorID(), a, b)
#        z = (a.Z() + b.Z())/2.
#
#        t = [hit.GetTime(i)*TDC2ns for i in range(2)]
#        if not isMC :
#            for i in range(2) :
#                t[i] = muFilterDet.GetCorrectedtime(hit.GetDetectorID(), i, t[i], 0)
#        t = np.mean(t)
#        t_muFilter.append(t - dsL/dsV*0.5 -z/u.speedOfLight)
#        
#    t_min_scifi = np.min(t_scifi)
#    t_min_muFilter = np.min(t_muFilter)
#
#    delta = t_min_muFilter - t_min_scifi
#
#    dts.append(delta)
#
#    if delta < delta_t_cut :
#        return False
#    return True
cuts.append(["Event direction", direction])

################################################################################
# Track intercepts first SciFi plane within 5 cm of the edge
################################################################################
d_scifi_fiducial = 5
def trackInScifiFiducial(event) :
    track = event.Reco_MuonTracks.At(0)

    track_start = track.getStart()
    track_stop = track.getStop()

    slope_X = (track_stop.X() - track_start.X())/(track_stop.Z() - track_start.Z())
    slope_Y = (track_stop.Y() - track_start.Y())/(track_stop.Z() - track_start.Z())

    z_first_scifi = 300.

    track_first_scifi_extrap = ROOT.TVector3(track_start.X() + slope_X*(z_first_scifi - track_start.Z()), track_start.Y() + slope_Y*(z_first_scifi - track_start.Z()), z_first_scifi)

    if track_first_scifi_extrap.Y() < 15.5 + d_scifi_fiducial :
        return False
    if track_first_scifi_extrap.Y() > 15.5+39 - d_scifi_fiducial :
        return False
    if track_first_scifi_extrap.X() > -8 - d_scifi_fiducial :
        return False
    if track_first_scifi_extrap.X() < -8 -39 + d_scifi_fiducial:
        return False
    return True
cuts.append(["Track intercepts first SciFi plane < {0} cm from edge".format(d_scifi_fiducial), trackInScifiFiducial])

################################################################################
# SciFi hit to DS track DOCA per plane < 3 cm for both projections
################################################################################
sum_min_dca_cut = 3
def sum_min_dca(event) :
    
    track = event.Reco_MuonTracks.At(0)

    track_start = track.getStart()
    track_stop = track.getStop()

    slope_X = (track_stop.X() - track_start.X())/(track_stop.Z() - track_start.Z())
    slope_Y = (track_stop.Y() - track_start.Y())/(track_stop.Z() - track_start.Z())

    a = ROOT.TVector3()
    b = ROOT.TVector3()

    min_dca_ver = [1e10]*5
    min_dca_hor = [1e10]*5
    
    for hit in event.Digi_ScifiHits :
        if not hit.isValid() :
            continue
            
        scifiDet.GetSiPMPosition(hit.GetDetectorID(), a, b)

        if hit.isVertical() :
            slope_XY = slope_X
            track_start_XY = track_start.X()
            track_stop_XY = track_stop.X()
            hit_pos = [(a.Z() + b.Z())/2., (a.X() + b.X())/2.]
        else :
            slope_XY = slope_Y
            track_start_XY = track_start.Y()
            track_stop_XY = track_stop.Y()
            hit_pos = [(a.Z() + b.Z())/2., (a.Y() + b.Y())/2.]

        if hit.isVertical() :
            this_d = np.abs(track_start_XY + slope_XY*(hit_pos[0]-track_start.Z()) - hit_pos[1])
            if this_d < min_dca_ver[hit.GetStation()-1] :
                min_dca_ver[hit.GetStation()-1] = this_d
        else :
            this_d = np.abs(track_start_XY + slope_XY*(hit_pos[0]-track_start.Z()) - hit_pos[1])
            if this_d < min_dca_hor[hit.GetStation()-1] :
                min_dca_hor[hit.GetStation()-1] = this_d

    min_dca_ver = np.array(min_dca_ver)
    min_dca_hor = np.array(min_dca_hor)
    
    n_ver = 5 - np.sum(min_dca_ver > 9.9e9)
    n_hor = 5 - np.sum(min_dca_hor > 9.9e9)

    min_dca_ver[min_dca_ver > 9.9e9] = 0
    min_dca_hor[min_dca_hor > 9.9e9] = 0

    avg_ver = np.sum(min_dca_ver)/n_ver
    avg_hor = np.sum(min_dca_hor)/n_hor

    return np.max([avg_ver, avg_hor]) <= sum_min_dca_cut
    
    min_dca_ver[min_dca_ver > 9.9e9] = 0
    min_dca_hor[min_dca_hor > 9.9e9] = 0

    return np.sum([min_dca_ver, min_dca_hor]) <= sum_min_dca_cut
cuts.append(["Sum of min DOCA per station < {0} cm".format(sum_min_dca_cut), sum_min_dca])

################################################################################
# At least 35 SciFi hits
################################################################################
min_scifi_hits_cut = 35
def min_scifi_hits(event) :
    n_hits = 0
    for hit in event.Digi_ScifiHits :
        if hit.isValid() :
            n_hits += 1
            if n_hits > 35 :
                return True
    return False
cuts.append(["More than {0} SciFi hits".format(min_scifi_hits_cut), min_scifi_hits])

################################################################################
# Min QDC
################################################################################
min_QDC_data = 600
min_QDC_MC = 700
def min_US_QDC(event) :
    US_QDC = 0
    for hit in event.Digi_MuFilterHits :
        if hit.GetSystem() != 2 :
            continue
        if not hit.isValid() :
            continue
        for key, value in hit.GetAllSignals() :
            US_QDC += value
            if isMC and (US_QDC > min_QDC_MC) :
                return True
            if (not isMC) and (US_QDC > min_QDC_data) :
                return True
    return False
cuts.append(["US QDC larger than {0} ({1}) for data (MC)".format(min_QDC_data, min_QDC_MC), min_US_QDC])

################################################################################
# Max DS activity < 10
################################################################################
max_DS_hits_cut = 10
def max_DS_hits(event) :
    DS_hits = 0.
    for hit in event.Digi_MuFilterHits :
        if hit.GetSystem() != 3 :
            continue
        if not hit.isValid() :
            continue
        
        if hit.GetPlane() == 3 :
            DS_hits += 1.
        else :
            DS_hits += 0.5
    if DS_hits > max_DS_hits_cut :
        return False
    return True
cuts.append(["Number of DS hits per projection is < {0}".format(max_DS_hits_cut), max_DS_hits])

################################################################################
# END CUT DEFINITONS
################################################################################


# Set up cut flow histogram
ch.GetEntry(0)
f = ch.GetFile()

cut_flow = f.Get("cutFlow")
cut_flow_extended = ROOT.TH1D(cut_flow.GetName()+"_extended", cut_flow.GetTitle(), cut_flow.GetNbinsX()+len(cuts), 0, cut_flow.GetNbinsX()+len(cuts))

for i in range(1, cut_flow.GetNbinsX()+1) :
    cut_flow_extended.SetBinContent(i, cut_flow.GetBinContent(i))
    cut_flow_extended.GetXaxis().SetBinLabel(i, cut_flow.GetXaxis().GetBinLabel(i))

for i in range(len(cuts)) :
    cut_flow_extended.GetXaxis().SetBinLabel(i+1+cut_flow.GetNbinsX(), cuts[i][0])

# Set up output file
output_file = ROOT.TFile(args.outputFile, "RECREATE")
output_tree = ch.CloneTree(0)

# Copy branch list
branch_list = f.Get("BranchList")
branch_list_copy = branch_list.Clone()
branch_list_copy.Write("BranchList", 1)

import matplotlib.pyplot as plt

#dca_for_hist = []
i_pass = 0
for event in ch :
    passes_cut = True
    for i_cut, cut in enumerate(cuts) :
        if cut[1](event) :
            cut_flow_extended.Fill(cut_flow.GetNbinsX()+i_cut)
        else :
            passes_cut = False
            break
    if passes_cut :
        output_tree.Fill()
    else :
        continue

    print("EVENT {0}".format(i_pass))
    i_pass +=1
 #   dca_for_hist.append(sum_min_dca(event))
 #   print(dca_for_hist[-1])

#plt.figure()
#plt.hist(dca_for_hist, bins = 100)
#plt.show()

plt.figure()
plt.hist(dts, bins = 100)
plt.show()

cut_flow_extended.Write()
output_file.Write()
output_file.Close()
