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

#if ch.GetEntries() == 0 :
#    print("Chain is empty. Exitting")
#    exit(-1)

ch_tracks = ROOT.TChain(treeName)
ch_tracks.Add(args.trackFile)

n_events = min(ch.GetEntries(), ch_tracks.GetEntries())
ch.AddFriend(ch_tracks)

# Set up cuts
cuts = []

################################################################################
# No hits in second SciFi plane
################################################################################
#def no_hits_second_sf(event) :
#    n_hits_2_sf = 0
#    for hit in event.Digi_ScifiHits :
#        if hit.GetStation() == 2 :
#            if hit.isValid() :
#                n_hits_2_sf += 1
#    return n_hits_2_sf == 0, n_hits_2_sf
#cuts.append(["No hits in second Scifi plane", no_hits_second_sf, "no_sf2", 100, 0, 100])

################################################################################
# Event has one reconstructed DS track
################################################################################
def eventHasOneTrack(event) :
    ret = False
    n_tracks = event.Reco_MuonTracks.GetEntries()
    if n_tracks >= 1 :
        ret =  True
    return ret, n_tracks 
cuts.append(["Event has one reconstructed DS track", eventHasOneTrack, "n_DS_tracks", 3, 0, 3])
################################################################################
# Event direction cut
################################################################################
def direction(event) :
    if not isMC :
        scifiDet.InitEvent(event.EventHeader)
        muFilterDet.InitEvent(event.EventHeader)

    t_scifi = []
    t_muon = []

    for hit in event.Digi_ScifiHits :
        if not hit.isValid() : 
            continue

        t = hit.GetTime()*TDC2ns
#        if not isMC :
#            t = scifiDet.GetCorrectedTime(hit.GetDetectorID(), t, 0)
        
        t_scifi.append(t)

    for hit in event.Digi_MuFilterHits :
        if hit.GetSystem() != 3 :
            continue
        if not hit.isValid() :
            continue

        if hit.isVertical() :
            t = hit.GetTime(0)*TDC2ns
#            if not isMC :
#                t = muFilterDet.GetCorrectedTime(hit.GetDetectorID(), 0, t, 0)
        else :
            t = [hit.GetTime(i)*TDC2ns for i in range(2)]

 #           if not isMC :
 #               for i in range(2) :
 #                   t[i] = muFilterDet.GetCorrectedTime(hit.GetDetectorID(), i, t[i], 0)
            t = np.mean(t)
            t_muon.append(t)

    delta_t = np.sort(t_scifi)[0] - np.sort(t_muon)[-1]
    ret = delta_t < 0

    return ret, delta_t

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
#    t_muFilter = []
#
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
#                t[i] = muFilterDet.GetCorrectedTime(hit.GetDetectorID(), i, t[i], 0)
#        t = np.array(t)
#        t = t[np.abs(t) < 500]
#        t = np.mean(t)
#        t_muFilter.append(t - dsL/dsV*0.5 -z/u.speedOfLight)
#        
#    t_min_scifi = np.mean(t_scifi)
#    t_min_muFilter = np.mean(t_muFilter)
#
#    delta = t_min_muFilter - t_min_scifi
#
#    dts.append(delta)
#
#    if delta < delta_t_cut :
#        return False
#    return True
cuts.append(["Event direction", direction, "event_dir", 60, -30, 30])

################################################################################
# Track intercepts first SciFi plane within 5 cm of the edge
################################################################################
d_scifi_fiducial = 5
def trackInScifiFiducial(event) :
    if len(event.Reco_MuonTracks) < 1 :
        return False, -999
    track = event.Reco_MuonTracks.At(0)

    track_start = track.getStart()
    track_stop = track.getStop()

    slope_X = (track_stop.X() - track_start.X())/(track_stop.Z() - track_start.Z())
    slope_Y = (track_stop.Y() - track_start.Y())/(track_stop.Z() - track_start.Z())

    z_first_scifi = 300.

    track_first_scifi_extrap = ROOT.TVector3(track_start.X() + slope_X*(z_first_scifi - track_start.Z()), track_start.Y() + slope_Y*(z_first_scifi - track_start.Z()), z_first_scifi)

    min_d = 1e6
    if -15.5 + track_first_scifi_extrap.Y() < min_d :
        min_d = -15.5 + track_first_scifi_extrap.Y() < min_d
    if + 15.5+39 - track_first_scifi_extrap.Y()  < min_d :
        min_d = + 15.5+39 - track_first_scifi_extrap.Y()
    if -8 - track_first_scifi_extrap.X() < min_d:
        min_d = -8  - track_first_scifi_extrap.X()
    if +8 +39 + track_first_scifi_extrap.X() < min_d :
        min_d = +8 +39 + track_first_scifi_extrap.X()
    ret = True
    if track_first_scifi_extrap.Y() < 15.5 + d_scifi_fiducial :
        ret =  False
    if track_first_scifi_extrap.Y() > 15.5+39 - d_scifi_fiducial :
        ret = False
    if track_first_scifi_extrap.X() > -8 - d_scifi_fiducial :
        ret = False
    if track_first_scifi_extrap.X() < -8 -39 + d_scifi_fiducial:
        ret = False
    return ret, min_d
cuts.append(["Track intercepts first SciFi plane < {0} cm from edge".format(d_scifi_fiducial), trackInScifiFiducial, "ds_track_extrap_fiducial", 100, -50, 50])

################################################################################
# SciFi hit to DS track DOCA per plane < 3 cm for both projections
################################################################################
sum_min_dca_cut = 3
def sum_min_dca(event) :
    if len(event.Reco_MuonTracks) < 1 :
        return False, 999
    
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

    ret_value = np.max([avg_ver, avg_hor])
    ret = ret_value <= sum_min_dca_cut

    return ret, ret_value
cuts.append(["Sum of min DOCA per station < {0} cm".format(sum_min_dca_cut), sum_min_dca, "doca", 30, 0, 30])

################################################################################
# At least 35 SciFi hits
################################################################################
min_scifi_hits_cut = 35
def min_scifi_hits(event) :
    n_hits = 0
    ret = False
    for hit in event.Digi_ScifiHits :
        if hit.isValid() :
            n_hits += 1
            if n_hits > 35 :
                ret = True
    return ret, n_hits
cuts.append(["More than {0} SciFi hits".format(min_scifi_hits_cut), min_scifi_hits, "scifi_nhits", 100, 0, 3000])

################################################################################
# Min QDC
################################################################################
min_QDC_data = 600
min_QDC_MC = 700
def min_US_QDC(event) :
    US_QDC = 0
    ret = False
    for hit in event.Digi_MuFilterHits :
        if hit.GetSystem() != 2 :
            continue
        if not hit.isValid() :
            continue
        for key, value in hit.GetAllSignals() :
#            if (key >= 8) and (hit.GetPlane() == 1) : 
#                continue
            US_QDC += value
            if isMC and (US_QDC > min_QDC_MC) :
                ret = True
            if (not isMC) and (US_QDC > min_QDC_data) :
                ret = True
    if isMC :
        return ret, US_QDC - 100
    else :
        return ret, US_QDC
cuts.append(["US QDC larger than {0} ({1}) for data (MC)".format(min_QDC_data, min_QDC_MC), min_US_QDC, "US_QDC", 100, 0, 4000])

################################################################################
# Max DS activity < 10
################################################################################
max_DS_hits_cut = 10
def max_DS_hits(event) :
    DS_hits = 0.
    ret = True
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
        ret =  False
    return ret, DS_hits
cuts.append(["Number of DS hits per projection is < {0}".format(max_DS_hits_cut), max_DS_hits, "DS_hits", 100, 0, 100])

################################################################################
# END CUT DEFINITONS
################################################################################


ch.GetEntry(0)
f = ch.GetFile()

# Set up output file
output_file = ROOT.TFile(args.outputFile, "RECREATE")
output_tree = ch.CloneTree(0)

# Copy branch list
branch_list = f.Get("BranchList")
branch_list_copy = branch_list.Clone()
branch_list_copy.Write("BranchList", 1)

# Set up cut flow histogram
cut_flow = f.Get("cutFlow")
cut_flow_extended = ROOT.TH1D(cut_flow.GetName()+"_extended", cut_flow.GetTitle(), cut_flow.GetNbinsX()+len(cuts), 0, cut_flow.GetNbinsX()+len(cuts))

for i in range(1, cut_flow.GetNbinsX()+1) :
    cut_flow_extended.SetBinContent(i, cut_flow.GetBinContent(i))
    cut_flow_extended.GetXaxis().SetBinLabel(i, cut_flow.GetXaxis().GetBinLabel(i))

for i in range(len(cuts)) :
    cut_flow_extended.GetXaxis().SetBinLabel(i+1+cut_flow.GetNbinsX(), cuts[i][0])

# Cut-by-cut histograms
cut_by_cut_var_histos = []
for i_cut in range(-1, len(cuts)) :
    this_cut_by_cut_var_histos = []
    for this_cut_name, cut_function, short_name, nbins, range_start, range_end in cuts :
        print("Initializing", short_name, nbins, range_start, range_end)
        this_cut_by_cut_var_histos.append(ROOT.TH1D(str(i_cut)+"_"+short_name+"_0",
                                                    short_name,
                                                    nbins, range_start, range_end))
    cut_by_cut_var_histos.append(this_cut_by_cut_var_histos)

if isMC :
    suffix_dict = {}
    suffix_dict[0] = "nueCC"
    suffix_dict[1] = "numuCC"
    suffix_dict[2] = "nutauCC"
    suffix_dict[3] = "NC"
    suffix_dict[4] = "Other"
    
    cut_by_cut_truth_histos = []
    for i_species in range(5) :
        this_species_histos = []
        for i_cut in range(-1, len(cuts)) :
            this_cut_by_cut_truth_histos = []
            this_cut_by_cut_truth_histos.append(ROOT.TH1D(suffix_dict[i_species]+"_"+str(i_cut)+"_Enu", "Enu", 300, 0, 3000))
            this_cut_by_cut_truth_histos.append(ROOT.TH1D(suffix_dict[i_species]+"_"+str(i_cut)+"_EEM", "ELep", 300, 0, 3000))
            this_cut_by_cut_truth_histos.append(ROOT.TH1D(suffix_dict[i_species]+"_"+str(i_cut)+"_EHad", "EHad", 300, 0, 3000))
            this_cut_by_cut_truth_histos.append(ROOT.TH1D(suffix_dict[i_species]+"_"+str(i_cut)+"_vtxX", "vtxX", 200, -100, 0))
            this_cut_by_cut_truth_histos.append(ROOT.TH1D(suffix_dict[i_species]+"_"+str(i_cut)+"_vtxY", "vtxY", 200, 0, 100))
            this_cut_by_cut_truth_histos.append(ROOT.TH1D(suffix_dict[i_species]+"_"+str(i_cut)+"_vtxZ", "vtxZ", 200, 280, 380))
            
            this_species_histos.append(this_cut_by_cut_truth_histos)
        cut_by_cut_truth_histos.append(this_species_histos)

# N-1
n_minus_1_var_histos = []
for this_cut_name, cut_function, short_name, nbins, range_start, range_end in cuts :
    n_minus_1_var_histos.append(ROOT.TH1D("n_minus_1_"+short_name+"_0", short_name,
                                          nbins, range_start, range_end))


#import matplotlib.pyplot as plt

#dca_for_hist = []
i_pass = 0
passes_cut = [False]*len(cuts)
cut_var = [0.]*len(cuts)
# Set up output summary file for emulsion matching
import csv
header = ["RunNo", "EventNo", "wallNo", "trackSFinterceptX", "trackSFinterceptY", "trackTX", "trackTY", "SFX", "SFY"]
with open("nu_candidates_summary.txt", "w") as output_nu_cand_summary :
    csv_writer = csv.writer(output_nu_cand_summary)
    csv_writer.writerow(header)
    for i_event, event in enumerate(ch) :
        if i_event >= n_events :
            break
        n_cuts_passed = 0
        accept_event = True
        for i_cut, cut in enumerate(cuts) :
            this_cut_passed, this_cut_var = cut[1](event)
            passes_cut[i_cut] = this_cut_passed
            if this_cut_passed :
                n_cuts_passed += 1
            cut_var[i_cut] = this_cut_var
            if accept_event and this_cut_passed :
                cut_flow_extended.Fill(cut_flow.GetNbinsX()+i_cut)
            else :
                accept_event = False
        if accept_event :
            output_tree.Fill()
    
        # Fill histograms
        # Sequential
        for seq_cut in range(-1, len(passes_cut)) :
            if seq_cut >= 0 :
                if not passes_cut[seq_cut] :
                    break
    #        for i_hist in range(len(cut_by_cut_var_histos[seq_cut+1])) :
            for i_cut_var, this_cut_var in enumerate(cut_var) :
                cut_by_cut_var_histos[seq_cut+1][i_cut_var].Fill(this_cut_var)
            if isMC :
                this_species = -1
                if len(event.MCTrack) < 2 :
                    this_species = 4
                else :
                    pdgIn = abs(event.MCTrack[0].GetPdgCode())
                    pdgOut = abs(event.MCTrack[1].GetPdgCode())
                if pdgIn == (pdgOut+1) :
                    if pdgIn == 12 :
                        this_species = 0
                    if pdgIn == 14 :
                        this_species = 1
                    if pdgIn == 16 :
                        this_species = 2
                elif pdgIn == pdgOut :
                    this_species = 3
                else :
                    this_species = 4
                
                cut_by_cut_truth_histos[this_species][seq_cut+1][0].Fill(event.MCTrack[0].GetEnergy())
                if this_species < 4 :
                    cut_by_cut_truth_histos[this_species][seq_cut+1][1].Fill(event.MCTrack[1].GetEnergy())
                    cut_by_cut_truth_histos[this_species][seq_cut+1][2].Fill(event.MCTrack[0].GetEnergy()-event.MCTrack[1].GetEnergy())
                cut_by_cut_truth_histos[this_species][seq_cut+1][3].Fill(event.MCTrack[0].GetStartX())
                cut_by_cut_truth_histos[this_species][seq_cut+1][4].Fill(event.MCTrack[0].GetStartY())
                cut_by_cut_truth_histos[this_species][seq_cut+1][5].Fill(event.MCTrack[0].GetStartZ())
                
        # N-1
        for i_cut in range(len(passes_cut)) :
            if ((not passes_cut[i_cut] and n_cuts_passed == (len(passes_cut)-1))) or (n_cuts_passed == len(passes_cut)) :
                n_minus_1_var_histos[i_cut].Fill(cut_var[i_cut])
    
        if (n_cuts_passed == len(passes_cut)) :
            print("EVENT {0}".format(i_pass))
            SFX = 0.
            SFY = 0.
    
            sfPerStation = np.array([0]*5)
            a = ROOT.TVector3()
            b = ROOT.TVector3()
            
            for hit in event.Digi_ScifiHits :
                if not hit.isValid() :
                    continue
                
                sfPerStation[hit.GetStation()-1] += 1
                scifiDet.GetSiPMPosition(hit.GetDetectorID(), a, b)
                
                SFX += (a.X() + b.X())/2.
                SFY += (a.Y() + b.Y())/2.
    
            SFX /= sfPerStation.sum()
            SFY /= sfPerStation.sum()
    
            wallNo = int(np.argmax(np.cumsum(sfPerStation)/sfPerStation.sum() >= 0.05) + 1)
            print(sfPerStation, wallNo)
    
            wallZ = 0.
            for hit in event.Digi_ScifiHits :
                if not hit.isValid() :
                    continue
                if hit.GetStation() == wallNo :
                    scifiDet.GetSiPMPosition(hit.GetDetectorID(), a, b)
                    wallZ += (a.Z() + b.Z())/2.
            wallZ /= sfPerStation[wallNo-1]

            print(wallZ)
    
            track = event.Reco_MuonTracks.At(0)
    
            track_start = track.getStart()
            track_stop = track.getStop()
    
            slope_X = (track_stop.X() - track_start.X())/(track_stop.Z() - track_start.Z())
            slope_Y = (track_stop.Y() - track_start.Y())/(track_stop.Z() - track_start.Z())
    
            trackSFinterceptX = track_start.X() + slope_X * (wallZ - track_start.Z())
            trackSFinterceptY = track_start.Y() + slope_Y * (wallZ - track_start.Z())
            
    
            summary = []
            try :
                summary.append(event.EventHeader.GetRunId())
                summary.append(event.EventHeader.GetEventNumber())
            except AttributeError :
                summary.append(-1)
                summary.append(-1)
            summary.append(wallNo)
            summary.append(trackSFinterceptX)
            summary.append(trackSFinterceptY)
            summary.append(slope_X)
            summary.append(slope_Y)
            summary.append(SFX)
            summary.append(SFY)
            
            csv_writer.writerow(map( lambda t: t if isinstance(t, int) else "%.3f" % t, summary))

            i_pass +=1
 #   dca_for_hist.append(sum_min_dca(event))
 #   print(dca_for_hist[-1])

#plt.figure()
#plt.hist(dca_for_hist, bins = 100)
#plt.show()

#plt.figure()
#dts = np.array(dts)
#dts = dts[np.abs(dts) < 100]
#plt.hist(dts, bins = 100)
#plt.show()

cut_flow_extended.Write()
output_file.Write()
output_file.Close()
