import ROOT
import numpy as np

import SndlhcGeo

# Code snippet from Simona. This should go somewhere where it can be used by several different pieces of code (monitoring, analysis, etc)
import pickle
p = open("/eos/experiment/sndlhc/convertedData/commissioning/TI18/FSdict.pkl",'rb')
FSdict = pickle.load(p)

def bunchXtype(eventTime, runN):
    if runN in FSdict: 
       fsdict = FSdict[runN]      
    else: fsdict = False
    if fsdict:
             bunchNumber = int(eventTime%(4*3564)/4+0.5)
             nb1 = (3564 + bunchNumber - fsdict['phaseShift1'])%3564
             nb2 = (3564 + bunchNumber - fsdict['phaseShift1']- fsdict['phaseShift2'])%3564
             if not "B1" in fsdict: b1 = False
             else: b1 = nb1 in fsdict['B1']
             if not "B2" in fsdict: b2 = False
             else: b2 = nb2 in fsdict['B2']
             IP1 = False
             IP2 = False
             B2noB1 = False
             b1Only = False
             noBeam = False
             if b1:
                IP1 =  fsdict['B1'][nb1]['IP1']
             if b2:
                IP2 =  fsdict['B2'][nb2]['IP2']
             if b2 and not b1:
                B2noB1 = True             
             if b1 and not b2 and not IP1:
                b1Only = True
             if not b1 and not b2: noBeam = True
             return IP1, IP2, b1Only, B2noB1, noBeam
# End snippet

import argparse
parser = argparse.ArgumentParser()

parser.add_argument("-f", dest = "inputFile", required = True)
parser.add_argument("-o", dest = "outputFile", required = True)
parser.add_argument("-t", dest = "trackFile", required = True)
parser.add_argument("-g", dest = "geoFile", required = True)

args = parser.parse_args()

snd_geo = SndlhcGeo.GeoInterface(args.geoFile)
scifiDet = ROOT.gROOT.GetListOfGlobals().FindObject('Scifi')

# Set up TTrees
ch = ROOT.TChain("rawConv")
ch.Add(args.inputFile)

ch_tracks = ROOT.TChain("rawConv")
ch_tracks.Add(args.trackFile)
ch.AddFriend(ch_tracks)

# Set up cuts
cuts = []

# Event in time with IP1 bunch crossing
def eventInBunchCrossing(event) :
    IP1, IP2, b1Only, B2noB1, noBeam = bunchXtype(event.EventHeader.GetEventTime(), event.EventHeader.GetRunId())
    return IP1
#cuts.append(["Event in time with IP1 collision", eventInBunchCrossing])
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
cuts.append(["Event direction", direction])

# Event has one reconstructed DS track
def eventHasOneTrack(event) :
    if event.Reco_MuonTracks.GetEntries() >= 1 :
        return True
    else :
        return False
cuts.append(["Event has one reconstructed DS track", eventHasOneTrack])

def trackInterceptShower(event) : 
    n_ver = [0]*5
    n_hor = [0]*5
    x_sta = [0.]*5
    y_sta = [0.]*5

    a = ROOT.TVector3()
    b = ROOT.TVector3()

    for hit in event.Digi_ScifiHits :
        if not hit.isValid() :
            continue
            
        scifiDet.GetSiPMPosition(hit.GetDetectorID(), a, b)

        if hit.isVertical() :
            n_ver[hit.GetStation()-1] += 1
            x_sta[hit.GetStation()-1] += (a.X() + b.X())/2.
        else :
            n_hor[hit.GetStation()-1] += 1
            y_sta[hit.GetStation()-1] += (a.Y() + b.Y())/2.
    
    fracsum = np.cumsum(np.add(n_ver, n_hor)/(np.sum(n_ver)+np.sum(n_hor)))
    station = next(x[0] for x in enumerate(fracsum) if x[1] > 0.05)
    
    x_first_sta = None
    y_first_sta = None

    for i in range(station, 5) :
        if n_ver[i] >= 2 :
            x_first_sta = x_sta[i]/n_ver[i]
            break

    for i in range(station, 5) :
        if n_hor[i] >= 2 :
            y_first_sta = y_sta[i]/n_hor[i]
            break
    
    track = event.Reco_MuonTracks.At(0)

    track_start = track.getStart()
    track_stop = track.getStop()

    slope_X = (track_stop.X() - track_start.X())/(track_stop.Z() - track_start.Z())
    slope_Y = (track_stop.Y() - track_start.Y())/(track_stop.Z() - track_start.Z())
    
    station_z = 298.97+13*station
        
    station_vertex = ROOT.TVector3(track_start.X() + slope_X*(station_z - track_start.Z()), track_start.Y() + slope_Y*(station_z - track_start.Y()), station_z)
    try :
        if ((station_vertex.X() - x_first_sta)**2 + (station_vertex.Y() - y_first_sta)**2)**0.5 < 5 :
            return True
        else :
            return False
    except :
        print(n_ver)
        print(n_hor)
        exit()
#cuts.append(["Track intercepts first layer < 5 cm from shower center", trackInterceptShower])

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

def SciFiRMSEdge(event) :
    x_hits = []
    y_hits = []

    for hit in event.Digi_ScifiHits :
        if not hit.isValid() :
            continue
        mat = hit.GetMat()
        sipm = hit.GetSiPM()
        channel = hit.GetSiPMChan()
        
        x = channel + sipm*128 + mat*4*128

        if hit.isVertical() :
            x_hits.append(x)
        else :
            y_hits.append(x)

    x_hits = np.array(x_hits)
    y_hits = np.array(y_hits)

    x_mean = x_hits.mean()
    y_mean = y_hits.mean()

    rms_x = (np.square(x_hits-x_mean).sum()/len(x_hits))**0.5
    rms_y = (np.square(y_hits-y_mean).sum()/len(y_hits))**0.5

    x_edge_rms =  (768 - np.abs(x_mean-768))/rms_x
    y_edge_rms = (768 - np.abs(y_mean-768))/rms_y

    if x_edge_rms < 4 :
        return False
    if y_edge_rms < 4 :
        return False
    return True

#cuts.append(["SciFi hit center more than 4*RMS away from edge", SciFiRMSEdge])

def SciFiStationProd(event) :
    hits_per_station = [0]*5

    for hit in event.Digi_ScifiHits :
        if not hit.isValid() :
            continue
        station = hit.GetStation()-1
        
        hits_per_station[station] += 1

    hits_per_station = np.array(hits_per_station)
    prod_nhits_layer = np.prod(hits_per_station[hits_per_station>0]/2.)

    if prod_nhits_layer < 500 :
        return False
    return True

#cuts.append(["SciFi hits per station product >= 500", SciFiStationProd])

def nSciFiHitsWeightedByTrackProximity(event) :
    
    track = event.Reco_MuonTracks.At(0)

    track_start = track.getStart()
    track_stop = track.getStop()

    slope_X = (track_stop.X() - track_start.X())/(track_stop.Z() - track_start.Z())
    slope_Y = (track_stop.Y() - track_start.Y())/(track_stop.Z() - track_start.Z())

    a = ROOT.TVector3()
    b = ROOT.TVector3()

    weighted_n_hits_hor = 0.
    weighted_n_hits_ver = 0.


    hits_per_sta = [0]*5


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

        hits_per_sta[hit.GetStation()-1] += 1

#        print(hit.isVertical(), track_start_XY + slope_XY*(hit_pos[0]-track_start.Z()), hit_pos[1])
#        print(np.abs(track_start_XY + slope_XY*(hit_pos[0]-track_start.Z()) - hit_pos[1]))
        if hit.isVertical() :
            weighted_n_hits_ver += np.square(track_start_XY + slope_XY*(hit_pos[0]-track_start.Z()) - hit_pos[1])
        else :
            weighted_n_hits_hor += np.square(track_start_XY + slope_XY*(hit_pos[0]-track_start.Z()) - hit_pos[1])

#    print(min([weighted_n_hits_ver, weighted_n_hits_hor]), weighted_n_hits_ver, weighted_n_hits_hor)
#    print(hits_per_sta)
    return min([weighted_n_hits_ver, weighted_n_hits_hor])

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

dca_for_hist = []
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
    dca_for_hist.append(sum_min_dca(event))
    print(dca_for_hist[-1])

plt.figure()
plt.hist(dca_for_hist, bins = 100)
plt.show()

cut_flow_extended.Write()
output_file.Write()
output_file.Close()
