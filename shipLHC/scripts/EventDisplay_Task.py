#!/usr/bin/env python
import ROOT, os, sys
import rootUtils as ut
import subprocess
import json
import time
from array import array

A, B = ROOT.TVector3(), ROOT.TVector3()


class twod(ROOT.FairTask):
    "2d event display X and Y projections"

    def Init(self, options, monitor):
        self.trans2local = False
        self.minSipmMult = 1
        self.detSize = {}
        self.M = monitor
        self.systems = monitor.sdict
        run = ROOT.FairRunAna.Instance()
        self.trackTask = run.GetTask("simpleTracking")
        self.trackColor = {1: ROOT.kBlue, 3: ROOT.kRed}
        self.OT = run.GetSink().GetOutTree()
        si = self.M.snd_geo.snd_geo.Scifi
        self.detSize[0] = [si.channel_width, si.channel_width, si.scifimat_z]
        mi = self.M.snd_geo.snd_geo.MuFilter
        self.detSize[1] = [mi.VetoBarX / 2, mi.VetoBarY / 2, mi.VetoBarZ / 2]
        self.detSize[2] = [
            mi.UpstreamBarX / 2,
            mi.UpstreamBarY / 2,
            mi.UpstreamBarZ / 2,
        ]
        self.detSize[3] = [
            mi.DownstreamBarX_ver / 2,
            mi.DownstreamBarY / 2,
            mi.DownstreamBarZ / 2,
        ]
        self.options = options
        h = self.M.h
        if "simpleDisplay" not in h:
            ut.bookCanvas(
                h,
                key="simpleDisplay",
                title="2d event display",
                nx=1200,
                ny=1600,
                cx=1,
                cy=2,
            )
        h["simpleDisplay"].cd(1)
        zStart = 250.0  # TI18 coordinate system
        if self.M.snd_geo.snd_geo.Floor.z == 0:
            zStart = 60.0
        ut.bookHist(
            h, "xz", "; z [cm]; x [cm]", 500, zStart, zStart + 350.0, 100, -100.0, 10.0
        )
        ut.bookHist(
            h, "yz", "; z [cm]; y [cm]", 500, zStart, zStart + 350.0, 100, -30.0, 80.0
        )

        h["xz"].SetStats(0)
        h["yz"].SetStats(0)
        self.proj = {1: "xz", 2: "yz"}
        self.Proj = {"X": 0, "Y": 1}
        self.xNodes = {"UpstreamBar", "VetoBar", "hor"}
        self.drawLogo = True
        self.drawText = True
        self.startTimeOfRun = {}

        self.nodes = {"volMuFilter_1/volFeBlockEnd_1": ROOT.kGreen - 6}
        nodes = self.nodes
        for i in range(2):
            nodes["volVeto_1/volVetoPlane_{}_{}".format(i, i)] = ROOT.kRed
            for j in range(7):
                nodes[
                    "volVeto_1/volVetoPlane_{}_{}/volVetoBar_1{}{:0>3d}".format(
                        i, i, i, j
                    )
                ] = ROOT.kRed
            nodes["volVeto_1/subVetoBox_{}".format(i)] = ROOT.kGray + 1
        for i in range(4):
            nodes["volMuFilter_1/volMuDownstreamDet_{}_{}".format(i, i + 7)] = (
                ROOT.kBlue + 1
            )
            for j in range(60):
                nodes[
                    "volMuFilter_1/volMuDownstreamDet_{}_{}/volMuDownstreamBar_ver_3{}{:0>3d}".format(
                        i, i + 7, i, j + 60
                    )
                ] = ROOT.kBlue + 1
                if i < 3:
                    nodes[
                        "volMuFilter_1/volMuDownstreamDet_{}_{}/volMuDownstreamBar_hor_3{}{:0>3d}".format(
                            i, i + 7, i, j
                        )
                    ] = ROOT.kBlue + 1
        for i in range(4):
            nodes["volMuFilter_1/subDSBox_{}".format(i + 7)] = ROOT.kGray + 1
        for i in range(5):
            nodes["volTarget_1/ScifiVolume{}_{}000000".format(i + 1, i + 1)] = (
                ROOT.kBlue + 1
            )
            nodes["volTarget_1/volWallborder_{}".format(i)] = ROOT.kGray
            nodes["volMuFilter_1/subUSBox_{}".format(i + 2)] = ROOT.kGray + 1
            nodes["volMuFilter_1/volMuUpstreamDet_{}_{}".format(i, i + 2)] = (
                ROOT.kBlue + 1
            )
            for j in range(10):
                nodes[
                    "volMuFilter_1/volMuUpstreamDet_{}_{}/volMuUpstreamBar_2{}00{}".format(
                        i, i + 2, i, j
                    )
                ] = ROOT.kBlue + 1
            nodes["volMuFilter_1/volFeBlock_{}".format(i)] = ROOT.kGreen - 6
        for i in range(7, 10):
            nodes["volMuFilter_1/volFeBlock_{}".format(i)] = ROOT.kGreen - 6
        self.passNodes = {"Block", "Wall"}
        self.xNodes = {"UpstreamBar", "VetoBar", "hor"}

        self.Enodes = {}
        Enodes = self.Enodes
        for i in range(2):
            for j in range(7):
                Enodes[
                    "volVeto_1/volVetoPlane_{}_{}/volVetoBar_1{}{:0>3d}".format(
                        i, i, i, j
                    )
                ] = ROOT.kRed
        for i in range(3):
            for j in range(60):
                Enodes[
                    "volMuFilter_1/volMuDownstreamDet_{}_{}/volMuDownstreamBar_hor_3{}{:0>3d}".format(
                        i, i + 7, i, j
                    )
                ] = ROOT.kBlue + 2
                Enodes[
                    "volMuFilter_1/volMuDownstreamDet_{}_{}/volMuDownstreamBar_ver_3{}{:0>3d}".format(
                        i, i + 7, i, j + 60
                    )
                ] = ROOT.kBlue + 2
        for j in range(60):
            Enodes[
                "volMuFilter_1/volMuDownstreamDet_3_10/volMuDownstreamBar_ver_33{:0>3d}".format(
                    j + 60
                )
            ] = ROOT.kBlue + 2
        for i in range(5):
            for j in range(10):
                Enodes[
                    "volMuFilter_1/volMuUpstreamDet_{}_{}/volMuUpstreamBar_2{}00{}".format(
                        i, i + 2, i, j
                    )
                ] = ROOT.kBlue + 2

    def ExecuteEvent(self, event):
        h = self.M.h
        proj = self.proj
        geo = self.M.snd_geo
        detSize = self.detSize
        N = self.M.EventNumber
        options = self.options

        nav = ROOT.gGeoManager.GetCurrentNavigator()
        if options.goodEvents and not self.goodEvent(event):
            return
        if options.withTrack:
            self.trackTask.ExecuteTask()
            ntracks = self.M.Reco_MuonTracks.GetEntries()
            uniqueTracks = self.M.Reco_MuonTracks  # self.cleanTracks()
            if len(uniqueTracks) < options.nTracks:
                return

            for aTrack in self.M.Reco_MuonTracks:
                mom = aTrack.getFittedState().getMom()
                pos = aTrack.getFittedState().getPos()
                if aTrack.GetUniqueID() == 3:
                    tt = "DS"
                else:
                    tt = "Scifi"
                print(
                    tt
                    + " track direction:  X %5.2Fmrad  Y %5.2Fmrad"
                    % (mom.X() / mom.Z() * 1000, mom.Y() / mom.Z() * 1000)
                )
                print("   track position:   X %5.2Fcm  Y %5.2Fcm " % (pos.X(), pos.Y()))

        digis = []
        if event.FindBranch("Digi_ScifiHits"):
            digis.append(event.Digi_ScifiHits)
        if event.FindBranch("Digi_MuFilterHits"):
            digis.append(event.Digi_MuFilterHits)
        empty = True
        for x in digis:
            if x.GetEntries() > 0:
                if empty:
                    print("event -> %i" % N)
                empty = False
        if empty:
            return
        h["hitCollectionX"] = {
            "Scifi": [0, ROOT.TGraphErrors()],
            "DS": [0, ROOT.TGraphErrors()],
        }
        h["hitCollectionY"] = {
            "Veto": [0, ROOT.TGraphErrors()],
            "Scifi": [0, ROOT.TGraphErrors()],
            "US": [0, ROOT.TGraphErrors()],
            "DS": [0, ROOT.TGraphErrors()],
        }
        h["firedChannelsX"] = {"Scifi": [0, 0, 0], "DS": [0, 0, 0]}
        h["firedChannelsY"] = {
            "Veto": [0, 0, 0, 0],
            "Scifi": [0, 0, 0],
            "US": [0, 0, 0, 0],
            "DS": [0, 0, 0, 0],
        }
        systems = {1: "Veto", 2: "US", 3: "DS", 0: "Scifi"}
        for collection in ["hitCollectionX", "hitCollectionY"]:
            for c in h[collection]:
                rc = h[collection][c][1].SetName(c)
                rc = h[collection][c][1].Set(0)

        for p in proj:
            rc = h["simpleDisplay"].cd(p)
            h[proj[p]].Draw("b")
            self.emptyNodes()
            self.drawDetectors()
        for D in digis:
            for digi in D:
                detID = digi.GetDetectorID()
                sipmMult = 1
                if digi.GetName() == "MuFilterHit":
                    system = digi.GetSystem()
                    geo.modules["MuFilter"].GetPosition(detID, A, B)
                    sipmMult = len(digi.GetAllSignals())
                    if sipmMult < self.minSipmMult and (system == 1 or system == 2):
                        continue
                else:
                    geo.modules["Scifi"].GetSiPMPosition(detID, A, B)
                    system = 0
                curPath = nav.GetPath()
                tmp = curPath.rfind("/")
                nav.cd(curPath[:tmp])
                globA, locA = (
                    array("d", [A[0], A[1], A[2]]),
                    array("d", [A[0], A[1], A[2]]),
                )
                if self.trans2local:
                    nav.MasterToLocal(globA, locA)
                Z = A[2]
                if digi.isVertical():
                    collection = "hitCollectionX"
                    Y = locA[0]
                    sY = detSize[system][0]
                else:
                    collection = "hitCollectionY"
                    Y = locA[1]
                    sY = detSize[system][1]
                c = h[collection][systems[system]]
                rc = c[1].SetPoint(c[0], Z, Y)
                rc = c[1].SetPointError(c[0], detSize[system][2], sY)
                c[0] += 1

                self.fillNode(curPath)

                if digi.isVertical():
                    F = "firedChannelsX"
                else:
                    F = "firedChannelsY"
                ns = max(1, digi.GetnSides())
                for side in range(ns):
                    for m in range(digi.GetnSiPMs()):
                        qdc = digi.GetSignal(m + side * digi.GetnSiPMs())
                        if qdc < 0 and qdc > -900:
                            h[F][systems[system]][1] += 1
                        elif not qdc < 0:
                            h[F][systems[system]][0] += 1
                            # h[F][systems[system]][2+side]+=qdc
        h["hitCollectionY"]["Scifi"][1].SetMarkerColor(ROOT.kBlue + 2)
        h["hitCollectionX"]["Scifi"][1].SetMarkerColor(ROOT.kBlue + 2)
        k = 1
        for collection in ["hitCollectionX", "hitCollectionY"]:
            h["simpleDisplay"].cd(k)
            self.drawInfo(h["simpleDisplay"], k, N, event)
            k += 1
            for c in h[collection]:
                F = collection.replace("hitCollection", "firedChannels")
                pj = collection.split("ion")[1]
                if pj == "X" or c == "Scifi":
                    print(
                        "%1s %5s %3i  +:%3i -:%3i qdc :%5.1F"
                        % (
                            pj,
                            c,
                            h[collection][c][1].GetN(),
                            h[F][c][0],
                            h[F][c][1],
                            h[F][c][2],
                        )
                    )
                else:
                    print(
                        "%1s %5s %3i  +:%3i -:%3i qdcL:%5.1F qdcR:%5.1F"
                        % (
                            pj,
                            c,
                            h[collection][c][1].GetN(),
                            h[F][c][0],
                            h[F][c][1],
                            h[F][c][2],
                            h[F][c][3],
                        )
                    )
                if h[collection][c][1].GetN() < 1:
                    continue
                if c == "Scifi":
                    h[collection][c][1].SetMarkerStyle(20)
                    h[collection][c][1].SetMarkerSize(1.5)
                    rc = h[collection][c][1].Draw("sameP")
                    h["display:" + c] = h[collection][c][1]
        if options.withTrack:
            self.addTrack()
        h["simpleDisplay"].Update()
        if options.save:
            h["simpleDisplay"].Print("event_" + "{:04d}".format(N) + ".png")
        if options.interactive:
            rc = input("hit return for next event ")
        else:
            self.M.myPrint(
                h["simpleDisplay"], "event" + str(N % 10), subdir="eventdisplay"
            )

    def Plot(self):
        if self.M.options.save:
            os.system("convert -delay 60 -loop 0 event*.png animated.gif")

    def goodEvent(self, event):
        # can be replaced by any user selection
        maxScifiOcc = 25
        minScifi = 3
        minMufi = 5
        stations = {"Scifi": {}, "Mufi": {}}
        if event.Digi_ScifiHits.GetEntries() > maxScifiOcc:
            return False
        for d in event.Digi_ScifiHits:
            stations["Scifi"][d.GetDetectorID() // 1000000] = 1
        for d in event.Digi_MuFilterHits:
            plane = d.GetDetectorID() // 1000
            stations["Mufi"][plane] = 1
        totalN = len(stations["Mufi"]) + len(stations["Scifi"])
        if len(stations["Scifi"]) > minScifi or len(stations["Mufi"]) > minMufi:
            return True
        return False

    def addTrack(self):
        h = self.M.h
        xax = h["xz"].GetXaxis()
        nTrack = 0
        for aTrack in self.M.Reco_MuonTracks:
            for p in [0, 1]:
                h["aLine" + str(nTrack * 10 + p)] = ROOT.TGraph()
            zEx = xax.GetBinCenter(1)
            mom = aTrack.getFittedState().getMom()
            pos = aTrack.getFittedState().getPos()
            lam = (zEx - pos.z()) / mom.z()
            Ex = [pos.x() + lam * mom.x(), pos.y() + lam * mom.y()]
            for p in [0, 1]:
                h["aLine" + str(nTrack * 10 + p)].SetPoint(0, zEx, Ex[p])

            for i in range(aTrack.getNumPointsWithMeasurement()):
                state = aTrack.getFittedState(i)
                pos = state.getPos()
                for p in [0, 1]:
                    h["aLine" + str(nTrack * 10 + p)].SetPoint(i + 1, pos[2], pos[p])

            zEx = xax.GetBinCenter(xax.GetLast())
            mom = aTrack.getFittedState().getMom()
            pos = aTrack.getFittedState().getPos()
            lam = (zEx - pos.z()) / mom.z()
            Ex = [pos.x() + lam * mom.x(), pos.y() + lam * mom.y()]
            for p in [0, 1]:
                h["aLine" + str(nTrack * 10 + p)].SetPoint(i + 2, zEx, Ex[p])

            for p in [0, 1]:
                tc = h["simpleDisplay"].cd(p + 1)
                h["aLine" + str(nTrack * 10 + p)].SetLineColor(
                    self.trackColor[aTrack.GetUniqueID()]
                )
                h["aLine" + str(nTrack * 10 + p)].SetLineWidth(2)
                h["aLine" + str(nTrack * 10 + p)].Draw("same")
                tc.Update()
                h["simpleDisplay"].Update()
            nTrack += 1

    def cleanTracks(self):
        listOfDetIDs = {}
        n = 0
        for aTrack in self.M.Reco_MuonTracks:
            listOfDetIDs[n] = []
            for i in range(aTrack.getNumPointsWithMeasurement()):
                M = aTrack.getPointWithMeasurement(i)
                R = M.getRawMeasurement()
                listOfDetIDs[n].append(R.getDetId())
                if R.getDetId() > 0:
                    listOfDetIDs[n].append(R.getDetId() - 1)
                listOfDetIDs[n].append(R.getDetId() + 1)
            n += 1
        uniqueTracks = []
        for n1 in range(len(listOfDetIDs)):
            unique = True
            for n2 in range(len(listOfDetIDs)):
                if n1 == n2:
                    continue
                I = set(listOfDetIDs[n1]).intersection(listOfDetIDs[n2])
                if len(I) > 0:
                    unique = False
            if unique:
                uniqueTracks.append(n1)
        if len(uniqueTracks) > 1:
            for n1 in range(len(listOfDetIDs)):
                print(listOfDetIDs[n1])
        return uniqueTracks

    def drawDetectors(self):
        h = self.M.h
        proj = self.Proj
        nodes = self.nodes
        xNodes = self.xNodes
        passNodes = self.passNodes
        nav = ROOT.gGeoManager.GetCurrentNavigator()

        for node_ in nodes:
            node = "/cave_1/Detector_0/" + node_
            for p in proj:
                if node + p not in h:
                    nav.cd(node)
                    N = nav.GetCurrentNode()
                    S = N.GetVolume().GetShape()
                    dx, dy, dz = S.GetDX(), S.GetDY(), S.GetDZ()
                    ox, oy, oz = S.GetOrigin()[0], S.GetOrigin()[1], S.GetOrigin()[2]
                    P = {}
                    M = {}
                    if p == "X" and not any(xNode in node for xNode in xNodes):
                        P["LeftBottom"] = array("d", [-dx + ox, oy, -dz + oz])
                        P["LeftTop"] = array("d", [dx + ox, oy, -dz + oz])
                        P["RightBottom"] = array("d", [-dx + ox, oy, dz + oz])
                        P["RightTop"] = array("d", [dx + ox, oy, dz + oz])
                    elif p == "Y" and "ver" not in node:
                        P["LeftBottom"] = array("d", [ox, -dy + oy, -dz + oz])
                        P["LeftTop"] = array("d", [ox, dy + oy, -dz + oz])
                        P["RightBottom"] = array("d", [ox, -dy + oy, dz + oz])
                        P["RightTop"] = array("d", [ox, dy + oy, dz + oz])
                    else:
                        continue
                    for C in P:
                        M[C] = array("d", [0, 0, 0])
                        nav.LocalToMaster(P[C], M[C])
                    h[node + p] = ROOT.TPolyLine()
                    X = h[node + p]
                    c = proj[p]
                    X.SetPoint(0, M["LeftBottom"][2], M["LeftBottom"][c])
                    X.SetPoint(1, M["LeftTop"][2], M["LeftTop"][c])
                    X.SetPoint(2, M["RightTop"][2], M["RightTop"][c])
                    X.SetPoint(3, M["RightBottom"][2], M["RightBottom"][c])
                    X.SetPoint(4, M["LeftBottom"][2], M["LeftBottom"][c])
                    X.SetLineColor(nodes[node_])
                    X.SetLineWidth(1)
                    h["simpleDisplay"].cd(c + 1)
                    if any(passNode in node for passNode in passNodes):
                        X.SetFillColorAlpha(nodes[node_], 0.5)
                        X.Draw("f&&same")
                    X.Draw("same")
                else:
                    X = h[node + p]
                    c = proj[p]
                    h["simpleDisplay"].cd(c + 1)
                    if any(passNode in node for passNode in passNodes):
                        X.Draw("f&&same")
                    X.Draw("same")

    def fillNode(self, node):  #
        h = self.M.h
        xNodes = self.xNodes
        proj = self.Proj
        color = ROOT.kBlack
        thick = 5
        for p in proj:
            if node + p in h:
                X = h[node + p]
                if "Veto" in node:
                    color = ROOT.kRed + 1
                if "Downstream" in node:
                    thick = 5
                c = proj[p]
                h["simpleDisplay"].cd(c + 1)
                X.SetFillColor(color)
                X.SetLineColor(color)
                X.SetLineWidth(thick)
                X.Draw("f&&same")
                X.Draw("same")

    def emptyNodes(self):
        Enodes = self.Enodes
        proj = self.Proj
        h = self.M.h
        for node_ in Enodes:
            node = "/cave_1/Detector_0/" + node_
            for p in proj:
                if node + p in h:
                    X = h[node + p]
                    if X.GetFillColor() != 19:
                        c = proj[p]
                        h["simpleDisplay"].cd(c + 1)
                        X.SetFillColorAlpha(19, 0)
                        X.SetLineColor(Enodes[node_])
                        X.SetLineWidth(1)
                        X.Draw("f&&same")
                        X.Draw("same")
                else:
                    notFilled = 1

    def drawInfo(self, pad, k, N, event):
        if self.drawLogo:
            padLogo = ROOT.TPad("logo", "logo", 0.1, 0.1, 0.2, 0.3)
            padLogo.SetFillStyle(4000)
            padLogo.SetFillColorAlpha(0, 0)
            padLogo.Draw()
            logo = ROOT.TImage.Open("$SNDSW_ROOT/shipLHC/Large__SND_Logo_black_cut.png")
            logo.SetConstRatio(True)
            logo.DrawText(0, 0, "SND", 98)
            padLogo.cd()
            logo.Draw()
            pad.cd(k)

        if self.drawText:
            runNumber = event.EventHeader.GetRunId()
            if self.options.online:
                time_event = time.ctime()
                timestamp_start = True
            else:
                timestamp = event.EventHeader.GetEventTime()
                timestamp_start = self.getStartTime(runNumber)
                if timestamp_start:
                    TDC2ns = 6.23768  # conversion factor from 160MHz clock to ns
                    timestamp_s = timestamp * TDC2ns * 1e-9
                    timestamp_event = int(timestamp_start + timestamp_s)
                    time_event = datetime.fromtimestamp(timestamp_event)
            padText = ROOT.TPad("info", "info", 0.19, 0.1, 0.6, 0.3)
            padText.SetFillStyle(4000)
            padText.Draw()
            padText.cd()
            textInfo = ROOT.TLatex()
            textInfo.SetTextAlign(11)
            textInfo.SetTextFont(42)
            textInfo.SetTextSize(0.15)
            textInfo.DrawLatex(0, 0.6, "SND@LHC Experiment, CERN")
            textInfo.DrawLatex(
                0, 0.4, "Run / Event: " + str(runNumber) + " / " + str(N)
            )
            if timestamp_start:
                textInfo.DrawLatex(0, 0.2, "Time (GMT): {}".format(time_event))
            pad.cd(k)

    def getStartTime(self, runNumber):
        if runNumber in self.startTimeOfRun:
            return self.startTimeOfRun[runNumber]
        runDir = "/eos/experiment/sndlhc/raw_data/commissioning/TI18/data/run_" + str(
            runNumber
        ).zfill(6)
        jname = "run_timestamps.json"
        dirlist = str(
            subprocess.check_output(
                "xrdfs " + options.server + " ls " + runDir, shell=True
            )
        )
        if not jname in dirlist:
            return False
        with client.File() as f:
            f.open(options.server + runDir + "/run_timestamps.json")
            status, jsonStr = f.read()
            f.close()
        date = json.loads(jsonStr)
        time_str = date["start_time"]
        time_obj = time.strptime(time_str, "%Y-%m-%dT%H:%M:%S")
        self.startTimeOfRun[runNumber] = time.mktime(time_obj)
        return self.startTimeOfRun[runNumber]
