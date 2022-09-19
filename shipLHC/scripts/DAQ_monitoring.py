#!/usr/bin/env python
import ROOT,os,sys
import rootUtils as ut
import reverseMapping
from rootpyPickler import Unpickler

class DAQ_boards(ROOT.FairTask):
   " produce hitmaps as function of boards"
   def Init(self,options,monitor):
       self.M = monitor
       h = self.M.h
       ut.bookHist(h,'scifiboard','scifi hits per board; n',100,-0.5,99.5,100,0.5,500.)
       ut.bookHist(h,'mufiboard','mufi hits per board; n',100,-0.5,99.5,100,0.5,500.)
       ut.bookHist(h,'scifiboard0','scifi hits per board, unbiased; n',100,-0.5,99.5,100,0.5,500.)
       ut.bookHist(h,'mufiboard0','mufi hits per board, unbiased; n',100,-0.5,99.5,100,0.5,500.)
       self.R = reverseMapping.reversChannelMapping()
       runNr = str(options.runNumber).zfill(6)
       if options.online:
          self.R.Init(options.server+options.path+'run_'+ runNr+'/')
       else:
          self.R.Init(options.server+options.path.replace("convertedData","raw_data")+"/data/run_"+runNr+'/')
   def ExecuteEvent(self,event):
       h = self.M.h
       mult = {'scifi': [0]*100,'mufi': [0]*100}
       for aHit in event.Digi_ScifiHits:
          daq = self.R.daqChannel(aHit)
          boardN = int(daq["board"].split('_')[1])
          mult['scifi'][boardN] += 1
       for aHit in event.Digi_MuFilterHits:
          allChannels = self.M.map2Dict(aHit,'GetAllSignals')
          for c in allChannels:
              daq = self.R.daqChannel(aHit,c)
              boardN = int(daq["board"][0].split('_')[1])
              mult['mufi'][boardN] += 1
       # check for a scifi board with more than 10hits to have unbiased info for others
       bTriggered = []
       for b in range(len(mult['scifi'])):
          if mult['scifi'][b]>10:
               bTriggered.append(b)
       trigger = -1
       if len(bTriggered)>0: 
            trigger = bTriggered[int(ROOT.gRandom.Rndm()*len(bTriggered))]
       for x in mult:
            for b in range(len(mult[x])):
               rc = h[x+'board'].Fill(b,mult[x][b])
               if trigger >0 and not b==trigger:  rc = h[x+'board0'].Fill(b,mult[x][b])

   def Plot(self):
       h = self.M.h
       ut.bookCanvas(h,'boardmaps',' ',1024,768,2,3)
       h['boardmaps'].cd(1)
       h['scifiboard'].Draw('colz')
       h['boardmaps'].cd(2)
       h['mufiboard'].Draw('colz')
       h['boardmaps'].cd(3)
       h['scifiboard0'].Draw('colz')
       h['boardmaps'].cd(4)
       h['mufiboard0'].Draw('colz')
       h['boardmaps'].cd(5)
       h['scifiboard0'].ProfileX().Draw()
       h['boardmaps'].cd(6)
       h['mufiboard0'].ProfileX().Draw()
       self.M.myPrint(h['boardmaps'],'boardmaps',subdir='daq')

class Time_evolution(ROOT.FairTask):
   " time evolution of run"
   def Init(self,options,monitor):
       self.M = monitor
       h = self.M.h
       self.gtime = {'all':{0:[],1:[],3:[]}, 'B1only':{0:[],1:[],3:[]}, 'B2noB1':{0:[],1:[],3:[]}, 'noBeam':{0:[],1:[],3:[]}}
       self.QDCtime = {0:ROOT.TGraph(),1:ROOT.TGraph(),2:ROOT.TGraph(),3:ROOT.TGraph()}

       # 8*2*2*8 + 10*5*2*8 + 60*3*2 + 60*4
       self.offsets = {1:[0,8*16,16],2:[8*2*2*8,10*16,16],3:[8*2*2*8+10*5*2*8,60*2,2],4:[8*2*2*8+10*5*2*8+3*60*2,60,1]}
       ut.bookHist(h,'ctime','delta event time per channel; dt [s]'  ,1000,0.0,10.,1700,-0.5,1699.5)
       ut.bookHist(h,'ctimeZ','delta event time per channel; dt [us]',10000,0.0,100.,1700,-0.5,1699.5)
       ut.bookHist(h,'ctimeM','delta event time per channel; dt [ms]',1000,0.0,10.,1700,-0.5,1699.5)
       ut.bookHist(h,'btime','delta timestamp per channel; ',3564*4+200,-0.5,3564*4-0.5+200,1700,-0.5,1699.5)
       ut.bookHist(h,'bnr','bunch number; ',3564,-0.5,3564-0.5)
       ut.bookHist(h,'bnrF','bunch number forward tracks; ',3564,-0.5,3564-0.5)
       ut.bookHist(h,'bnrB','bunch number backward tracks; ',3564,-0.5,3564-0.5)
# type of crossing, check for b1only,b2nob1,nobeam
       self.xing = {'all':True,'B1only':False,'B2noB1':False,'noBeam':False}
       for x in self.xing:
           ut.bookHist(h,'bnr'+x,'bunch number; ',3564,-0.5,3564-0.5)

       ut.bookHist(h,'trackDir','track direction;',300,-0.5,0.25)
       ut.bookHist(h,'trackDirSig','track direction significance;',100,-20,10)
# get filling scheme
       fg  = ROOT.TFile.Open(options.server+options.path+'FSdict.root')
       pkl = Unpickler(fg)
       FSdict = pkl.load('FSdict')
       fg.Close()
       if options.runNumber in FSdict: self.fsdict = FSdict[options.runNumber]
       else:  self.fsdict = False

       self.boardsVsTime = {}
                       
       self.Nevent = -1
       self.Tprev = [-1]*1700
   def ExecuteEvent(self,event):
       xing = self.xing
       self.Nevent +=1
       h = self.M.h
       T   = event.EventHeader.GetEventTime()
       xing = {'all':True,'B1only':False,'B2noB1':False,'noBeam':False}
       if self.fsdict:
             bunchNumber = int(T%(4*3564)/4+0.5)
             nb1 = (3564 + bunchNumber - self.fsdict['phaseShift1'])%3564
             nb2 = (3564 + bunchNumber - self.fsdict['phaseShift1']- self.fsdict['phaseShift2'])%3564
             b1 = nb1 in self.fsdict['B1']
             b2 = nb2 in self.fsdict['B2']
             IP1 = False
             IP2 = False
             if b1:
                IP1 =  self.fsdict['B1'][nb1]['IP1']
             if b2:
                IP2 =  self.fsdict['B2'][nb2]['IP2']
             xing['B1only']   = b1 and not IP1 and not b2
             xing['B2noB1']  = b2 and not b1
             xing['noBeam'] = not b1 and not b2

       if xing['B1only']  and xing['B2noB1']  or xing['B1only'] and xing['noBeam'] : print('error with b1only assignment',xing)
       if xing['B2noB1']  and xing['noBeam'] : print('error with b2nob1 assignment',xing)
       Tsec = int(T/self.M.freq)

       trackTask = self.M.FairTasks['simpleTracking']
       direction = 0
       DStrack = False
       SFtrack = False
       for theTrack in self.M.Reco_MuonTracks:
            if not theTrack.getFitStatus().isFitConverged() and theTrack.GetUniqueID()==1: continue
            if theTrack.GetUniqueID()!=1: 
                DStrack = True
                continue
            SL = trackTask.trackDir(theTrack)
            if not SL: continue
            SFtrack = True
            rc = h['trackDir'].Fill(SL[0])
            rc = h['trackDirSig'].Fill(SL[1])
            if abs(SL[0])<0.03:  direction = 1
            elif SL[0]<-0.07:     direction = -1
       bn = (T%(4*3564))/4
       rc = h['bnr'].Fill( bn )
       for x in xing:
            if xing[x]: 
                 if x=='all' or (DStrack or SFtrack):  rc = h['bnr'+x].Fill( bn )

       if direction >0: rc = h['bnrF'].Fill( (T%(4*3564))/4)
       elif direction <0: rc = h['bnrB'].Fill( (T%(4*3564))/4)
       for x in xing:
          if xing[x]:
             self.gtime[x][0].append(T/self.M.freq)
             if DStrack: self.gtime[x][3].append(T/self.M.freq)
             if SFtrack: self.gtime[x][1].append(T/self.M.freq)

       qdc = {0:0,1:0,2:0,3:0}
       
       for aHit in event.Digi_MuFilterHits:
          if not aHit.isValid(): continue
          detID = aHit.GetDetectorID()
          s  = detID//10000
          p = (detID//1000)%10
          b = detID%1000
          allChannels = self.M.map2Dict(aHit,'GetAllSignals')
          for c in allChannels:
             qdc[s]+= allChannels[c]
             if b<60:                  cNr = self.offsets[s][0] + self.offsets[s][1]*p + self.offsets[s][2]*b + c 
             else:                     cNr = self.offsets[s+1][0] + self.offsets[s+1][1]*p + self.offsets[s+1][2]*(b-60) + c 
             if self.Tprev[cNr]>0:
                dT = (T - self.Tprev[cNr])/self.M.freq
                if dT<5E-9: print('something wrong',self.Nevent,s,p,b,c,dT,T,self.Tprev[cNr])
                rc = h['ctimeZ'].Fill(dT*1E6,cNr)
                rc = h['btime'].Fill(T-self.Tprev[cNr],cNr)
                rc = h['ctimeM'].Fill(dT*1E3,cNr)
                rc = h['ctime'].Fill(dT,cNr)
             nb = aHit.GetBoardID(c)
             if not nb in self.boardsVsTime: self.boardsVsTime[nb]={}
             if not Tsec in self.boardsVsTime[nb]: self.boardsVsTime[nb][Tsec]=0
             self.boardsVsTime[nb][Tsec]+=1
             self.Tprev[cNr] = T
       for aHit in event.Digi_ScifiHits:
          if not aHit.isValid(): continue
          qdc[0]+=1   
          nb = aHit.GetBoardID(0)
          if not nb in self.boardsVsTime: self.boardsVsTime[nb]={}
          if not Tsec in self.boardsVsTime[nb]: self.boardsVsTime[nb][Tsec]=0
          self.boardsVsTime[nb][Tsec]+=1
       for s in range(4):
          self.QDCtime[s].SetPoint(self.Nevent,self.Nevent,qdc[s])

   def Plot(self):
       h = self.M.h
       gtime = self.gtime
       if self.M.TStart < 0:
          T0       = gtime['all'][0][0]
          tmax   = gtime['all'][0][len(gtime['all'][0])-1] - T0
       else:
           T0 = self.M.TStart/self.M.freq
           tmax = self.M.TEnd/self.M.freq - T0

       nbins  = int(tmax)
       yunit = "events per s"
       systems = {0:'Scifi',1:'Veto',2:'US',3:'DS'}
       if self.fsdict:  Xi = ['','B1only', 'B2noB1','noBeam']
       else: Xi = ['']

       if 'time' in h:
          if self.M.TStart < 0:
             for x in Xi:
                for y in ['time','timeWt','timeWtDS']: 
                   h.pop(y+x).Delete()
             h.pop('boardVStime').Delete()
       if not 'time' in h:
         for x in Xi:
                ut.bookHist(h,'time'+x,'elapsed time from start '+x+'; t [s];'+yunit,nbins,0,tmax)
                ut.bookHist(h,'timeWt'+x,'events with Scifi(red) DS(cyan) tracks; elapsed time from start t [s];'+yunit,nbins,0,tmax)
                ut.bookHist(h,'timeWtDS'+x,'elapsed time from start, events with DS tracks; t [s];'+yunit,nbins,0,tmax)
# time evolution of boards
         ut.bookHist(h,'boardVStime','board vs time; t [s];'+yunit,nbins,0,tmax,len(self.boardsVsTime),0.5,len(self.boardsVsTime)+0.5)
       if not 'Etime' in h:
         ut.bookHist(h,'Etime','delta event time; dt [s]',100,0.0,1.)
         ut.bookHist(h,'EtimeZ','delta event time; dt [us]',10000,0.0,100.)

       for n in range(1,len(gtime['all'][0])):
           dT = gtime['all'][0][n]-gtime['all'][0][n-1]
           rc = h['Etime'].Fill( dT )
           rc = h['EtimeZ'].Fill( dT*1E6)
       K = {0:'',1:'Wt',3:'WtDS'}
       for k in K:
         for n in range(len(gtime['all'][k])):
             rc = h['time'+K[k]].Fill(gtime['all'][k][n]-T0)
       if self.fsdict:
            for x in ['B1only', 'B2noB1','noBeam']:
               for k in K:
                  for n in range(len(gtime[x][k])):
                     rc = h['time'+K[k]+x].Fill(gtime[x][k][n]-T0)

       boards = list(self.boardsVsTime.keys())
       boards.sort()
       i = 1
       yAx = h['boardVStime'].GetYaxis()
       for nb in boards:
          snb = str(nb)
          yAx.SetBinLabel(i,snb)
          for t in self.boardsVsTime[nb]:
             rc = h['boardVStime'].Fill(t-T0,i,self.boardsVsTime[nb][t]) 
          i+=1     
       ut.bookCanvas(h,'bT','board nr vs time',2000,1600,1,1)
       h['bT'].cd()
       h['boardVStime'].Draw('colz')
       self.M.myPrint(h['bT'],"board nr versus time",subdir='daq')
       
# analyse splash events
       withTGraph = False
       anaSplash = False
       if anaSplash and withTGraph:
        splashBins = []
        av = h['time'].GetEntries()/nbins
        for i in range(1,nbins+1):
           B = h['time'].GetBinContent(i)
           if B>5*av: 
              tmin = h['time'].GetBinLowEdge(i)
              tmax = tmin+h['time'].GetBinWidth(i)
              if 'splash'+str(i) in h: 
                    h.pop('splash'+str(i)).Delete()
                    h.pop('Qsplash'+str(i)).Delete() 
              ut.bookHist(h,'splash'+str(i),'; t [us];events per usec',1000000,0,(tmax-tmin)*1E6)
              ut.bookHist(h,'Qsplash'+str(i),' qdc; t [us];sum qdc per usec',1000000,0,(tmax-tmin)*1E6)
              for sy in systems:
                   # if systems[sy]+'splash'+str(i) in h: h.pop(systems[sy]+'splash'+str(i)).Delete() 
                   h[systems[sy]+'splash'+str(i)] = ROOT.TGraph()
              splashBins.append( [i,tmin,tmax] )
        for n in range(1,len(gtime['all'][0])):
           T = gtime['all'][0][n-1]-T0
           for s in splashBins:
                if T>s[1] and T<s[2]: 
                     rc = h['splash'+str(s[0])].Fill((T-s[1])*1E6)
                     for sy in systems:
                           N = h[systems[sy]+'splash'+str(s[0])].GetN()
                           h[systems[sy]+'splash'+str(s[0])].SetPoint(N,(T-s[1])*1E6,self.QDCtime[sy].GetPointY(n-1))
        N = len(splashBins)
        if N>0:
          iy = int(ROOT.TMath.Sqrt(N))
          ix = N//iy
          if N>ix*iy: ix+=1
          ut.bookCanvas(h,'Tsplash','rates',1800,1200,ix,iy)
          for sy in systems: ut.bookCanvas(h,systems[sy]+'splash','qdc sum',1800,1200,ix,iy)
          n=1
          for s in splashBins:
            h['Tsplash'].cd(n)
            h['splash'+str(s[0])].Draw('hist')
            for sy in systems: 
               tc = h[systems[sy]+'splash'].cd(n)
               tc.SetLogy(1)
               hist = h['Qsplash'+str(s[0])]
               hist.Reset()
               g = h[systems[sy]+'splash'+str(s[0])]
               for i in range(g.GetN()):
                   T,QDC =  g.GetPointX(i),g.GetPointY(i)
                   rc = hist.Fill(T,QDC)
               h[systems[sy]+'Qsplash'+str(s[0])] = ROOT.TGraph()
               for i in range(hist.GetNbinsX()):
                    QDC = hist.GetBinContent(i+1)
                    if not QDC>0: continue
                    T = hist.GetBinCenter(i+1)
                    h[systems[sy]+'Qsplash'+str(s[0])].AddPoint(T,QDC)
               h['QXsplash'+str(s[0])] = h['Qsplash'+str(s[0])].RebinX(10000,'QXsplash'+str(s[0]))
               xmax = hist.GetMaximum()
               h['QXsplash'+str(s[0])].Reset()
               h['QXsplash'+str(s[0])].SetMaximum(xmax)
               h['QXsplash'+str(s[0])].DrawClone()
               h[systems[sy]+'Qsplash'+str(s[0])].SetName(systems[sy]+'Gsplash'+str(s[0]))
               h[systems[sy]+'Qsplash'+str(s[0])].Draw('Bsame')
            n+=1
          self.M.myPrint(h['Tsplash'],"Splashes",subdir='daq')   
          for sy in systems: self.M.myPrint(h[systems[sy]+'splash'],systems[sy]+" qdc sum",subdir='daq')   

       elif anaSplash:
# analyse splash events
        splashBins = []
        av = h['time'].GetEntries()/nbins
        for i in range(1,nbins+1):
           B = h['time'].GetBinContent(i)
           if B>5*av: 
              tmin = h['time'].GetBinLowEdge(i)
              tmax = tmin+h['time'].GetBinWidth(i)
              if 'splash'+str(i) in h: h.pop('splash'+str(i)).Delete()
              ut.bookHist(h,'splash'+str(i),'; t [#mus];events per #mus',1000000,0,(tmax-tmin)*1E6)
              for sy in systems:
                   if systems[sy]+'splash'+str(i) in h: h.pop(systems[sy]+'splash'+str(i)).Delete() 
                   ut.bookHist(h,systems[sy]+'splash'+str(i),systems[sy]+'sum qdc / N; t [1#mus];average sum qdc per event #mus',1000000,0,(tmax-tmin)*1E6)
              splashBins.append( [i,tmin,tmax] )
        for n in range(1,len(gtime['all'][0])):
           T = gtime['all'][0][n-1]-T0
           for s in splashBins:
                if T>s[1] and T<s[2]: 
                     rc = h['splash'+str(s[0])].Fill((T-s[1])*1E6)
                     for sy in systems: 
                           rc = h[systems[sy]+'splash'+str(s[0])].Fill((T-s[1])*1E6,self.QDCtime[sy].GetPointY(n-1))
        N = len(splashBins)
        if N>0:
          iy = int(ROOT.TMath.Sqrt(N))
          ix = N//iy
          if N>ix*iy: ix+=1
          ut.bookCanvas(h,'Tsplash','rates',1800,1200,ix,iy)
          for sy in systems: ut.bookCanvas(h,systems[sy]+'splash','qdc sum',1800,1200,ix,iy)
          n=1
          for s in splashBins:
            h['Tsplash'].cd(n)
            h['splash'+str(s[0])].Draw('hist')
            for sy in systems: 
               tc = h[systems[sy]+'splash'].cd(n)
               tc.SetLogy(1)
               h[systems[sy]+'splash'+str(s[0])].Divide(h['splash'+str(s[0])])
               h[systems[sy]+'splash'+str(s[0])].Draw('hist')
            n+=1
          self.M.myPrint(h['Tsplash'],"Splashes",subdir='daq')
          for sy in systems: self.M.myPrint(h[systems[sy]+'splash'],systems[sy]+" qdc sum",subdir='daq')   

       ut.bookCanvas(h,'T','rates',1024,3*768,1,3)
       ut.bookCanvas(h,'Txing','rates per xing type',1024,3*768,1,4)

       tc = h['T'].cd(1)
       h['time'].SetStats(0)
       h['time'].Draw()
       tc = h['T'].cd(2)
       h['timeWt'].SetStats(0)
       h['timeWt'].SetLineColor(ROOT.kRed)
       h['timeWtDS'].SetStats(0)
       h['timeWtDS'].SetLineColor(ROOT.kCyan)
       mx = max( h['timeWt'].GetBinContent(h['timeWtDS'].GetMaximumBin()),\
                         h['timeWtDS'].GetBinContent(h['timeWtDS'].GetMaximumBin()) )
       h['timeWt'].SetMaximum(mx)
       h['timeWt'].Draw()
       h['timeWtDS'].Draw('same')
       tc = h['T'].cd(3)
       tc.SetLogy(1)
       h['EtimeZ'].Draw()
       #rc = h['EtimeZ'].Fit('expo','S','',0.,250.)
       h['T'].Update()
       stats = h['EtimeZ'].FindObject('stats')
       stats.SetOptFit(1111111)
       tc = h['T'].cd(4)
       tc.SetLogy(1)
       h['Etime'].Draw()
       #rc = h['Etime'].Fit('expo','S')
       h['T'].Update()
       stats = h['Etime'].FindObject('stats')
       stats.SetOptFit(1111111)
       h['T'].Update()
       self.M.myPrint(h['T'],"Rates",subdir='daq')
       if self.fsdict:
            j = 1
            for x in ['B1only', 'B2noB1','noBeam']:
                 tc = h['Txing'].cd(j)
                 hmax = 0
                 for ty in ['','Wt','WtDS']:
                     mx = h['time'+ty+x].GetMaximum()
                     if mx>hmax: hmax = mx
                 for ty in ['','Wt','WtDS']:
                     h['time'+ty+x].SetMaximum(mx)
                     h['time'+ty+x].SetStats(0)
                 h['time'+x].Draw()
                 h['timeWt'+x].SetLineColor(ROOT.kRed)
                 h['timeWtDS'+x].SetLineColor(ROOT.kCyan)
                 h['timeWt'+x].Draw('same')
                 h['timeWtDS'+x].Draw('same')
                 j+=1
            tc = h['Txing'].cd(j)
            hmax = 0
            for x in ['B1only', 'B2noB1','noBeam']:
               mx = h['bnr'+x].GetMaximum()
               if mx>hmax: hmax = mx
            h['bnrB1only'].SetMaximum(hmax)
            h['bnrB1only'].SetStats(0)
            h['bnrB1only'].SetLineColor(ROOT.kBlue)
            h['bnrB1only'].Draw()
            h['bnrB2noB1'].SetStats(0)
            h['bnrB2noB1'].SetLineColor(ROOT.kGreen)
            h['bnrB2noB1'].Draw('same')
            h['bnrnoBeam'].SetStats(0)
            h['bnrnoBeam'].SetLineColor(ROOT.kOrange)
            h['bnrnoBeam'].Draw('same')
            h['Txing'].Update()
            self.M.myPrint(h['Txing'],"RatesXing",subdir='daq')

       ut.bookCanvas(h,'TD',' ',1024,768,2,1)
       h['TD'].cd(1)
       h['trackDir'].Draw()
       h['TD'].cd(2)
       h['trackDirSig'].Draw()
       self.M.myPrint(h['TD'],'trackdirections',subdir='daq')

       ut.bookCanvas(h,'bunchNumber','bunch nr',2048,1600,1,3)
       tc = h['bunchNumber'].cd(1)
       h['bnr'].SetStats(0)
       h['bnr'].Draw()
       tc = h['bunchNumber'].cd(2)
       h['bnrF'].SetStats(0)
       h['bnrF'].Draw()
       tc = h['bunchNumber'].cd(3)
       h['bnrB'].SetStats(0)
       h['bnrB'].Draw()
       self.M.myPrint(h['bunchNumber'],"BunchNr",subdir='daq')

       ut.bookCanvas(h,'channels',' channel dt',1024,4*768,1,4)
       tc = h['channels'].cd(1)
       h['ctimeZ'].Draw('colz')
       tc = h['channels'].cd(2)
       h['ctimeM'].Draw('colz')
       tc = h['channels'].cd(3)
       h['ctime'].Draw('colz')
       tc = h['channels'].cd(4)
       h['btime'].Draw('colz')
       self.M.myPrint(h['channels'],"mufilter channel dT",subdir='daq')





