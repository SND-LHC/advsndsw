#!/usr/bin/env python
import ROOT,os,sys
import rootUtils as ut
import reverseMapping

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
