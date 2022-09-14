import ROOT,os
import rootUtils as ut
from urllib.request import urlopen
import numpy
import time,calendar
import pickle
import subprocess
from XRootD import client

fromElog = {4361:7902,4362: 7921, 4363: 7921, 4364: 7921, 4365: 7921, 4366: 7922, 4367: 7923, 4368: 7923, 4369: 7923, 4370: 7924, 4371: 7924, 4372: 7925, 4373: 7926, 4374: 7927, 4375: 7928, 4376: 7929, 4377: 7930, 4378: 7931, 4379: 7931, 4380: 7932, 4381: 7933, 4382: 7934, 4383: 7935, 4384: 7936, 4385: 7937, 4386: 7938, 4387: 7939, 4388: 7940, 4389: 7941, 4390: 7942, 4391: 7943, 4392: 7944, 4393: 7945, 4394: 7946, 4395: 7947, 4396: 7948, 4397: 7949, 4398: 7950, 4399: 7951, 4400: 7952, 4401: 7953, 4402: 7954, 4403: 7955, 4404: 7956, 4405: 7957, 4406: 7958, 4407: 7959, 4408: 7959, 4409: 7960, 4410: 7960, 4411: 7961, 4412: 7961, 4413: 7962, 4414: 7963, 4415: 7963, 4416: 7964, 4418: 7964, 4419: 7965, 4420: 7966, 4421: 7966, 4422: 7967, 4423: 7967, 4424: 7967, 4425: 7967, 4426: 7968, 4427: 7969, 4428: 7969, 4429: 7970, 4430: 7971, 4431: 7971, 4432: 7972, 4433: 7973, 4434: 7974, 4435: 7974, 4436: 7975, 4437: 7976, 4438: 7976, 4440: 7977, 4441: 7977, 4443: 7977, 4446: 7977, 4447: 7977, 4448: 7977, 4449: 7978, 4450: 7979, 4451: 7979, 4452: 7979, 4461: 7979, 4462: 7982, 4463: 7983, 4464: 7984, 4465: 7985, 4466: 7986, 4467: 7987, 4468: 7988, 4469: 7989, 4470: 7990, 4471: 7991, 4472: 7992, 4473: 7993, 4474: 7994, 4475: 7995, 4476: 7996, 4477: 7997, 4478: 7998, 4479: 7999, 4480: 7999, 4481: 8000, 4482: 8001, 4483: 8002, 4484: 8003, 4485: 8004, 4486: 8005, 4487: 8006, 4488: 8007, 4489: 8008, 4493: 8008, 4494: 8008, 4495: 8009, 4496: 8010, 4497: 8010, 4498: 8011, 4499: 8012, 4500: 8013, 4501: 8014, 4503: 8016, 4504: 8017, 4505: 8018, 4506: 8018, 4507: 8019, 4508: 8019, 4509: 8020, 4510: 8021, 4511: 8021, 4512: 8022, 4513: 8022, 4514: 8023, 4515: 8023, 4516: 8024, 4523: 8025, 4524: 8025, 4525: 8025, 4526: 8026, 4527: 8027, 4528: 8028, 4529: 8028, 4530: 8029, 4531: 8029, 4532: 8030, 4533: 8031, 4534: 8032, 4535: 8032, 4536: 8033, 4537: 8033, 4539: 8033, 4540: 8033, 4541: 8033, 4542: 8034, 4543: 8034, 4544: 8034, 4557: 8035, 4558: 8035, 4559: 8036, 4561: 8038, 4562: 8039, 4563: 8040, 4564: 8041, 4568: 8043, 4569: 8044, 4570: 8044, 4571: 8045, 4572: 8046, 4573: 8047, 4574: 8047, 4575: 8047, 4578: 8047, 4579: 8047, 4580: 8050, 4581: 8051, 4582: 8052, 4583: 8052, 4585: 8053, 4586: 8053, 4587: 8054, 4588: 8055, 4589: 8056, 4590: 8056, 4591: 8057, 4592: 8058, 4593: 8058, 4594: 8059, 4595: 8059, 4598: 8061, 4601: 8061, 4602: 8061, 4603: 8061, 4604: 8062, 4606: 8063, 4612: 8063, 4613: 8064, 4614: 8064, 4615: 8065, 4616: 8066, 4617: 8067, 4619: 8068, 4620: 8069, 4621: 8069, 4622: 8070, 4623: 8071, 4624: 8071, 4625: 8072, 4626: 8072, 4627: 8073, 4628: 8073, 4629: 8073, 4630: 8073, 4631: 8073, 4632: 8073, 4635: 8073, 4636: 8074, 4637: 8075, 4638: 8076, 4639: 8076, 4641: 8076, 4642: 8076, 4643: 8076, 4644: 8076, 4645: 8076, 4646: 8076, 4647: 8076, 4648: 8077, 4649: 8078, 4650: 8079, 4653: 8080, 4654: 8081, 4657: 8082, 4658: 8082, 4659: 8082, 4660: 8082, 4662: 8083, 4665: 8083, 4666: 8083, 4667: 8083, 4669: 8084, 4670: 8084, 4672: 8084, 4673: 8084, 4675: 8084, 4676: 8084, 4677: 8084, 4678: 8084, 4679: 8084, 4680: 8084, 4681: 8084, 4682: 8084, 4684: 8084, 4686: 8084, 4687: 8084, 4688: 8084, 4689: 8084, 4690: 8084, 4693: 8084, 4560: 8037, 4661: 8083}


class fillingScheme():

   def Init(self,options):
     self.options = options
     self.h={}
     self.path = options.path
     self.content = ''
     self.phaseShift1 =  3564-404  # IP1
     self.phaseShift2 =  129          # IP2
     self.FSdict = {}
     self.LumiInt = {}
     if "FSdict.pkl" in os.listdir():
         p = open("FSdict.pkl",'rb')
         self.FSdict = pickle.load(p)
         p = open("RunInfodict.pkl",'rb')
         self.runInfo = pickle.load(p)

     self.startTimes = { }
     self.date = {}
     self.runInfo = {}
   def readStartTimes(self):
       months = {'Jan.':'01-','Feb.':'02-','Mar.':'03-','Apr.':'04-','May.':'05-','Jun.':'06-','Jul.':'07-','Aug.':'08-','Sep-':'09-','Oct.':'10-','Nov.':'11-','Dec.':'12-'}
       s = open('/mnt/hgfs/microDisk/SND@LHC/StartTimes.log')
       L = s.readlines()
       for i in range(len(L)):
            l = L[i]
            if not l.find("(Info,Run)")>0: continue
            offset = 0
            if l.find('PM')>0: offset = 12*3600
            offset -= 2*3600  # stupid time zone
            for x in ['AM','PM']:
                 l = l.replace(x,'')
            tmp = l.split("(Info,Run)")
            tmp[0] = tmp[0].replace(' ','')
            for m in months: tmp[0] = tmp[0].replace(m,months[m])
            time_obj = time.strptime(tmp[0], '%m-%d,%Y-%H:%M:%S')
            r = tmp[1].rfind('Run')
            if r > 4570: continue #not needed, extracted from json file in run directory
            run = int(tmp[1][r+3:].split('Started')[0])
            self.startTimes[run] = calendar.timegm(time_obj)+offset

   def getFillNrFromElog(self):
       # starts with 4362 and ends with 4693
       # from 4663 - 4687 are not valid
       e = open('/mnt/hgfs/microDisk/SND@LHC/ELOG DAQ Test.html')

       self.fromElogX = { }
       L = e.readlines()
       for i in range(len(L)):
            l = L[i]
            k = l.find('- New Fill')
            if k<0: continue
            fillNr = int(l[k+10:].split(' ')[1])
            for j in range(i+1,len(L)):
               k = L[j].find('Started with')
               if k>0:
                   m = L[j][:k].rfind('Run')
                   runNr =  int(L[j][m+3:k])
                   self.fromElogX[runNr] = fillNr
               if L[j].find('- New Fill')>0: break

# by hand manipulations:
# 4560, only stops no start
# 4661, only stops no start
       self.fromElogX[4560]=8037
       self.fromElogX[4661]=8083

   def getFillNrFromRunNr(self,runNumber):
       if runNumber in fromElog:
            return fromElog[runNumber]
       FILL_NUMBER_MASK = 0x000000000000FFFF
       R = ROOT.TFile.Open(os.environ['EOSSHIP']+\
       "/eos/experiment/sndlhc/raw_data/commissioning/TI18/data/run_"+str(runNumber).zfill(6)+"/data_0000.root")
       try:
         if R.Get('event'):
            rc = R.event.GetEvent(0)
            flags = R.event.flags
         else:
            event = R.data
            event.GetEvent(0)
            flags = event.evt_flags
         fillNumber = numpy.bitwise_and(FILL_NUMBER_MASK,flags)
         if fillNumber<1: fillNumber = False
         print('fill number =',fillNumber )
       except:
         fillNumber = False
       return fillNumber

   def getLumiAtIP1(self,fillnr):
       try:
          with urlopen('https://lpc.web.cern.ch/cgi-bin/fillAnalysis.py?year=2022&action=fillData&exp=ATLAS&fillnr='+str(fillnr)) as webpage:
              tmp = webpage.read().decode()
       except:
          print('lumi info not avaible',fillnr)
          return -1
       exec("self.content = "+tmp)
       """ Each lumi file contains:
            time stab l dl sl dsl
            where
            time =  UNIX time, i.e. in seconds since UTC Jan 1, 1970, 00:00:00, not counting leap seconds.
            stab =  stable beam flag: a float value between 0 and 1 which corresponds to the fraction of time spent in stable beams for this time bin.
            l    =  luminosity in Hz/ub
            dl   =  point-to-point error on luminosity in Hz/ub
            sl   =  specific luminosity in Hz/ub  (see below)
            dsl  =  point-to-point error on specific luminosity in Hz/ub
            self.content['data']['fillData']['data'][0]
            [1660892004.2119756, 1.0, 9907.5947265625, 0.0, 3.2100153151839246e-22, 0.0]
            time.ctime(1660892004.2119756)   -> 'Fri Aug 19 08:53:24 2022'
       """
       self.lumiAtIP1 = {'startTime':self.content['data']['fillData']['data'][0][0],'lumiTime':ROOT.TGraph()}
       X = self.content['data']['fillData']['data']
       t0 =  self.lumiAtIP1['startTime']
       for n in range(len(X)):
            self.lumiAtIP1['lumiTime'].AddPoint(X[n][0]-t0,X[n][2])
       return 0
       
   def drawLumi(self,runNumber):
       R = ROOT.TFile.Open(os.environ['EOSSHIP']+"/eos/experiment/sndlhc/www/offline/run"+str(runNumber).zfill(6)+".root")
       ROOT.gROOT.cd()
       bCanvas = R.daq.Get('T')
       self.h['time'] = bCanvas.FindObject('time').Clone('time')
       nbins = self.h['time'].GetNbinsX()
       endTime = self.h['time'].GetBinCenter(nbins)  # in seconds
       ROOT.gROOT.cd()
       fillNumber = self.getFillNrFromRunNr(runNumber)
       rc = self.getLumiAtIP1(fillNumber)
       if rc<0: return
       # for overlay, need SND@LHC startTime to adjust with lumiTime 
       runDir = "/eos/experiment/sndlhc/raw_data/commissioning/TI18/data/run_"+str(runNumber).zfill(6)
       jname = "run_timestamps.json"
       dirlist  = str( subprocess.check_output("xrdfs "+os.environ['EOSSHIP']+" ls "+runDir,shell=True) )
       startTime=0
       if jname in dirlist:
         with client.File() as f:
            f.open(os.environ['EOSSHIP']+runDir+"/run_timestamps.json")
            status, jsonStr = f.read()
         exec("self.date = "+jsonStr.decode())
         time_str = self.date['start_time'].replace('Z','')
         time_obj = time.strptime(time_str, '%Y-%m-%dT%H:%M:%S')
         self.startTime = calendar.timegm(time_obj)
       else:  # try reading ecs log file
         if len(self.startTimes)==0: self.readStartTimes()
         if runNumber in self.startTimes:
             self.startTime = self.startTimes[runNumber]
         else: return
       self.date['start_time'] = time.ctime(self.startTime)
       self.lumiAtlas = ROOT.TGraph()
       deltaT = self.lumiAtIP1['startTime']  - self.startTime # account for timezone/summertime
       self.Lmax = 0
       self.Lint = 0
       self.Lsnd = 0
       tprev = [-1,0]
       for n in range(self.lumiAtIP1['lumiTime'].GetN()):
               t = self.lumiAtIP1['lumiTime'].GetPointX(n) + deltaT
               l = self.lumiAtIP1['lumiTime'].GetPointY(n)
               self.lumiAtlas.AddPoint(t,l)
               if l>self.Lmax: self.Lmax = l
               if tprev[0] < 0: 
                   tprev = [t,l]
               else:
                   dt = t - tprev[0]
                   X =  (tprev[1]+l)/2.
                   self.Lint +=X*dt
                   if t>0 and t<endTime: self.Lsnd += X*dt
                   tprev[0] = t
                   tprev[1] = l
       self.LumiInt[runNumber] = [self.Lint,self.Lsnd]
       self.h['c1'].cd()
       self.h['time10'] = self.h['time'].Clone('time10')
       self.h['time10'].Rebin(10)
       self.h['time10'].SetMinimum(0)
       self.h['time10'].Scale(self.options.postScale/10.)
       if self.date['start_time'].find('Tu')<0 and self.date['start_time'].find('Th')<0:
           tmp = self.date['start_time'].replace('T',' ').replace('Z','')
       else: tmp = self.date['start_time']
       self.h['time10'].SetTitle('Run '+str(runNumber)+'  Fill '+str(fillNumber)+'  '+tmp)
# what to do with spikes? 
       mx = self.h['time10'].GetMaximumBin()
       side = self.h['time10'].GetBinContent(mx-1)+self.h['time10'].GetBinContent(mx+1)
       if self.h['time10'].GetBinContent(mx+1)<0.2*self.h['time10'].GetBinContent(mx-1): side = 2*self.h['time10'].GetBinContent(mx-1) # happens at end of fill
       newMx = max(10,0.75*side)
# special runs
       if runNumber==4423: newMx = 100
       if runNumber==4415: newMx = 150
       if runNumber==4362: newMx = 7
       if self.h['time10'].GetBinContent(mx) > side: self.h['time10'].SetMaximum(newMx)
       self.h['time10'].Draw('hist')
       self.h['c1'].Update()
# work with second axis
       rightmax = 1.1*self.Lmax/1000
       self.scale = ROOT.gPad.GetUymax()/rightmax
       self.lumiAtlasS = ROOT.TGraph()
       self.lumiAtlasS.SetLineColor(ROOT.kMagenta)
       self.lumiAtlasS.SetLineWidth(3)

       for n in range(self.lumiAtIP1['lumiTime'].GetN()):
             t = self.lumiAtIP1['lumiTime'].GetPointX(n) + deltaT
             l = self.lumiAtIP1['lumiTime'].GetPointY(n)
             self.lumiAtlasS.AddPoint(t,l*self.scale/1000)
       self.lumiAtlasS.Draw('same')
       self.h['ax1'] = ROOT.TGaxis(ROOT.gPad.GetUxmax(), ROOT.gPad.GetUymin(),
                   ROOT.gPad.GetUxmax(), ROOT.gPad.GetUymax(),
                   0, rightmax, 510, "+L")
       l = self.Lsnd/1E9
       ul = 'fb'
       if l < 0.01:
           l = self.Lsnd/1E9
           ul = 'pb'
       self.h['ax1'].SetTitle('L [Hz/nb]   Integral '+"%5.2F %s^{-1}"%(l,ul))
       self.h['ax1'].SetTextFont(42)
       self.h['ax1'].SetLabelFont(42)
       self.h['ax1'].SetTextColor(ROOT.kMagenta)
       self.h['ax1'].Draw()
       self.h['time10'].Draw('histsame')
       self.h['c1'].Print(options.path+'Lumi-run'+str(runNumber).zfill(6)+'.root')
       self.h['c1'].Print(options.path+'Lumi-run'+str(runNumber).zfill(6)+'.pdf')
       self.h['c1'].Print(options.path+'Lumi-run'+str(runNumber).zfill(6)+'.png')

   def extractFillingScheme(self,fillNr):
       with urlopen('https://lpc.web.cern.ch/cgi-bin/schemeInfo.py?fill='+fillNr+'&fmt=json') as webpage:
           tmp = webpage.read().decode()
       exec("self.content = "+tmp)
       if len(self.content['fills']) < 1: 
              print('Filling scheme not yet known',fillNr,self.options.runNumbers)
              return -1
       csv = self.content['fills'][fillNr]['csv'].split('\n')
       nB1 = csv.index('B1 bucket number,IP1,IP2,IP5,IP8')
       F = ROOT.TFile(self.path+'fillingScheme-'+fillNr+'.root','recreate')
       nt = ROOT.TNtuple('fill'+fillNr,'b1 IP1 IP2','B1:IP1:IP2:IsB2')
       while nB1>0:
           tmp = csv[nB1+1].split(',')
           if len(tmp)!=5: break
           nB1+=1
           rc = nt.Fill(int(tmp[0]),int(tmp[1].replace('-','-1')),int(tmp[2].replace('-','-1')),0)
       nB2 = csv.index('B2 bucket number,IP1,IP2,IP5,IP8')
       while nB2>0:
           tmp = csv[nB2+1].split(',')
           if len(tmp)!=5: break
           nB2+=1
           rc = nt.Fill(int(tmp[0]),int(tmp[1].replace('-','-1')),int(tmp[2].replace('-','-1')),1)
       nt.Write()
       F.Close()
       return 0

   def extractPhaseShift(self,fillNr,runNumber):
         self.F = ROOT.TFile(self.path+'fillingScheme-'+fillNr+'.root')
         self.fs = self.F.Get('fill'+fillNr)
# convert to dictionary
         self.FSdict[runNumber] = {"fillNumber":fillNr,'phaseShift1':0,'phaseShift2':0,"B1":{},"B2":{}}
         fsdict = self.FSdict[runNumber]
         for x in self.fs:
              if x.IsB2>0:
                   fsdict['B2'][(x.B1-1)/10]={'IP1':x.IP1>0,'IP2':x.IP2>0}
              else:
                   fsdict['B1'][(x.B1-1)/10]={'IP1':x.IP1>0,'IP2':x.IP2>0}
         R = ROOT.TFile.Open(os.environ['EOSSHIP']+"/eos/experiment/sndlhc/www/offline/run"+str(runNumber).zfill(6)+".root")
         ROOT.gROOT.cd()
         self.h['bnr'] = R.daq.Get('bunchNumber').FindObject('bnr').Clone('bnr')
         R.Close()
         self.matches = {}
         for phase1 in range(0,3564):
               self.matches[phase1]=0
               for n in range(0,3564):
                   if not n in fsdict['B1']: continue
                   j = (n+phase1)%3564 + 1
                   if fsdict['B1'][n]['IP1']: self.matches[phase1]+=self.h['bnr'].GetBinContent(j)
         self.phaseShift1 = max(self.matches,key=self.matches.get)
         print('phaseShift1 found:',self.phaseShift1,3564-self.phaseShift1)
         self.matches = {}
         for phase2 in range(0,3564):
               self.matches[phase2]=0
               for n in range(0,3564):
                   if not n in fsdict['B2']: continue
                   j = (n+self.phaseShift1+phase2)%3564 + 1    # bin number
                   ip1 = (j-1+3564-self.phaseShift1)%3564
                   # take only bins which are not associated to collisions in IP1
                   if ip1 in fsdict['B1']:
                       if fsdict['B1'][ip1]['IP1']: continue
                   if fsdict['B2'][n]['IP2'] or 1>0: 
                      self.matches[phase2]+=self.h['bnr'].GetBinContent(j)
         self.phaseShift2 = max(self.matches,key=self.matches.get)
         print('phaseShift2 found:',self.phaseShift2,3564-self.phaseShift2)
         fsdict['phaseShift1'] = self.phaseShift1
         fsdict['phaseShift2'] = self.phaseShift2
         if not (3564-self.phaseShift2) == 129:
            print('There is a problem with phaseshift2 for run',runNumber,3564-self.phaseShift2)
         fsdict['phaseShift2'] = 3564 - 129
         self.phaseShift2 = fsdict['phaseShift2']

   def plotBunchStructure(self,fillNr,runNumber):
         h=self.h
         self.F = ROOT.TFile(self.path+'fillingScheme-'+fillNr+'.root')
         self.fs = self.F.Get('fill'+fillNr)
         ut.bookHist(h,'b1','b1',35640,-0.5,35639.5)
         ut.bookHist(h,'IP1','IP1',35640,-0.5,35639.5)
         ut.bookHist(h,'IP2','IP2',35640,-0.5,35639.5)
         ut.bookHist(h,'b1z','b1',3564,-0.5,3563.5)
         ut.bookHist(h,'b2z','b2',3564,-0.5,3563.5)
         ut.bookHist(h,'IP1z','IP1',3564,-0.5,3563.5)
         ut.bookHist(h,'IP2z','IP2',3564,-0.5,3563.5)
         h['b1'].Draw()
         h['b1z'].SetLineColor(ROOT.kBlue)
         h['b2z'].SetLineColor(ROOT.kCyan)
         h['IP1z'].SetLineColor(ROOT.kRed)
         h['IP2z'].SetLineColor(ROOT.kOrange)
         h['b1z'].SetStats(0)
         h['IP1z'].SetStats(0)
         h['IP2z'].SetStats(0)

         R = ROOT.TFile.Open(os.environ['EOSSHIP']+"/eos/experiment/sndlhc/www/offline/run"+str(runNumber).zfill(6)+".root")
         ROOT.gROOT.cd()
         bCanvas = R.daq.Get('bunchNumber')
         h['bnr']= bCanvas.FindObject('bnr').Clone('bnr')
         ROOT.gROOT.cd()

         self.Draw()

   def Draw(self):
         h = self.h
         h['c1'].cd()
         self.fs.Draw('(( (B1-1)/10+'+str(self.phaseShift1)+')%3564)>>b1z','!(IsB2>0)','hist')
         self.fs.Draw('(( (B1-1)/10+'+str(self.phaseShift1)+')%3564)>>IP1z','IP1>-0.6&&(!(IsB2>0))','hist')
         self.fs.Draw('(( (B1-1)/10+'+str(self.phaseShift1+ self.phaseShift2)+')%3564)>>IP2z','IP2>-0.6&&IsB2>0','hist')
         self.fs.Draw('(( (B1-1)/10+'+str(self.phaseShift1+ self.phaseShift2)+')%3564)>>b2z','IsB2>0','hist')
         norm = h['bnr'].GetBinContent(h['bnr'].GetMaximumBin())
         h['b1z'].Scale(norm*1.5)
         h['IP1z'].Scale(norm*1.0)
         h['b2z'].Scale(norm*0.5)
         h['IP2z'].Scale(norm*0.3)
         h['bnr'].SetStats(0)
         h['bnr'].SetLineColor(ROOT.kBlack)
         h['b1z'].Draw('hist')
         h['bnr'].Draw('histsame')
         h['b1z'].Draw('histsame')
         txt = 'phase shift B1, B2: '+str(3564-self.phaseShift1)+','+str(3564-self.phaseShift2)+' for run '+str(options.runNumbers)
         txt += " fill nr "+options.fillNumbers
         h['b1z'].SetTitle(txt)
         if self.options.withIP2: 
               h['IP2z'].Draw('histsame')
               h['b2z'].Draw('histsame')
         h['IP1z'].Draw('histsame')

   def Extract(self):
        if self.options.fillNumbers=='':
           fillNumber = self.getFillNrFromRunNr(int(options.runNumbers))
           if not fillNumber: 
               print('Fill number not found')
           else:
               rc = self.extractFillingScheme(str(fillNumber))
               if not rc<0:
                 self.options.fillNumbers = str(fillNumber)
                 self.extractPhaseShift(self.options.fillNumbers,self.options.runNumbers)
                 r = int(self.options.runNumbers)
                 self.plotBunchStructure(self.options.fillNumbers,r)
                 self.h['c1'].Print(self.options.path+'FS-run'+str(r).zfill(6)+'.root')
                 self.h['c1'].Print(self.options.path+'FS-run'+str(r).zfill(6)+'.pdf')
                 self.h['c1'].Print(self.options.path+'FS-run'+str(r).zfill(6)+'.png')
        else:
           for r in options.fillNumbers.split(','):
              self.extractFillingScheme(r)

   def test(self,runnr,I=True):
       h = self.h
       p = open("FSdict.pkl",'rb')
       self.FSdict = pickle.load(p)
       if runnr in self.FSdict: fsdict = self.FSdict[runnr]
       else: 
            print('run number not known. Exit.')
            return()
       options.runNumbers = str(runnr)
       fillNr = fsdict['fillNumber']
       if I:
          f=open('fill'+str(fillNr)+'.csv')
          X = f.readlines()
          keys = {'A6R4.B1':0,'A6R4.B2':0,'B6R4.B1':0,'B6R4.B2':0}
          for k in keys:
             n=0
             for l in X:
                if l.find(k)>0: keys[k]=n
                n+=1
          print(keys)
          A = X[keys['A6R4.B2']//2].replace('\n','').split(',')
          print(len(A))
          for bunchNumber in range(1,3565):
             b1 = (bunchNumber-1) in fsdict['B1']
             if b1 and float(A[bunchNumber])>1: continue
             if not b1 and float(A[bunchNumber])<1: continue
             print(bunchNumber,A[bunchNumber],b1)
          return
       FS.phaseShift1 = fsdict['phaseShift1']
       self.plotBunchStructure(fsdict['fillNumber'],runnr)
       w = {}
       for x in ['b1z','IP1z','b2z','IP2z']:
          w[x]=h[x].GetMaximum()
          h[x].Reset()
       for bunchNumber in range(0,3564):
             nb1 = ( 3564+bunchNumber - fsdict['phaseShift1'])%3564
             if nb1 in fsdict['B1']:
                 rc = h['b1z'].Fill(bunchNumber,w['b1z'])
                 if fsdict['B1'][nb1]['IP1']: rc = h['IP1z'].Fill(bunchNumber,w['IP1z'])
             nb2 = ( 3564 + bunchNumber - fsdict['phaseShift1'] - fsdict['phaseShift2'])%3564
             if nb2 in fsdict['B2']:
                 rc = h['b2z'].Fill(bunchNumber,w['b2z'])
                 if fsdict['B2'][nb2]['IP2']: rc = h['IP2z'].Fill(bunchNumber,w['IP2z'])

   def b1b2(self,runnr,b):
     if runnr in self.FSdict: 
            fsdict = self.FSdict[runnr]
            nb1 = ( 3564 + b - fsdict['phaseShift1'])%3564
            nb2 = ( 3564 + b - fsdict['phaseShift1'] - fsdict['phaseShift2'])%3564
            print('b1 bunch number',nb1,nb2)

   def merge(self):
        h = self.h
        for fname in os.listdir():
            if fname.find('FS')==0 and fname.find('.root')>0:
                rname = fname.split('-')[1].split('.')[0]
                F = ROOT.TFile(fname)
                h[rname] = F.c1.Clone(rname)
                h[rname].SetName(rname)
                h[rname].SetTitle(rname)
        F = ROOT.TFile('BunchStructure.root','recreate')
        keys = list(h.keys())
        keys.sort(reverse=True)
        for r in keys:
           if r.find('run')==0: h[r].Write()

   def mergeLumi(self):
        h = self.h
        for fname in os.listdir():
            if fname.find('Lumi')==0 and fname.find('.root')>0:
                if fname=="Lumi.root" : continue
                rname = fname.split('-')[1].split('.')[0]
                F = ROOT.TFile(fname)
                h[rname] = F.c1.Clone(rname)
                h[rname].SetName(rname)
                h[rname].SetTitle(rname)
        F = ROOT.TFile('Lumi.root','recreate')
        keys = list(h.keys())
        keys.sort(reverse=True)
        for r in keys:
           if r.find('run')==0: h[r].Write()

   def lhcNumbering(self):
        h = self.h
        F = ROOT.TFile('BunchStructure.root')
        p = open("FSdict.pkl",'rb')
        self.FSdict = pickle.load(p)
        for k in F.GetListOfKeys():
                rname = k.GetName()
                newname = 'lhc-'+rname
                h[newname] = F.Get(rname).Clone(newname)
                h[newname].SetName(rname)
                h[newname].SetTitle(rname)
        F.Close()
        F = ROOT.TFile('BunchStructureLHC.root','recreate')
        for newname in h:
                if not newname.find('lhc-') == 0:continue
                X = {}
                for p in h[newname].GetListOfPrimitives():
                   X[p.GetName()] = p
                a = X['b1z'].GetTitle().replace(' ','')
                X['bnr'].SetLineColor(ROOT.kBlack)
                runnr = int(a.split('run')[1].split('fill')[0])
                phaseShift1 = self.FSdict[runnr]['phaseShift1']
                phaseShift2 = self.FSdict[runnr]['phaseShift2']
                for hname in X:
                   histo=X[hname]
                   if not hasattr(histo,'GetBinContent'): continue
                   tmp = {}
                   for i in range(3564):
                      tmp[i+1] = histo.GetBinContent(i+1)
                   for i in range(3564):
                        newBin = (3564 - phaseShift1 + i)%3564
                        histo.SetBinContent(newBin,tmp[i+1])
                h[newname].Update()
                h[newname].Write()

   def getEntriesPerRun(self,r):
# check for partitions
          runNr = str(r).zfill(6)
          partitions = []
          path = "/eos/experiment/sndlhc/convertedData/commissioning/TI18/"
          dirlist  = str( subprocess.check_output("xrdfs "+os.environ['EOSSHIP']+" ls "+path+"run_"+runNr,shell=True) )
          for x in dirlist.split('\\n'):
             ix = x.find('sndsw_raw-')
             if ix<0: continue
             partitions.append(x[ix:])
          eventChain = ROOT.TChain('rawConv')
          for p in partitions:
             eventChain.Add(os.environ['EOSSHIP']+path+'run_'+runNr+'/'+p)
          return eventChain.GetEntries()

   def getTotalStat(self):
      L = 0
      N = 0
      for r in self.runInfo:
           L+=self.runInfo[r]['lumiAtIP1withSNDLHC']
           N+=self.runInfo[r]['Entries']

   def makeLatex(self):
        # make latex code 
        # filling schemes:  FS-run004626.pdf  and  Lumi-run004626.pdf
        L,N= 0,0
        for r in self.runInfo:
            L+=self.runInfo[r]['lumiAtIP1withSNDLHC']
            N+=self.runInfo[r]['Entries']
        lines = []
        lines.append("\documentclass{beamer}")
        lines.append("\mode<presentation>")
        lines.append("{\\usetheme{Singapore}}")
        lines.append("\\usepackage{graphicx}")
        lines.append("\\usepackage[space]{grffile}")
        lines.append("\\usepackage[english]{babel}")
        lines.append("\\usepackage[latin1]{inputenc}")
        lines.append("\\usepackage[T1]{fontenc}")
        lines.append("\\title[Short Paper Title] % (optional, use only with long paper titles)")
        lines.append("{SND@LHC Run Summary July - August 2022}")
        lines.append("\date[Short Occasion] % (optional)")
        lines.append("{ 1 September 2022}")
        lines.append("\\begin{document}")
        lines.append("\\begin{frame}{}")
        lines.append("1 September 2022")
        lines.append("\\newline  ")
        lines.append("\\newline  ")
        lines.append("Run Summary for July - August 2022")
        nTXT = "$%5.2F\\times 10^9 $"%(N/1E9)
        lines.append("\\begin{itemize}")
        lines.append("\item total number of events: "+nTXT)
        lines.append("\item integrated luminosity (lower limit): $%5.2F\mathrm{fb}^{-1}$"%(L/1E9))
        lines.append("\end{itemize}")
        lines.append("\end{frame}")
#
        R = list(self.runInfo.keys())
        R.sort(reverse=True)
        lines.append("\\begin{frame}{}")
        for i in range(len(R)//3):
           r = str(R[3*i]).zfill(6)
           lines.append("\includegraphics[width = 0.3\\textwidth]{Lumi-run"+r+".pdf}")
           if 3*i+1<len(R):
              r = str(R[3*i+1]).zfill(6)
              lines.append("\includegraphics[width = 0.3\\textwidth]{Lumi-run"+r+".pdf}")
           if 3*i+2<len(R):
              r = str(R[3*i+2]).zfill(6)
              lines.append("\includegraphics[width = 0.3\\textwidth]{Lumi-run"+r+".pdf}")
           if (3*i+3)%12==0: 
                lines.append("\end{frame}")
                lines.append("\\begin{frame}{}")
        lines.append("\end{frame}")
#
        lines.append(" ")
        lines.append("\\begin{frame}{}")
        for i in range(len(R)//3):
           r = str(R[3*i]).zfill(6)
           lines.append("\includegraphics[width = 0.3\\textwidth]{FS-run"+r+".pdf}")
           if 3*i+1<len(R):
              r = str(R[3*i+1]).zfill(6)
              lines.append("\includegraphics[width = 0.3\\textwidth]{FS-run"+r+".pdf}")
           if 3*i+2<len(R):
              r = str(R[3*i+2]).zfill(6)
              lines.append("\includegraphics[width = 0.3\\textwidth]{FS-run"+r+".pdf}")
           if (3*i+3)%12==0: 
                lines.append("\end{frame}")
                lines.append("\\begin{frame}{}")
        lines.append("\end{frame}")
        lines.append(" ")
#
        lines.append("\end{document}")
        outFile = open('LumiSummary.tex','w')
        for l in lines:
             rc = outFile.write(l+"\n")
        outFile.close()

   def plotLumiPerTime(self):
        h = self.h
        runInfo = self.runInfo
        h['LumiT']=ROOT.TGraph()
        h['ILumiT']=ROOT.TGraph()
        Lint = 0
        runList = list(runInfo.keys())
        runList.sort()
        lmax = 0
        for r in runList:
           h['LumiT'].AddPoint(runInfo[r]['StartTime'],0)
           X = runInfo[r]['lumiAtIP1withSNDLHC']/1E6
           h['LumiT'].AddPoint(runInfo[r]['StartTime']+1,X)
           if X>lmax : lmax = X
           h['LumiT'].AddPoint(runInfo[r]['StartTime']+3600,0)
           Lint+=X
           h['ILumiT'].AddPoint(runInfo[r]['StartTime'],Lint/1000)
        tstart = h['LumiT'].GetPointX(0)
        tend = h['LumiT'].GetPointX(h['LumiT'].GetN()-1)
        ut.bookHist(h,'LumiTime','; ; pb^{-1}',100,tstart,tend)
        h['LumiTime'].GetXaxis().SetTimeFormat("%d-%m")
        h['LumiTime'].GetXaxis().SetTimeOffset(0,'gmt')
        h['LumiTime'].SetMaximum(lmax*1.2)
        h['LumiTime'].SetStats(0)
        h['LumiTime'].Draw()
        h['c1'].Update()
# work with second axis
        rightmax = 1.1*Lint/1000.
        scale = ROOT.gPad.GetUymax()/rightmax
        h['ILumiT'].Scale(scale)
        h['ax2'] = ROOT.TGaxis(ROOT.gPad.GetUxmax(), ROOT.gPad.GetUymin(),
                   ROOT.gPad.GetUxmax(), ROOT.gPad.GetUymax(),
                   0, rightmax, 510, "+L")
        h['ax2'].SetTitle('L  Integral; fb^{-1}')
        h['ax2'].SetTextFont(42)
        h['ax2'].SetLabelFont(42)
        h['ax2'].SetTextColor(ROOT.kRed)
        h['ax2'].Draw()
        h['LumiT'].SetLineColor(ROOT.kBlue)
        h['LumiT'].SetLineWidth(2)
        h['LumiT'].Draw('Lsame')
        h['ILumiT'].SetLineColor(ROOT.kRed)
        h['ILumiT'].SetLineWidth(3)
        h['ILumiT'].Draw('same')
        h['c1'].Print('Lumi-time.root')
        h['c1'].Print('Lumi-time.pdf')
        h['c1'].Print('Lumi-run-time.png')

   def LumiIntegral(self,rmin,rmax):
        L = 0
        for r in self.runInfo:
           if r>= rmin and r<=rmax:
             L+=self.runInfo[r]['lumiAtIP1withSNDLHC']
        print('Lumi for run ',rmin,' to run ',rmax,' = ',L)
   def LumiPerFill(self):
      lumiPerFill={}
      for r in  self.runInfo:
         lumiPerFill[ self.runInfo[r]['Fillnumber'] ]=self.runInfo[r]['lumiAtIP1']
         I = 0
      for x in lumiPerFill:
         I+=lumiPerFill[x]
      print('total:',l)

if __name__ == '__main__':

    from argparse import ArgumentParser
    parser = ArgumentParser()

    parser.add_argument("-r", "--runNumbers", dest="runNumbers", help="list of run numbers",required=False, type=str,default="")
    parser.add_argument("-F", "--fillNumbers",  dest="fillNumbers",   help="corresponding fill number",type=str,required=False,default="")
    parser.add_argument("-c", "--command", dest="command",       help="command", default=None)
    parser.add_argument("-p", dest="path",       help="path to filling scheme", default="/mnt/hgfs/microDisk/SND@LHC/TI18/FillingSchemes/")
    parser.add_argument("-ip2", dest="withIP2", help="with IP2",default=True)

    options = parser.parse_args()
    FS = fillingScheme()
    FS.Init(options)
    ut.bookCanvas(FS.h,'c1','c1',1800,900,1,1)
    FS.h['c1'].cd()

    if options.command == "extract":
        FS.Extract()
    elif options.command == "draw":
        options.withIP2=True
        FS.plotBunchStructure(options.fillNumbers,int(options.runNumbers))
    elif options.command == "makeAll":
           offline =os.environ['EOSSHIP']+"/eos/experiment/sndlhc/www/offline.html"
           with client.File() as f:
               f.open(offline)
               status, L = f.read()
               Lhtml = L.decode().split('\n')
               f.close()
           runs = []
           for x in Lhtml:
                if x.find('#event')<0: continue
                ir = x.find("run ")
                if ir<0: continue
                k = x[ir:].find(" ")
                runnr = int(x[ir+4:ir+4+k+1])
                j = x.find("#events")
                tmp=x[j:].split('=')[1].split(' ')
                tmp.sort()
                for y in tmp:
                    if y=='': continue
                    nev = int( y )
                    break
                if nev < 100000: continue
                scale = 1
                k = x.find('post scaled:')
                if k>0: scale = int(x[k+12:])
                if runnr == 4425: continue   # no beam, but with fill number
                if runnr > 4360: runs.append([runnr,scale])
           for z in runs:
                 r = z[0]
                 options.postScale = z[1]
                 FS.options.fillNumbers = ""
                 FS.options.runNumbers = r
                 FS.Extract()
                 FS.drawLumi(r)
# fill dictionary with useful info
                 if not r in FS.FSdict: continue
                 FS.runInfo[r] = {'Fillnumber':FS.options.fillNumbers,'phaseShift1':FS.FSdict[r]['phaseShift1'],'phaseShift2':FS.FSdict[r]['phaseShift2'],
                                          'StartTime':FS.startTime,'StartTimeC':time.ctime(FS.startTime),'Entries':FS.getEntriesPerRun(r),
                                          'lumiAtIP1':FS.LumiInt[r][0],'lumiAtIP1withSNDLHC':FS.LumiInt[r][1],
                                          'OfflineMonitoring':"https://snd-lhc-monitoring.web.cern.ch/offline/run.html?file=run"+str(r).zfill(6)+".root&lastcycle"}
           FS.merge()
           import pickle
           p = open(options.path+"FSdict.pkl",'wb')
           pickle.dump(FS.FSdict,p)
           p.close()
           FS.mergeLumi()
           p = open(options.path+"Lumidict.pkl",'wb')
           pickle.dump(FS.LumiInt,p)
           p.close()
           p = open(options.path+"RunInfodict.pkl",'wb')
           pickle.dump(FS.runInfo,p)
           p.close()

           L,N= 0,0
           for r in FS.runInfo:
              L+=FS.runInfo[r]['lumiAtIP1withSNDLHC']
              N+=FS.runInfo[r]['Entries']
           print('total nr of events',N,'  integrated luminosity ',L/1E9,'fb')
           FS.makeLatex()
           # p = open(options.path+"FSdict.pkl",'rb')
           # FSdict = pickle.load(p)
