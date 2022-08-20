import ROOT,os
import rootUtils as ut
from urllib.request import urlopen
import numpy

class fillingScheme():

   def Init(self,options):
     self.options = options
     self.h={}
     self.path = options.path
     self.content = ''
     self.phaseShift1 = 3564-404  # IP1
     self.phaseShift2 = 3564-404+764          # IP2

   def getFillNrFromRunNr(self,runNumber):
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
         if fillNumber<0: fillNumber = False
         print('fill number =',fillNumber )
       except:
         fillNumber = False
       return fillNumber

   def extractFillingScheme(self,fillNr):
       with urlopen('https://lpc.web.cern.ch/cgi-bin/schemeInfo.py?fill='+fillNr+'&fmt=json') as webpage:
           tmp = webpage.read().decode()
       content = ""
       exec("self.content = "+tmp)
       csv = self.content['fills'][fillNr]['csv'].split('\n')
       nB1 = csv.index('B1 bucket number,IP1,IP2,IP5,IP8')
       F = ROOT.TFile(self.path+'fillingScheme-'+fillNr+'.root','recreate')
       nt = ROOT.TNtuple('fill'+fillNr,'b1 IP1 IP2','B1:IP1:IP2')
       while nB1>0:
           tmp = csv[nB1+1].split(',')
           if len(tmp)!=5: break
           nB1+=1
           rc = nt.Fill(int(tmp[0]),int(tmp[1].replace('-','-1')),int(tmp[2].replace('-','-1')))
       nB2 = csv.index('B2 bucket number,IP1,IP2,IP5,IP8')
       while nB2>0:
           tmp = csv[nB2+1].split(',')
           if len(tmp)!=5: break
           nB2+=1
           rc = nt.Fill(-int(tmp[0]),int(tmp[1].replace('-','-1')),int(tmp[2].replace('-','-1')))
       nt.Write()
       F.Close()

   def extractPhaseShift(self,fillNr,runNumber):
         self.F = ROOT.TFile(self.path+'fillingScheme-'+fillNr+'.root')
         self.fs = self.F.Get('fill'+fillNr)
# convert to dictionary
         self.fsdict = {"B1":{},"B2":{}}
         for x in self.fs:
              if x.B1<0:
                   self.fsdict['B2'][(-x.B1-1)/10]={'IP1':x.IP1>0,'IP2':x.IP2>0}
              else:
                   self.fsdict['B1'][(x.B1-1)/10]={'IP1':x.IP1>0,'IP2':x.IP2>0}
         R = ROOT.TFile.Open(os.environ['EOSSHIP']+"/eos/experiment/sndlhc/www/offline/run"+str(runNumber).zfill(6)+".root")
         ROOT.gROOT.cd()
         self.h['bnr'] = R.daq.Get('bunchNumber').FindObject('bnr').Clone('bnr')
         R.Close()
         self.matches = {}
         for phase1 in range(0,3564):
               self.matches[phase1]=0
               for n in range(0,3564):
                   if not n in self.fsdict['B1']: continue
                   j = (n+phase1)%3564 + 1
                   if self.fsdict['B1'][n]['IP1']: self.matches[phase1]+=self.h['bnr'].GetBinContent(j)
         self.phaseShift1 = max(self.matches,key=self.matches.get)
         print('phaseShift1 found:',self.phaseShift1,3564-self.phaseShift1)
         self.matches = {}
         for phase2 in range(0,3564):
               self.matches[phase2]=0
               for n in range(0,3564):
                   if not n in self.fsdict['B2']: continue
                   j = (n+phase2)%3564 + 1
                   ip1 = (j-1+3564-self.phaseShift1)%3564
                   # take only bins which are not associated to collisions in IP1
                   if ip1 in self.fsdict['B1']:
                       if self.fsdict['B1'][ip1]['IP1']: continue
                   if self.fsdict['B2'][n]['IP2'] or 1>0: 
                      self.matches[phase2]+=self.h['bnr'].GetBinContent(j)
         self.phaseShift2 = max(self.matches,key=self.matches.get)
         print('phaseShift2 found:',self.phaseShift2,3564-self.phaseShift2)

   def plotBunchStructure(self,fillNr,runNumber,withIP2=False):
         h=self.h
         self.F = ROOT.TFile(self.path+'fillingScheme-'+fillNr+'.root')
         self.fs = self.F.Get('fill'+fillNr)
         ut.bookHist(h,'b1','b1',35640,-0.5,35639.5)
         ut.bookHist(h,'IP1','IP1',35640,-0.5,35639.5)
         ut.bookHist(h,'IP2','IP2',35640,-0.5,35639.5)
         ut.bookHist(h,'b1z','b1',3564,-0.5,3563.5)
         ut.bookHist(h,'IP1z','IP1',3564,-0.5,3563.5)
         ut.bookHist(h,'IP2z','IP2',3564,-0.5,3563.5)
         h['b1'].Draw()
         h['b1z'].SetLineColor(ROOT.kBlue)
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

         self.Draw(withIP2)

   def Draw(self,withIP2=False):
         h = self.h
         h['c1'].cd()
         self.fs.Draw('(( (B1-1)/10+'+str(self.phaseShift1)+')%3564)>>b1z','B1>-0.5','hist')
         self.fs.Draw('(( (B1-1)/10+'+str(self.phaseShift1)+')%3564)>>IP1z','IP1>-0.6&&B1>-0.5','hist')
         self.fs.Draw('(( (IP2-1)/10+'+str(self.phaseShift2)+')%3564)>>IP2z','IP2>-0.6&&B1<0.','hist')
         norm = h['bnr'].GetBinContent(h['bnr'].GetMaximumBin())
         h['b1z'].Scale(norm*1.5)
         h['IP1z'].Scale(norm*1.0)
         if withIP2: h['IP2z'].Scale(norm*0.3)
         h['bnr'].SetStats(0)
         h['b1z'].Draw('hist')
         if withIP2: h['IP2z'].Draw('histsame')
         h['IP1z'].Draw('histsame')
         h['bnr'].Draw('histsame')

if __name__ == '__main__':

    from argparse import ArgumentParser
    parser = ArgumentParser()

    parser.add_argument("-r", "--runNumbers", dest="runNumbers", help="list of run numbers",required=False, type=str,default="")
    parser.add_argument("-F", "--fillNumbers",  dest="fillNumbers",   help="corresponding fill number",type=str,required=False,default="")
    parser.add_argument("-c", "--command", dest="command",       help="command", default=None)
    parser.add_argument("-p", dest="path",       help="path to filling scheme", default="/mnt/hgfs/microDisk/SND@LHC/TI18/FillingSchemes/")

    options = parser.parse_args()
    FS = fillingScheme()
    FS.Init(options)

    if options.command == "extract":
        ut.bookCanvas(FS.h,'c1','c1',1800,900,1,1)
        FS.h['c1'].cd()
        if options.fillNumbers=='':
           fillNumber = FS.getFillNrFromRunNr(int(options.runNumbers))
           if not fillNumber: 
               print('Fill number not found')
           else:
               FS.extractFillingScheme(str(fillNumber))
               options.fillNumbers = str(fillNumber)
               FS.extractPhaseShift(options.fillNumbers,options.runNumbers)
               FS.plotBunchStructure(options.fillNumbers,int(options.runNumbers),True)
        else:
           for r in options.fillNumbers.split(','):
              FS.extractFillingScheme(r)
    elif options.command == "draw":
        FS.plotBunchStructure(options.fillNumbers,int(options.runNumbers),withIP2=True)

