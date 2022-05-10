import os,subprocess,time,multiprocessing
import pwd
import ROOT
ncpus = multiprocessing.cpu_count()

""" input runlist
     for each run in runlist. look for number of partitions
     for each partition, start conversion job.
     check consistency of output file
     copy optional file to EOS
     run as many jobs in parallel as cpus are available
""" 

# use cases: H6, TI18
from argparse import ArgumentParser
parser = ArgumentParser()

parser.add_argument("-r", "--runNumbers", dest="runNumbers", help="list of run numbers",required=False, type=str,default="")
parser.add_argument("-P", "--production", dest="prod",       help="H6 / epfl / TI18"   ,required=False,default="TI18")
parser.add_argument("-c", "--command", dest="command",       help="command", default=None)
parser.add_argument("-o", "--overwrite", dest="overwrite",   action='store_true', help="overwrite EOS", default=False)
parser.add_argument("-cpp", "--convRawCPP", action='store_true', dest="FairTask_convRaw", help="convert raw data using ConvRawData FairTask", default=False)
parser.add_argument("--latest", dest="latest", help="last fully converted run", default=0,type=int)
parser.add_argument("-A", "--auto", dest="auto", help="run in auto mode checking regularly for new files",default=False,action='store_true')

global options
options = parser.parse_args()

Mext = ''  
if options.FairTask_convRaw: 
   Mext = '_CPP'  

def count_python_processes(macroName):
    username = pwd.getpwuid(os.getuid()).pw_name
    callstring = "ps -f -u " + username
# only works if screen is wide enough to print full name!
    status = subprocess.check_output(callstring,shell=True)
    n=0
    for x in status.decode().split('\n'):
        if not x.find(macroName)<0 and not x.find('python') <0: n+=1
    return n

def getPartitions(runList,path):
 partitions = {}
 for r in runList:
   directory = path+"run_"+str(r).zfill(6)
   dirlist  = str( subprocess.check_output("xrdfs "+os.environ['EOSSHIP']+" ls "+directory,shell=True) )
   partitions[r] = []
   for x in dirlist.split('\\n'):
        ix = x.find('data_')
        if ix<0: continue
        partitions[r].append( int(x[ix+5:ix+9]) )
 return partitions

def convert(runList,path,partitions={}):
 if len(partitions)==0: partitions = getPartitions(runList,path)
 for r in partitions:
    runNr   = str(r).zfill(6)
    # find partitions
    for p in partitions[r]:
       os.system("python $SNDSW_ROOT/shipLHC/rawData/runProd.py -c  'runSinglePartition;"+path+";"+runNr+";"+str(p).zfill(4)+";EOScopy=False;check=True;'   &")
       time.sleep(10)
       while count_python_processes('convertRawData')>ncpus:
          time.sleep(200)
 return partitions         

def runSinglePartition(path,r,p,EOScopy=False,check=True):
     inFile = os.environ['EOSSHIP']+path+'run_'+ r+'/data_'+p+'.root'
     fI = ROOT.TFile.Open(inFile)
     if not fI:
        print('file not found',path,r,p)
        exit()
     if not fI.Get('event'):
        print('file corrupted',path,r,p)
        exit()
     if options.FairTask_convRaw:
        os.system("python $SNDSW_ROOT/shipLHC/rawData/convertRawData.py -cpp -b 100000 -p "+path+"  -r "+str(int(r))+ " -P "+str(int(p)) + "  >log_"+r+'-'+p)
     else: 
        os.system("python $SNDSW_ROOT/shipLHC/rawData/convertRawData.py -b 1000 -p "+path+"  -r "+str(int(r))+ " -P "+str(int(p)) + " -g ../geofile_sndlhc_TI18.root >log_"+r+'-'+p)
     if check:
        rc = checkFile(path,r,p)
        if rc<0: 
            print('converting failed',r,p,rc)
            exit()
     print('finished converting ',r,p)
     tmp = {int(r):[int(p)]}
     if EOScopy:  copy2EOS(path,tmp)

def checkFile(path,r,p):
     inFile = os.environ['EOSSHIP']+path+'run_'+ r+'/data_'+p+'.root'
     fI = ROOT.TFile.Open(inFile)
     Nraw = fI.event.GetEntries()
     outFile = 'sndsw_raw_'+r+'-'+p+Mext+'.root'
     fC = ROOT.TFile(outFile)
     test = fC.Get('rawConv')
     if not test:
          print('Error:',path,r,p,' rawConv not found')
          return -2       
     Nconv = fC.rawConv.GetEntries()
     if Nraw != Nconv: 
          print('Error:',path,r,p,':',Nraw,Nconv)
          return -1
     return 0

def copy2EOS(path,success,overwrite=options.overwrite):
  for i in success:
    r = str(i).zfill(6)
    for x in success[i]:
      p = str(x).zfill(4)
      outFile = 'sndsw_raw_'+r+'-'+p+Mext+'.root'
      tmp = path.split('raw_data')[1].replace('data/','')
      pathConv = os.environ['EOSSHIP']+"/eos/experiment/sndlhc/convertedData/"+tmp+"run_"+r+"/sndsw_raw-"+p+".root"
      print('copy '+outFile+' to '+tmp+"run_"+r+"/sndsw_raw-"+p+".root")
      os.system('xrdcp '+outFile+'  '+pathConv)

def check(path,partitions):
 success = {}
 for r in partitions:
    success[r] = []
    runNr   = str(r).zfill(6)
    for x in partitions[r]:
       p = str(x).zfill(4)
       rc = checkFile(path,runNr,p)
       if rc==0: success[r].append(x)
 return success      
 
def getFileList(p,latest,minSize=10E6):
    inventory = {}
    dirList = str( subprocess.check_output("xrdfs "+os.environ['EOSSHIP']+" ls "+p,shell=True) )
    for x in dirList.split('\\n'):
          aDir = x[x.rfind('/')+1:]
          if not aDir.find('run')==0:continue
          runNr = int(aDir.split('_')[1])
          if not runNr > latest: continue
          fileList = str( subprocess.check_output("xrdfs "+os.environ['EOSSHIP']+" ls -l "+p+"/"+aDir,shell=True) )
          for z in fileList.split('\\n'):
               k = max(z.find('data_'),z.find('sndsw'))
               if not k>0: continue
               j = z.split(' /eos')[0].rfind(' ')
               size = int(z.split(' /eos')[0][j+1:])
               if size<minSize: continue
               tmp = z.split(' ')
               theDay = tmp[1] 
               theTime = tmp[2]
               fname = z[k:]
               run = int(aDir.split('_')[1])
               if run>900000: continue     # not a physical run
               k = fname.find('.root')
               partition = int(fname[k-4:k])
               d = theDay.split('-')
               t = theTime.split(':')
               gmt = time.mktime( (int(d[0]), int(d[1]), int(d[2]),  int(t[0]), int(t[1]), int(t[2]), 0, 0, 0 ) )
               inventory[run*10000+partition] = [aDir+"/"+fname,gmt]
    return inventory

def check4NewFiles(latest):
      rawDataFiles = getFileList(path,latest,minSize=10E6)
      convDataFiles = getFileList(pathConv,latest,minSize=10E6)
      orderedRDF = list(rawDataFiles.keys())
      orderedCDF = list(convDataFiles.keys())
      orderedRDF.reverse(),orderedCDF.reverse()
      for x in orderedRDF: 
           if not x > orderedCDF[0]: continue
           r = x//10000 
           p = x%10000
           runSinglePartition(path,r,p,EOScopy=True,check=True)

def getConvStats(runList):
  for run in runList:
     try: 
        f=ROOT.TFile.Open("sndsw_raw_"+str(run).zfill(6)+'.root')
        print(run,':',f.rawConv.GetEntries())
     except:
        print('problem:',run)
         
def rawStats(runList):
   for run in runList:
     runNr   = str(run).zfill(6)
     r = ROOT.TFile.Open(os.environ['EOSSHIP']+path+"/run_"+runNr+"/data.root")
     raw = r.event.GetEntries()
     print(run,':',raw)

def makeHistos(runList):
  # use always DS tracks to survey US because of larger overlap
  for run in runList:
    command = "$SNDSW_ROOT/shipLHC/scripts/Survey-MufiScifi.py -r "+str(run)+" -p "+pathConv+" -g geofile_sndlhc_H6.root -c Mufi_Efficiency -n -1 -t DS"
    os.system("python "+command+ " &")
    while count_python_processes('Survey-MufiScifi')>ncpus:
       time.sleep(200)

def mips():
  for run in runs:
    command = "Survey-MufiScifi.py -r "+str(run)+" -p "+pathConv+" -g geofile_sndlhc_H6.root -c mips"
    os.system("python "+command+ " &")
    while count_python_processes('Survey-MufiScifi')>multiprocessing.cpu_count()-2:
       time.sleep(200)

runList = []
if options.prod == "TI18":
      path = "/eos/experiment/sndlhc/raw_data/commissioning/TI18/data/"
      pathConv = "/eos/experiment/sndlhc/convertedData/commissioning/TI18/"
      if options.runNumbers=="": 
          runList = [1,6,7,8,16,18,19,20,21,23,24,25,26,27]
          # 6,7,8   14,15,22 corrupted
          # 
elif options.prod == "H6":
      path     = "/eos/experiment/sndlhc/raw_data/commissioning/TB_H6/data/"
      pathConv = "/eos/experiment/sndlhc/convertedData/commissioning/TB_H6/"
      if options.runNumbers=="": 
          runList = [273,276,283,284,285,286,295,296,290,294,298,301,298,292,23,27,34,35,39,40, 41, 42, 43, 44, 45, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 123, 146, 165, 184] 
elif options.prod == "epfl":
      path = "/eos/experiment/sndlhc/raw_data/commissioning/scifi_cosmics_epfl/data/"   
      if options.runNumbers=="": 
          runList = [3,8]

else:
      print("production not known. you are on your own",options.prod)


if options.auto:
    while 1 > 0:
         check4NewFiles(options.latest)
         time.sleep(1800)
    exit(0)

if options.command == "convert":

   for r in options.runNumbers.split(','):
      if r!= '': runList.append(int(r))
      
   convert(runList,path)    
   
elif options.command:
    tmp = options.command.split(';')
    command = tmp[0]+"("
    for i in range(1,len(tmp)):
         if i<4: command+='"'+tmp[i]+'"'
         else:    command+=tmp[i]
         if i<len(tmp)-1: command+=","
    command+=")"
    print('executing ' + command )
    eval(command)
