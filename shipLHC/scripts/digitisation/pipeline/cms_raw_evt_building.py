import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import glob
import os
from pathlib import Path
  
options = VarParsing.VarParsing()  
  
options.register('rawDirectory',  
    '',  
    VarParsing.VarParsing.multiplicity.singleton,  
    VarParsing.VarParsing.varType.string,  
    "Directory containing raw directories"  
)  
  
options.register('convertedDirectory',  
    '',  
    VarParsing.VarParsing.multiplicity.singleton,  
    VarParsing.VarParsing.varType.string,  
    "Directory containing converted directories"  
)  
  
options.register('runNumber',  
    0,  
    VarParsing.VarParsing.multiplicity.singleton,  
    VarParsing.VarParsing.varType.int,  
    "Run number"  
)  
  
options.parseArguments()  
  
if not options.rawDirectory or not options.convertedDirectory or not options.runNumber:  
    raise ValueError("rawDirectory, convertedDirectory, and runNumber must all be provided")


process = cms.Process("Convert")
process.load("DQM.SiStripCommon.MessageLogger_cfi")

# process.MessageLogger = cms.Service("MessageLogger",  
#     cerr = cms.untracked.PSet(  
#         threshold = cms.untracked.string('WARNING'),  
#     ),  
# )

process.EvFDaqDirector = cms.Service("EvFDaqDirector",
    baseDir   = cms.untracked.string(options.rawDirectory),
    buBaseDir = cms.untracked.string(options.rawDirectory),
    runNumber = cms.untracked.uint32(options.runNumber),
)

process.FastMonitoringService = cms.Service("FastMonitoringService",
    filePerFwkStream = cms.untracked.bool(False),
    fastMonIntervals = cms.untracked.uint32(2),
)

raw_files = sorted(glob.glob(str(Path(options.rawDirectory) / f"run{options.runNumber:06d}" / "*index*.raw")))
raw_files = ["file:" + f for f in raw_files]
print("Converting", len(raw_files), ".raw files")

process.source = cms.Source("FedRawDataInputSource",
    fileNames        = cms.untracked.vstring(raw_files),
    verifyAdler32    = cms.untracked.bool(False),
    useL1EventID     = cms.untracked.bool(True),
    keepRawFiles  = cms.untracked.bool(True),
)

maxEvents = cms.PSet(input = cms.untracked.int32(-1))

output_dir = Path(options.convertedDirectory) / f"run{options.runNumber:06d}"  
os.makedirs(output_dir, exist_ok=True)  

output_file = str(output_dir / f"run{options.runNumber:06d}_converted.root")
if os.path.exists(output_file):  
    os.remove(output_file)
    print("Replaced", output_file)  

process.out = cms.OutputModule("PoolOutputModule",
    fileName       = cms.untracked.string(output_file),
    outputCommands = cms.untracked.vstring("drop *", "keep FEDRawDataCollection_*_*_*"),
)

process.e = cms.EndPath(process.out)