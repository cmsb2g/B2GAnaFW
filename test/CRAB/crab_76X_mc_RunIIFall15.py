import os
b2gana_test_dir = os.environ['CMSSW_BASE']+'/src/Analysis/B2GAnaFW/test/'

# Input .db and .txt files for JEC/JER
import glob
input_files = glob.glob(b2gana_test_dir+'Fall15_25nsV2_MC*')

# B2G Ana Fwk tag, change it to the latest one
# Please always include it in the outputDatasetTag together with the run era
tag = "B2GAnaFW_76X_V2p0" 

from WMCore.Configuration import Configuration
config = Configuration()

config.section_('General')
config.General.transferOutputs = True
config.General.transferLogs = True
config.General.requestName = 'TTJets_amcatnloFXFX'
config.General.workArea = tag

config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = b2gana_test_dir+'b2gedmntuples_cfg.py'
config.JobType.inputFiles = input_files
config.JobType.pyCfgParams = ['DataProcessing=MC25ns_MiniAOD_76X']

config.section_('Data')
config.Data.inputDBS = 'global'
config.Data.inputDataset = '/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
config.Data.publication = True
config.Data.publishDBS = 'phys03'
config.Data.outputDatasetTag = tag+'_RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.outLFNDirBase = '/store/user/jkarancs/SusyAnalysis/B2GEdmNtuple'

config.section_('Site')
config.Site.storageSite = 'T2_HU_Budapest'
