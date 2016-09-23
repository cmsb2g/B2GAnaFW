import os
b2gana_test_dir = os.environ['CMSSW_BASE']+'/src/Analysis/B2GAnaFW/test/'

# Input .db and .txt files for JEC/JER
import glob
input_files = glob.glob(b2gana_test_dir+'Fall15_25nsV2_DATA*')

# B2G Ana Fwk tag, change it to the latest one
# Please always include it in the outputDatasetTag together with the run era
tag = "B2GAnaFW_76X_V2p0" 

from WMCore.Configuration import Configuration
config = Configuration()

config.section_('General')
config.General.transferOutputs = True
config.General.transferLogs = True
config.General.requestName = 'JetHT_2015D'
config.General.workArea = tag

config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = b2gana_test_dir+'b2gedmntuples_cfg.py'
config.JobType.inputFiles = input_files
config.JobType.pyCfgParams = ['DataProcessing=Data25ns_76X']

config.section_('Data')
config.Data.inputDBS = 'global'
config.Data.inputDataset = '/JetHT/Run2015D-16Dec2015-v1/MINIAOD'
config.Data.publication = True
config.Data.publishDBS = 'phys03'
config.Data.outputDatasetTag = tag+'_Run2015D-16Dec2015-v1'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions15/13TeV/Reprocessing/Cert_13TeV_16Dec2015ReReco_Collisions15_25ns_JSON_v2.txt'
config.Data.outLFNDirBase = '/store/user/jkarancs/SusyAnalysis/B2GEdmNtuple'

config.section_('Site')
config.Site.storageSite = 'T2_HU_Budapest'
