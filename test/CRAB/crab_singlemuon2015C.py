from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'singlemuon2015C'
config.General.workArea = 'B2GAnaFW_v74x_V7_25ns'
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'b2gedmntuples_cfg.py'
config.JobType.pyCfgParams = ['maxEvents=-1', 'DataProcessing=Data25nsv2']
config.JobType.inputFiles = ['Cert_246908-251252_13TeV_PromptReco_Collisions15_JSON.txt']

config.section_("Data")
config.Data.inputDataset = '/SingleMuon/Run2015D-PromptReco-v3/MINIAOD' 
config.Data.splitting = 'LumiBased'
config.Data.inputDBS = 'global'
config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions15/13TeV/Cert_246908-256869_13TeV_PromptReco_Collisions15_25ns_JSON.txt' 
config.Data.unitsPerJob = 20
config.Data.ignoreLocality = False
config.Data.publication = True
# This string is used to construct the output dataset name
config.Data.publishDBS = 'phys03'
config.Data.publishDataName = 'B2GAnaFW_v74x_V7_25ns'
config.Data.outLFNDirBase = '/store/user/devdatta/B2GAnaFW/SingleMuon/Run2015D-PromptReco-v3' 

config.section_("Site")
config.Site.storageSite = 'T3_US_FNALLPC'

