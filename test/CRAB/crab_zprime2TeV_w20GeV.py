from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'zprime2TeV_w20GeV'
config.General.workArea = 'PHYS14'
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'b2gedmntuples_cfg.py'
config.JobType.pyCfgParams = ['LHE=0']

config.section_("Data")
config.Data.inputDataset = '/ZPrimeToTTJets_M2000GeV_W20GeV_Tune4C_13TeV-madgraph-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.ignoreLocality = True
config.Data.publication = True
# This string is used to construct the output dataset name
config.Data.publishDBS = 'phys03'
config.Data.publishDataName = 'B2GAnaFW_PHYS14'

config.section_("Site")
config.Site.storageSite = 'T3_US_FNALLPC'

