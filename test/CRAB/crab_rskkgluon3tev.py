from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'rskkgluon'
config.General.workArea = '741'
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'b2gedmntuples_cfg.py'
config.JobType.pyCfgParams = ['LHE=0']

config.section_("Data")
config.Data.inputDataset = '/RelValRSKKGluon_m3000GeV_13/CMSSW_7_4_1-MCRUN2_74_V9_gensim_740pre7-v1/MINIAODSIM'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.ignoreLocality = False
config.Data.publication = True
# This string is used to construct the output dataset name
config.Data.publishDBS = 'phys03'
config.Data.publishDataName = 'B2GAnaFW_741'

config.section_("Site")
config.Site.storageSite = 'T3_US_FNALLPC'

