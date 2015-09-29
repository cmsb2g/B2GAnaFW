from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'TT_TuneCUETP8M1_13TeV'
config.General.workArea = 'crab_projects'

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../b2gedmntuples_cfg.py'
config.JobType.allowUndistributedCMSSW = True
config.JobType.pyCfgParams = ['DataProcessing=MC25ns','lheLabel=externalLHEProducer']
config.JobType.inputFiles = ['Summer15_25nsV2_DATA.db','Summer15_25nsV2_MC.db','Summer15_50nsV2_MC.db','Summer15_50nsV4_DATA.db','Summer15_50nsV4_MC.db','Summer15_50nsV5_DATA.db','Summer15_50nsV5_MC.db']

config.Data.inputDataset = '/TT_TuneCUETP8M1_13TeV-powheg-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2/MINIAODSIM'

config.Data.splitting = 'LumiBased'
#config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions15/13TeV/Cert_246908-256869_13TeV_PromptReco_Collisions15_25ns_JSON.txt'
#config.Data.runRange = '251585-251883'
config.Data.inputDBS = 'global'
config.Data.unitsPerJob = 30
config.Data.outLFNDirBase = '/store/user/decosa/ttDM/CMSSW_7_4_X' # or '/store/group/<subdir>'
config.Data.publication = True
config.Data.publishDBS = 'phys03'
config.Data.publishDataName = 'TT_TuneCUETP8M1_13TeV'

config.Site.storageSite = 'T2_CH_CSCS'






