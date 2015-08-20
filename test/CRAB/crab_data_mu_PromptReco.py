from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'SingleMuon_Run2015B-PromptReco-v74x_V5_1'
config.General.workArea = 'crab_projects'

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../b2gedmntuples_cfg.py'
config.JobType.allowUndistributedCMSSW = True
config.JobType.pyCfgParams = ['globalTag=74X_dataRun2_Prompt_v0','isData=True', 'LHE=0', 'DataProcessing=PromptReco50ns']
config.JobType.inputFiles = ['Cert_246908-251883_13TeV_PromptReco_Collisions15_JSON_v2.txt', 'Summer15_50nsV4_DATA.db']

config.Data.inputDataset = '/SingleMuon/Run2015B-PromptReco-v1/MINIAOD'

config.Data.splitting = 'LumiBased'
config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions15/13TeV/Cert_246908-251883_13TeV_PromptReco_Collisions15_JSON_v2.txt'
config.Data.runRange = '251585-251883'
config.Data.inputDBS = 'global'
config.Data.unitsPerJob = 50
config.Data.outLFNDirBase = '/store/user/decosa/ttDM/Synchro' # or '/store/group/<subdir>'
config.Data.publication = True
config.Data.publishDBS = 'phys03'
config.Data.publishDataName = 'SingleMuon_Run2015B-PromptReco'

config.Site.storageSite = 'T2_CH_CSCS'






