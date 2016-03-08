header = """
### *****************************************************************************************
### Usage:
###    The globalTag is automatically chosen according to the input 'DataProcessing' value. 
###    However it can be explictily specified to override the default option.
###    Remember that the value of 'DataProcessing' is not set by default. The user has the choice of 'Data25ns_76X' or 'MC25ns_MiniAOD_76X' or 'MC25ns_MiniAODv2_FastSim' 
###
### Examples: 
###    Running on 25 ns MiniAODv1 and MiniAODv2 MC in 76X:
###        cmsRun b2gedmntuples_cfg.py maxEvents=1000 DataProcessing='MC25ns_MiniAOD_76X'
###    Running on 25 ns 16Dec reprocessed data:
###        cmsRun b2gedmntuples_cfg.py maxEvents=1000 DataProcessing='Data25ns_76X'
###    Running on 25 ns FastSim MC:
###        cmsRun b2gedmntuples_cfg.py maxEvents=1000 DataProcessing='MC25ns_MiniAODv2_FastSim'
###
### *****************************************************************************************
"""
print header

import sys
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as opts
import copy

options = opts.VarParsing ('analysis')

options.register('sample',
     #'/store/mc/RunIIFall15MiniAODv2/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/10000/029802B3-83B8-E511-A002-0025905C22AE.root',
     '/store/mc/RunIIFall15MiniAODv2/RPVStopStopToJets_UDD312_M-120_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/50000/06C089EF-2FCC-E511-9946-1CC1DE1CF69A.root',
		 #'/store/data/Run2015D/JetHT/MINIAOD/16Dec2015-v1/00000/3085A2EF-6BB0-E511-87ED-0CC47A4D75EE.root',
     opts.VarParsing.multiplicity.singleton,
     opts.VarParsing.varType.string,
     'Sample to analyze')

options.register('outputLabel',
    'B2GEDMNtuple.root',
    opts.VarParsing.multiplicity.singleton,
    opts.VarParsing.varType.string,
    'Output label')

options.register('DataProcessing',
    "",
    opts.VarParsing.multiplicity.singleton,
    opts.VarParsing.varType.string,
    'Data processing types. Options are: Data25ns_76X or MC25ns_MiniAOD_76X or MC25ns_MiniAODv2_FastSim')

options.register('lheLabel',
    "",
    opts.VarParsing.multiplicity.singleton,
    opts.VarParsing.varType.string,
    'LHE module label, MC sample specific. Can be: externalLHEProducer')

### Expert options, do not change.
options.register('useNoHFMET',
    False,
    opts.VarParsing.multiplicity.singleton,
    opts.VarParsing.varType.bool,
    'Adding met without HF and relative jets')

options.register('usePrivateSQLite',
    True,
    opts.VarParsing.multiplicity.singleton,
    opts.VarParsing.varType.bool,
    'Take Corrections from private SQL file')

options.register('forceResiduals',
    True,
    opts.VarParsing.multiplicity.singleton,
    opts.VarParsing.varType.bool,
    'Whether to force residuals to be applied')

options.register('globalTag',
    '',
    opts.VarParsing.multiplicity.singleton,
    opts.VarParsing.varType.string,
    'Global Tag')

options.register('wantSummary',
    True,
    opts.VarParsing.multiplicity.singleton,
    opts.VarParsing.varType.bool,
    'Want summary report')

### Events to process: 'maxEvents' is already registered by the framework
options.setDefault('maxEvents', 100)

options.parseArguments()

if options.DataProcessing == "":
  sys.exit("!!!!ERROR: Enter 'DataProcessing' period. Options are: 'MC25ns_MiniAOD_76X', 'Data25ns_76X', 'MC25ns_MiniAODv2_FastSim'.\n")


if options.globalTag != "": 
  print "!!!!WARNING: You have chosen globalTag as", options.globalTag, ". Please check if this corresponds to your dataset."
else: 
  if options.DataProcessing=="MC25ns_MiniAOD_76X":
    options.globalTag="76X_mcRun2_asymptotic_v12"
  elif options.DataProcessing=="Data25ns_76X":
    options.globalTag="76X_dataRun2_v15"
  else:
    sys.exit("!!!!ERROR: Enter 'DataProcessing' period. Options are: 'MC25ns_MiniAOD_76X', 'Data25ns_76X', 'MC25ns_MiniAODv2_FastSim'.\n")

    if "Data" in options.DataProcessing:
      print "!!!!WARNING: You have chosen to run over data. lheLabel will be unset.\n"
  options.lheLabel = ""

###inputTag labels
rhoLabel          	= "fixedGridRhoFastjetAll"
muLabel           	= 'slimmedMuons'
elLabel           	= 'slimmedElectrons'
phoLabel            = 'slimmedPhotons'
pvLabel           	= 'offlineSlimmedPrimaryVertices'
convLabel         	= 'reducedEgamma:reducedConversions'
particleFlowLabel 	= 'packedPFCandidates'    
metLabel 		        = 'slimmedMETs'
metNoHFLabel 	     	= 'slimmedMETsNoHF'

triggerResultsLabel 	  = "TriggerResults"
triggerSummaryLabel 	  = "hltTriggerSummaryAOD"
hltMuonFilterLabel      = "hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f40QL3crIsoRhoFiltered0p15"
hltPathLabel            = "HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL"
hltElectronFilterLabel  = "hltL1sL1Mu3p5EG12ORL1MuOpenEG12L3Filtered8"

triggerResultsLabel    	= "TriggerResults"
triggerSummaryLabel    	= "hltTriggerSummaryAOD"
hltMuonFilterLabel     	= "hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f40QL3crIsoRhoFiltered0p15"
hltPathLabel           	= "HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL"
hltElectronFilterLabel 	= "hltL1sL1Mu3p5EG12ORL1MuOpenEG12L3Filtered8"

if(options.DataProcessing in [ "Data25ns_PromptRecov4","Data25ns_ReReco", "Data25ns_76X" ]): metProcess = "RECO"
else: metProcess = "PAT"

print "\nRunning with DataProcessing option ", options.DataProcessing, " and with global tag", options.globalTag, "\n" 


process = cms.Process("b2gEDMNtuples")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.MessageLogger.categories.append('HLTrigReport')
### Output Report
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(options.wantSummary) )
### Number of maximum events to process
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )
### Source file
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
      options.sample
      )
    )
### Setting global tag 
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag.globaltag = options.globalTag 


#process.load("PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff")
process.load("Configuration.EventContent.EventContent_cff")
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load("RecoEgamma/PhotonIdentification/PhotonIDValueMapProducer_cfi")
process.load("RecoEgamma.ElectronIdentification.ElectronIDValueMapProducer_cfi")
#process.load('RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV51_cff')


### External JECs =====================================================================================================
corrections = ['L1FastJet', 'L2Relative', 'L3Absolute']
if ("Data" in options.DataProcessing and options.forceResiduals): corrections.extend(['L2L3Residual'])

if options.usePrivateSQLite:
    
    from CondCore.DBCommon.CondDBSetup_cfi import *
    import os
    if "Data" in options.DataProcessing: era = "Fall15_25nsV2_DATA"
    elif "MC" in options.DataProcessing: era = "Fall15_25nsV2_MC"
    ###>>>elif "Data25ns" in options.DataProcessing:
    ###>>>  era = "Summer15_25nsV7_DATA"
    ###>>>elif "MC25ns" in options.DataProcessing:
    ###>>>  era = "Summer15_25nsV7_MC"
    else: sys.exit("!!!!ERROR: Enter 'DataProcessing' period. Options are: 'MC25ns_MiniAOD_76X', 'Data25ns_76X', 'MC25ns_MiniAODv2_FastSim'.\n")

    dBFile = era+".db"
    print "\nUsing private SQLite file", dBFile, "\n"
    process.jec = cms.ESSource("PoolDBESSource",CondDBSetup,
		    connect = cms.string( "sqlite_file:"+dBFile ),
		    toGet =  cms.VPSet(
			    cms.PSet(
				    record = cms.string("JetCorrectionsRecord"),
				    tag = cms.string("JetCorrectorParametersCollection_"+era+"_AK4PF"),
				    label= cms.untracked.string("AK4PF")
				    ),
			    cms.PSet(
				    record = cms.string("JetCorrectionsRecord"),
				    tag = cms.string("JetCorrectorParametersCollection_"+era+"_AK4PFchs"),
				    label= cms.untracked.string("AK4PFchs")
				    ),
			    cms.PSet(
				    record = cms.string("JetCorrectionsRecord"),
				    tag = cms.string("JetCorrectorParametersCollection_"+era+"_AK4PFPuppi"),
				    label= cms.untracked.string("AK4PFPuppi")
				    ),
			    cms.PSet(
				    record = cms.string("JetCorrectionsRecord"),
				    tag = cms.string("JetCorrectorParametersCollection_"+era+"_AK8PF"),
				    label= cms.untracked.string("AK8PF")
				    ),
			    cms.PSet(
				    record = cms.string("JetCorrectionsRecord"),
				    tag = cms.string("JetCorrectorParametersCollection_"+era+"_AK8PFchs"),
				    label= cms.untracked.string("AK8PFchs")
				    ),
			    cms.PSet(
				    record = cms.string("JetCorrectionsRecord"),
				    tag = cms.string("JetCorrectorParametersCollection_"+era+"_AK8PFPuppi"),
				    label= cms.untracked.string("AK8PFPuppi")
				    ),
			    )
		    )

    process.es_prefer_jec = cms.ESPrefer("PoolDBESSource",'jec')

    ###>>>''' 
    ###>>>process.load("PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff")
    ###>>>from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import patJetCorrFactorsUpdated, patJetsUpdated
    ###>>>process.patJetCorrFactorsReapplyJEC = patJetCorrFactorsUpdated.clone(
    ###>>>    rho = cms.InputTag("fixedGridRhoFastjetAll"),
    ###>>>    src = cms.InputTag("slimmedJets"),

    ###>>>    levels = corrections )
    ###>>>process.updatedPatJetsAK4 = patJetsUpdated.clone(
    ###>>>    jetSource = cms.InputTag("slimmedJets"),
    ###>>>    jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactorsReapplyJEC"))
    ###>>>    )
    ###>>>process.patJetAK8CorrFactorsReapplyJEC = patJetCorrFactorsUpdated.clone(
    ###>>>    src = cms.InputTag("slimmedJetsAK8"),
    ###>>>    rho = cms.InputTag("fixedGridRhoFastjetAll"),
    ###>>>    levels = corrections )

    ###>>>process.updatedPatJetsAK8 = patJetsUpdated.clone(
    ###>>>  jetSource = cms.InputTag("slimmedJetsAK8"),
    ###>>>  jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetAK8CorrFactorsReapplyJEC"))
    ###>>>  )
    ###>>>'''

    ### =====================================================================================================


### ------------------------------------------------------------------
### Recluster jets and adding subtructure tools from jetToolbox 
### (https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetToolbox)
### ------------------------------------------------------------------
from JMEAnalysis.JetToolbox.jetToolbox_cff import jetToolbox
listBtagDiscriminators = [ 
		'pfJetProbabilityBJetTags',
		'pfCombinedInclusiveSecondaryVertexV2BJetTags',
		'pfCombinedMVAV2BJetTags',
		'pfBoostedDoubleSecondaryVertexAK8BJetTags',
		'pfCombinedCvsLJetTags',
		'pfCombinedCvsBJetTags'
		]

ak4Cut='pt > 25 && abs(eta) < 5.'
ak8Cut='pt > 100 && abs(eta) < 5.'
if "MC" in options.DataProcessing: 
	jetToolbox( process, 'ak4', 'analysisPath', 'edmNtuplesOut', addQGTagger=True, bTagDiscriminators=listBtagDiscriminators, Cut=ak4Cut )
	jetToolbox( process, 'ak8', 'analysisPath', 'edmNtuplesOut', addSoftDropSubjets=True, addTrimming=True, rFiltTrim=0.1, addPruning=True, addFiltering=True, addSoftDrop=True, addNsub=True, bTagDiscriminators=listBtagDiscriminators, Cut=ak8Cut )
	jetToolbox( process, 'ca8', 'analysisPath', 'edmNtuplesOut', addCMSTopTagger=True, bTagDiscriminators=listBtagDiscriminators, Cut=ak8Cut )
	jetToolbox( process, 'ak8', 'analysisPath', 'edmNtuplesOut', PUMethod='Puppi', addSoftDropSubjets=True, addTrimming=True, addPruning=True, addFiltering=True, addSoftDrop=True, addNsub=True, bTagDiscriminators=listBtagDiscriminators, Cut=ak8Cut )
else:
	jetToolbox( process, 'ak4', 'analysisPath', 'edmNtuplesOut', runOnMC=False, addQGTagger=True, bTagDiscriminators=listBtagDiscriminators, Cut=ak4Cut )
	jetToolbox( process, 'ak8', 'analysisPath', 'edmNtuplesOut', runOnMC=False, addSoftDropSubjets=True, addTrimming=True, rFiltTrim=0.1, addPruning=True, addFiltering=True, addSoftDrop=True, addNsub=True, bTagDiscriminators=listBtagDiscriminators, Cut=ak8Cut )
	jetToolbox( process, 'ca8', 'analysisPath', 'edmNtuplesOut', runOnMC=False, addCMSTopTagger=True, bTagDiscriminators=listBtagDiscriminators, Cut=ak8Cut )
	jetToolbox( process, 'ak8', 'analysisPath', 'edmNtuplesOut', runOnMC=False, PUMethod='Puppi', addSoftDropSubjets=True, addTrimming=True, addPruning=True, addFiltering=True, addSoftDrop=True, addNsub=True, bTagDiscriminators=listBtagDiscriminators, Cut=ak8Cut )



jLabelAK8	= 'selectedPatJetsAK8PFCHS'
jLabelAK8Puppi 	= 'selectedPatJetsAK8PFPuppi'
jLabel		= 'selectedPatJetsAK4PFCHS'
jLabelNoHF	= 'selectedPatJetsAK4PFCHS'



### ---------------------------------------------------------------------------
### Removing the HF from the MET computation as from 7 Aug 2015 recommendations
### ---------------------------------------------------------------------------
if options.useNoHFMET:
  process.noHFCands = cms.EDFilter("CandPtrSelector",
      src=cms.InputTag("packedPFCandidates"),
      cut=cms.string("abs(pdgId)!=1 && abs(pdgId)!=2 && abs(eta)<3.0")
      )

  #jets are rebuilt from those candidates by the tools, no need to do anything else
### =================================================================================

from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD

#For a full met computation, remove the pfCandColl input
runMetCorAndUncFromMiniAOD(process,
		isData=("Data" in options.DataProcessing),
		)

if options.useNoHFMET:
	runMetCorAndUncFromMiniAOD(process,
		  isData=("Data" in options.DataProcessing),
		  pfCandColl=cms.InputTag("noHFCands"),
		  recoMetFromPFCs=True,
		  reclusterJets=True, 
		  postfix="NoHF"
		  )
	jLabelNoHF = 'patJetsNoHF'

### -------------------------------------------------------------------
### the lines below remove the L2L3 residual corrections when processing data
### -------------------------------------------------------------------

if ("Data" in options.DataProcessing and  options.forceResiduals):
  #Take new pat jets as input of the entuples
  #process.patJetCorrFactors.levels = corrections 
  if options.useNoHFMET:
    process.patJetCorrFactorsNoHF.levels = corrections 
    process.patPFMetT1T2CorrNoHF.jetCorrLabelRes = cms.InputTag("L3Absolute")
    process.patPFMetT1T2SmearCorrNoHF.jetCorrLabelRes = cms.InputTag("L3Absolute")
    process.patPFMetT2CorrNoHF.jetCorrLabelRes = cms.InputTag("L3Absolute")
    process.patPFMetT2SmearCorrNoHF.jetCorrLabelRes = cms.InputTag("L3Absolute")
    process.shiftedPatJetEnDownNoHF.jetCorrLabelUpToL3Res = cms.InputTag("ak4PFCHSL1FastL2L3Corrector")
    process.shiftedPatJetEnUpNoHF.jetCorrLabelUpToL3Res = cms.InputTag("ak4PFCHSL1FastL2L3Corrector")
  else: 
    process.patPFMetT1T2Corr.jetCorrLabelRes = cms.InputTag("L3Absolute")
    process.patPFMetT1T2SmearCorr.jetCorrLabelRes = cms.InputTag("L3Absolute")
    process.patPFMetT2Corr.jetCorrLabelRes = cms.InputTag("L3Absolute")
    process.patPFMetT2SmearCorr.jetCorrLabelRes = cms.InputTag("L3Absolute")
    process.shiftedPatJetEnDown.jetCorrLabelUpToL3Res = cms.InputTag("ak4PFCHSL1FastL2L3Corrector")
    process.shiftedPatJetEnUp.jetCorrLabelUpToL3Res = cms.InputTag("ak4PFCHSL1FastL2L3Corrector")
### ------------------------------------------------------------------


### ------------------------------------------------------------------
### Configure UserData
### ------------------------------------------------------------------

### Selected leptons and jets


process.skimmedPatMuons = cms.EDFilter(
    "PATMuonSelector",
    src = cms.InputTag(muLabel),
    cut = cms.string("pt > 0.0 && abs(eta) < 2.4")
    )

process.skimmedPatPhotons = cms.EDFilter(
    "PATPhotonSelector",
    src = cms.InputTag(phoLabel),
    cut = cms.string("pt > 10.0 && abs(eta) < 2.4"),
)

process.skimmedPatElectrons = cms.EDFilter(
    "PATElectronSelector",
    src = cms.InputTag(elLabel),
    cut = cms.string("pt > 10 && abs(eta) < 2.5")
    )

process.skimmedPatMET = cms.EDFilter(
    "PATMETSelector",
    #    src = cms.InputTag(metLabel, "", "PAT"),
    src = cms.InputTag(metLabel, "", metProcess),
    cut = cms.string("")
    )


##### THERE IS NO slimmedMETsNoHF in miniAODv2
'''
process.skimmedPatMETNoHF = cms.EDFilter(
    "PATMETSelector",
    src = cms.InputTag(metNoHFLabel, "", metProcess),
    cut = cms.string("")
    )
'''

process.eventUserData = cms.EDProducer(
    'EventUserData',
    pileup = cms.InputTag("slimmedAddPileupInfo"),
    pvSrc = cms.InputTag("offlineSlimmedPrimaryVertices")
    )

process.muonUserData = cms.EDProducer(
    'MuonUserData',
    muonLabel = cms.InputTag("skimmedPatMuons"),
    pv        = cms.InputTag(pvLabel),
    packedPFCands = cms.InputTag("packedPFCandidates"),
    ### TTRIGGER ###
    triggerResults = cms.InputTag(triggerResultsLabel,"","HLT"),
    triggerSummary = cms.InputTag(triggerSummaryLabel,"","HLT"),
    hltMuonFilter  = cms.InputTag(hltMuonFilterLabel),
    hltPath            = cms.string("HLT_IsoMu40_eta2p1_v11"),
    hlt2reco_deltaRmax = cms.double(0.1),
    # mainROOTFILEdir    = cms.string("../data/")
    )

process.jetUserData = cms.EDProducer(
    'JetUserData',
    jetLabel  = cms.InputTag(jLabel),
    jecCorrection      = cms.string("AK4PFchs"),
    ### TTRIGGER ###
    triggerResults = cms.InputTag(triggerResultsLabel,"","HLT"),
    triggerSummary = cms.InputTag(triggerSummaryLabel,"","HLT"),
    hltJetFilter       = cms.InputTag("hltPFHT"),
    hltPath            = cms.string("HLT_PFHT800"),
    hlt2reco_deltaRmax = cms.double(0.2),
    )


'''
process.jetUserDataNoHF = cms.EDProducer(
    'JetUserData',
    jetLabel  = cms.InputTag(jLabelNoHF),
    ### TTRIGGER ###
    triggerResults = cms.InputTag(triggerResultsLabel,"","HLT"),
    triggerSummary = cms.InputTag(triggerSummaryLabel,"","HLT"),
    hltJetFilter       = cms.InputTag("hltSixCenJet20L1FastJet"),
    hltPath            = cms.string("HLT_QuadJet60_DiJet20_v6"),
    hlt2reco_deltaRmax = cms.double(0.2),
    jecCorrection      = cms.string("AK4PFchs"),
    )
'''

process.jetUserDataAK8 = cms.EDProducer(
    'JetUserData',
    jetLabel  = cms.InputTag(jLabelAK8),
    pv        = cms.InputTag(pvLabel),
    jecCorrection      = cms.string("AK8PFchs"),
    ### TTRIGGER ###
    triggerResults = cms.InputTag(triggerResultsLabel,"","HLT"),
    triggerSummary = cms.InputTag(triggerSummaryLabel,"","HLT"),
    hltJetFilter       = cms.InputTag("hltAK8PFJetsTrimR0p1PT0p03"),
    hltPath            = cms.string("HLT_AK8PFHT650_TrimR0p1PT0p03Mass50"),
    hlt2reco_deltaRmax = cms.double(0.2)
)

process.boostedJetUserDataAK8 = cms.EDProducer(
    'BoostedJetToolboxUserData',
    jetLabel  = cms.InputTag('jetUserDataAK8'),
    topjetLabel = cms.InputTag('patJetsCMSTopTagCHSPacked'),
    vjetLabel = cms.InputTag('selectedPatJetsAK8PFCHSSoftDropPacked'),
    distMax = cms.double(0.8)
)
process.jetUserDataAK8Puppi = cms.EDProducer(
    'JetUserData',
    jetLabel  = cms.InputTag( jLabelAK8Puppi ),
    pv        = cms.InputTag(pvLabel),
    jecCorrection = cms.string("AK8PFPuppi"),
    ### TTRIGGER ###
    triggerResults = cms.InputTag(triggerResultsLabel,"","HLT"),
    triggerSummary = cms.InputTag(triggerSummaryLabel,"","HLT"),
    hltJetFilter       = cms.InputTag("hltAK8PFJetsTrimR0p1PT0p03"),
    hltPath            = cms.string("HLT_AK8PFHT650_TrimR0p1PT0p03Mass50"),
    hlt2reco_deltaRmax = cms.double(0.2)
    )

process.boostedJetUserDataAK8Puppi = cms.EDProducer(
    'BoostedJetToolboxUserData',
    jetLabel  = cms.InputTag('jetUserDataAK8Puppi'),
    topjetLabel = cms.InputTag('selectedPatJetsAK8PFPuppiSoftDropPacked'),
    vjetLabel = cms.InputTag('selectedPatJetsAK8PFPuppiSoftDropPacked'),
    distMax = cms.double(0.8)
    )

process.electronUserData = cms.EDProducer(
    'ElectronUserData',
    eleLabel   = cms.InputTag("skimmedPatElectrons"),
    pv         = cms.InputTag(pvLabel),
    packedPFCands = cms.InputTag("packedPFCandidates"),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    conversion = cms.InputTag(convLabel),
    rho        = cms.InputTag(rhoLabel),
    triggerResults = cms.InputTag(triggerResultsLabel),
    triggerSummary = cms.InputTag(triggerSummaryLabel),
    hltElectronFilter  = cms.InputTag(hltElectronFilterLabel),  ##trigger matching code to be fixed!
    hltPath             = cms.string("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL"),
    electronVetoIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-veto"),
    electronLooseIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-loose"),
    electronMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-medium"),
    electronTightIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-tight"),
    electronHEEPIdMap = cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV60"), 
    eleMediumIdFullInfoMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-medium"),
    eleIdVerbose = cms.bool(False)
    )

process.photonUserData = cms.EDProducer(
    'PhotonUserData',
    rho                     = cms.InputTag(rhoLabel),
    pholabel                = cms.InputTag("slimmedPhotons"),
    phoLooseIdMap           = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring15-25ns-V1-standalone-loose"),
    phoMediumIdMap          = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring15-25ns-V1-standalone-medium"),
    phoTightIdMap           = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring15-25ns-V1-standalone-tight"),
    phoChgIsoMap            = cms.InputTag("photonIDValueMapProducer:phoChargedIsolation"),
    phoPhoIsoMap            = cms.InputTag("photonIDValueMapProducer:phoPhotonIsolation"),
    phoNeuIsoMap            = cms.InputTag("photonIDValueMapProducer:phoNeutralHadronIsolation"),
    full5x5SigmaIEtaIEtaMap = cms.InputTag("photonIDValueMapProducer:phoFull5x5SigmaIEtaIEta")
    )

process.photonJets = cms.EDProducer(
    'PhotonJets',
    phoLabel = cms.InputTag("skimmedPatPhotons"),
    pv        = cms.InputTag(pvLabel),
    rho               = cms.InputTag(rhoLabel),
    packedPFCands = cms.InputTag("packedPFCandidates"),
    jetLabel  = cms.InputTag("slimmedJetsAK8"),
    ebReducedRecHitCollection = cms.InputTag("reducedEgamma:reducedEBRecHits"),
    eeReducedRecHitCollection = cms.InputTag("reducedEgamma:reducedEERecHits")

    )

process.vertexInfo = cms.EDProducer(
    'VertexInfo',
    src  = cms.InputTag(pvLabel),
    )

#
# Set up photon ID (VID framework)
#

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
# turn on VID producer, indicate data format  to be
# DataFormat.AOD or DataFormat.MiniAOD, as appropriate

#useAOD = False
#if useAOD == True :
dataFormat = DataFormat.MiniAOD
#else :
#    dataFormat = DataFormat.MiniAOD

switchOnVIDPhotonIdProducer(process, dataFormat)
switchOnVIDElectronIdProducer(process, dataFormat)

# define which IDs we want to produce
my_phoid_modules = ['RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Spring15_25ns_V1_cff']

#add them to the VID producer
for idmod in my_phoid_modules:
  setupAllVIDIdsInModule(process,idmod,setupVIDPhotonSelection)

#

my_eid_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_25ns_V1_cff',
    'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV60_cff']
for idmod in my_eid_modules:
  setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)


process.egmGsfElectronIDs.physicsObjectSrc = 'skimmedPatElectrons'

from PhysicsTools.PatAlgos.tools.pfTools import *
## Adapt primary vertex collection
adaptPVs(process, pvCollection=cms.InputTag('offlineSlimmedPrimaryVertices'))

#################################################


from PhysicsTools.CandAlgos.EventShapeVars_cff import *
process.eventShapePFVars = pfEventShapeVars.clone()
process.eventShapePFVars.src = cms.InputTag(particleFlowLabel)

process.eventShapePFJetVars = pfEventShapeVars.clone()
process.eventShapePFJetVars.src = cms.InputTag( jLabel )

process.centrality = cms.EDProducer("CentralityUserData",
    src = cms.InputTag( jLabel )
    )                                    

process.TriggerUserData = cms.EDProducer(
    'TriggerUserData',
    bits = cms.InputTag("TriggerResults","","HLT"),
    prescales = cms.InputTag("patTrigger"),
    storePrescales = cms.untracked.bool(True), 
    hltProcName = cms.untracked.string("HLT"), 
    objects = cms.InputTag("selectedPatTrigger")
    )                                 

process.METUserData = cms.EDProducer(
  'TriggerUserData',
  bits = cms.InputTag("TriggerResults","", metProcess),
  prescales = cms.InputTag("patTrigger"),
  storePrescales = cms.untracked.bool(False), 
  hltProcName = cms.untracked.string(metProcess), 
  objects = cms.InputTag("selectedPatTrigger")
  )

process.load('CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi')
process.HBHENoiseFilterResultProducer.minZeros = cms.int32(99999)
process.HBHENoiseFilterResultProducer.IgnoreTS4TS5ifJetInLowBVRegion=cms.bool(False)
process.HBHENoiseFilterResultProducer.defaultDecision = cms.string("HBHENoiseFilterResultRun2Loose")

### Including ntuplizer 
process.load("Analysis.B2GAnaFW.b2gedmntuples_cff")
process.options.allowUnscheduled = cms.untracked.bool(True)
from Analysis.B2GAnaFW.b2gedmntuples_cff import basic, jetVars, jetToolboxAK8Vars, jetToolboxAK8PuppiVars, jetVarsForSys
process.subjetKeysAK8.jetLabel = cms.InputTag("selectedPatJetsAK8PFCHSSoftDropPacked", "SubJets")
process.subjetsAK8.src = cms.InputTag("selectedPatJetsAK8PFCHSSoftDropPacked", "SubJets")
process.subjetsCmsTopTag.src = cms.InputTag("patJetsCMSTopTagCHSPacked", "SubJets")
process.subjetsCmsTopTagKeys.jetLabel = cms.InputTag("patJetsCMSTopTagCHSPacked", "SubJets")
process.jetsAK8.src = 'boostedJetUserDataAK8'
process.jetsAK8.variables += jetToolboxAK8Vars

### Puppi
process.jetsAK8Puppi = copy.deepcopy(basic)
process.jetsAK8Puppi.variables += jetVars
process.jetsAK8Puppi.variables += jetVarsForSys
process.jetsAK8Puppi.variables += jetToolboxAK8PuppiVars 
process.jetsAK8Puppi.prefix = 'jetAK8Puppi'
process.jetsAK8Puppi.src = cms.InputTag( 'boostedJetUserDataAK8Puppi' )
process.subjetsAK8Puppi = process.subjetsAK8.clone( prefix = 'subjetAK8Puppi', src = cms.InputTag('selectedPatJetsAK8PFPuppiSoftDropPacked', "SubJets") )
process.jetKeysAK8Puppi = process.jetKeysAK8.clone( jetLabel = 'boostedJetUserDataAK8Puppi' )
process.subjetKeysAK8Puppi = process.subjetKeysAK8.clone( jetLabel = cms.InputTag('selectedPatJetsAK8PFPuppiSoftDropPacked', "SubJets") )

process.edmNtuplesOut = cms.OutputModule(
    "PoolOutputModule",
    fileName = cms.untracked.string('B2GEDMNtuple.root'),
    outputCommands = cms.untracked.vstring(
    "drop *",
    "keep *_muons_*_*",
    "keep *_vertexInfo_*_*",
    "keep *_electrons_*_*",
    "keep *_photons_*_*",
    "keep *_photonjets_*_*",
    "keep *_jetsAK4_*_*",
    "keep *_jetsAK8*_*_*",
    "keep *_eventShape*_*_*",
    "keep *_*_*centrality*_*",
    "keep *_metFull_*_*",
    "keep *_metNoHF_*_*",
    "keep *_METUserData*_trigger*_*",
    "keep *_eventInfo_*_*",
    "keep *_subjetsAK8*_*_*",
    "keep *_subjetsCmsTopTag*_*_*",
    "keep *_jetKeysAK4_*_*",
    "keep *_jetKeysAK8*_*_*",
    "keep *_subjetKeysAK8*_*_*",
    "keep *_subjetsCmsTopTagKeys_*_*",
    "keep *_electronKeys_*_*",   
    "keep *_photonKeys_*_*",   
    "keep *_muonKeys_*_*",
    "keep *_TriggerUserData*_trigger*_*",
    "keep *_fixedGridRhoFastjetAll_*_*",
    "keep *_eventUserData_*_*",
    "keep *_eeBadScFilter_*_*"
    ),
    dropMetaData = cms.untracked.string('ALL'),
    )

# Some collections are not available in the current FastSim
if not "FastSim" in options.DataProcessing:
  process.edmNtuplesOut.outputCommands+=(
      "keep *_HBHENoiseFilterResultProducer_*_*",
      )


  ### keep NoHF jets if needed:
if( options.useNoHFMET ):
  process.edmNtuplesOut.outputCommands+=('keep *_jetsAK4NoHF_*_*',)

### keep info from LHEProducts if they are stored in PatTuples
if(options.lheLabel != ""):
  process.LHEUserData = cms.EDProducer("LHEUserData",
      lheLabel = cms.InputTag(options.lheLabel)
      )
  #process.analysisPath+=process.LHEUserData
  process.edmNtuplesOut.outputCommands+=('keep *_*LHE*_*_*',)

if "MC" in options.DataProcessing: 
  process.edmNtuplesOut.outputCommands+=(
      'keep *_generator_*_*',
      "keep *_genPart_*_*",
      "keep *_genJetsAK8_*_*",
      "keep *_genJetsAK8SoftDrop_*_*",
      "keep LHEEventProduct_*_*_*",
      "keep LHERunInfoProduct_*_*_*"
      )

process.edmNtuplesOut.fileName=options.outputLabel

#process.edmNtuplesOut.SelectEvents = cms.untracked.PSet(
#    SelectEvents = cms.vstring('filterPath')
#    )

#process.analysisPath += process.photonIDValueMapProducer

#process.fullPath = cms.Schedule(
#     process.analysisPath
#    )

process.endPath = cms.EndPath(process.edmNtuplesOut)

#open('B2GEntupleFileDump.py','w').write(process.dumpPython())
