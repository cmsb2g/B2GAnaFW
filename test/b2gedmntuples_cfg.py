header = """
### *****************************************************************************************
### Usage:
###    The globalTag is automatically chosen according to the input 'DataProcessing' value. 
###    However it can be explictily specified to override the default option.
###    Remember that the value of 'DataProcessing' is not set by default. The user has the choice of
###        'Data_80X', 'MC_MiniAODv2_80X', MC_MiniAODv2_80X_reHLT or 'MC_MiniAODv2_80X_FastSim' 
###
### Examples: 
###    Running on 25 ns MiniAODv1 and MiniAODv2 MC in 76X:
###        cmsRun b2gedmntuples_cfg.py maxEvents=1000 DataProcessing='MC25ns_MiniAOD_76X'
###    Running on 25 ns 16Dec reprocessed data:
###        cmsRun b2gedmntuples_cfg.py maxEvents=1000 DataProcessing='Data25ns_76X'
###    Running on 25 ns MC in 80x:
###        cmsRun b2gedmntuples_cfg.py maxEvents=1000 DataProcessing='MC_MiniAODv2_80X'
###    Running on 25 ns MC in 80x:
###        cmsRun b2gedmntuples_cfg.py maxEvents=1000 DataProcessing='MC_MiniAODv2_80X_reHLT'
###    Running on 25 ns FastSim MC in 80x:
###        cmsRun b2gedmntuples_cfg.py maxEvents=1000 DataProcessing='MC_MiniAODv2_80X_FastSim'
###    Running on 25 ns data in 80x:
###        cmsRun b2gedmntuples_cfg.py maxEvents=1000 DataProcessing='Data_80X'
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
#                 'file:00181849-176A-E511-8B11-848F69FD4C94.root',
'/store/mc/RunIISpring16MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14_ext3-v1/00000/0064B539-803A-E611-BDEA-002590D0B060.root',
#'/store/mc/RunIISpring16MiniAODv2/TprimeBToTH_M-1800_LH_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUSpring16RAWAODSIM_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/70000/E639441D-AD25-E611-93E2-B499BAABD482.root',
    #'/store/mc/RunIISpring15MiniAODv2/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v3/60000/00181849-176A-E511-8B11-848F69FD4C94.root', 
     #'/store/mc/RunIIFall15MiniAODv2/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/10000/029802B3-83B8-E511-A002-0025905C22AE.root',
     #'/store/mc/RunIIFall15MiniAODv2/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/00000/002253C9-DFB8-E511-8B0A-001A648F1C42.root',
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
    '',
    opts.VarParsing.multiplicity.singleton,
    opts.VarParsing.varType.string,
    'Data processing types. Options are: Data_80X, MC_MiniAODv2_80X, MC_MiniAODv2_80X_reHLT or MC_MiniAODv2_80X_FastSim')

### Expert options, do not change.
options.register('useNoHFMET',
    False,
    opts.VarParsing.multiplicity.singleton,
    opts.VarParsing.varType.bool,
    'Adding met without HF and relative jets')

options.register('usePrivateSQLite',
    False,
    opts.VarParsing.multiplicity.singleton,
    opts.VarParsing.varType.bool,
    'Take Jet Energy Corrections from private SQL file')

options.register('usePrivateSQLiteForJER',
    False,
    opts.VarParsing.multiplicity.singleton,
    opts.VarParsing.varType.bool,
    'Take Jet Enery Resolution from private SQL file')

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
  sys.exit("!!!!ERROR: Enter 'DataProcessing' period. Options are: Data_80X, MC_MiniAODv2_80X, MC_MiniAODv2_80X_reHLT or MC_MiniAODv2_80X_FastSim.\n")


if options.globalTag != "": 
  print "!!!!WARNING: You have chosen globalTag as", options.globalTag, ". Please check if this corresponds to your dataset."
else: 
  if options.DataProcessing=="MC_MiniAODv2_80X":
    options.globalTag="80X_mcRun2_asymptotic_2016_miniAODv2_v1"
  elif options.DataProcessing=="MC_MiniAODv2_80X_reHLT":
    options.globalTag="80X_mcRun2_asymptotic_2016_miniAODv2_v1"    
  elif options.DataProcessing=="MC_MiniAODv2_80X_FastSim":
    options.globalTag="80X_mcRun2_asymptotic_2016_miniAODv2_v1"
    options.usePrivateSQLite = True
  elif options.DataProcessing=="Data_80X":
    options.globalTag="80X_dataRun2_Prompt_ICHEP16JEC_v0"
  else:
    sys.exit("!!!!ERROR: Enter 'DataProcessing' period. Options are: Data_80X, MC_MiniAODv2_80X, MC_MiniAODv2_80X_reHLT or MC_MiniAODv2_80X_FastSim.\n")

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

if(options.DataProcessing in [ "MC_MiniAODv2_80X", "MC_MiniAODv2_80X_reHLT", "MC_MiniAODv2_80X_FastSim" ]): metProcess = "PAT"
else: metProcess = "RECO"

if options.DataProcessing =="MC_MiniAODv2_80X_reHLT": hltProcess = "HLT2"
else: hltProcess = "HLT"

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


### External JEC =====================================================================================================
corrections = ['L1FastJet', 'L2Relative', 'L3Absolute']
if ("Data" in options.DataProcessing and options.forceResiduals): corrections.extend(['L2L3Residual'])

if options.usePrivateSQLite:
    if options.DataProcessing=="Data_80X":
      jec_era = "Spring16_25nsV6_DATA"
    elif options.DataProcessing=="MC_MiniAODv2_80X":
      jec_era = "Spring16_25nsV6_MC"
    elif options.DataProcessing=="MC_MiniAODv2_80X_FastSim":
      jec_era = "Spring16_25nsFastSimMC_V1"
    else: sys.exit("!!!!ERROR: Enter 'DataProcessing' period. Options are: Data_80X, MC_MiniAODv2_80X, MC_MiniAODv2_80X_reHLT or MC_MiniAODv2_80X_FastSim.\n")
    
    from CondCore.DBCommon.CondDBSetup_cfi import *
    import os
    
    # JEC
    dBFile = jec_era+".db"
    print "\nUsing private SQLite file", dBFile, "\n"
    process.jec = cms.ESSource("PoolDBESSource",CondDBSetup,
		    connect = cms.string( "sqlite_file:"+dBFile ),
		    toGet =  cms.VPSet(
			    cms.PSet(
				    record = cms.string("JetCorrectionsRecord"),
				    tag = cms.string("JetCorrectorParametersCollection_"+jec_era+"_AK4PF"),
				    label= cms.untracked.string("AK4PF")
				    ),
			    cms.PSet(
				    record = cms.string("JetCorrectionsRecord"),
				    tag = cms.string("JetCorrectorParametersCollection_"+jec_era+"_AK4PFchs"),
				    label= cms.untracked.string("AK4PFchs")
				    ),
			    cms.PSet(
				    record = cms.string("JetCorrectionsRecord"),
				    tag = cms.string("JetCorrectorParametersCollection_"+jec_era+"_AK8PF"),
				    label= cms.untracked.string("AK8PF")
				    ),
			    cms.PSet(
				    record = cms.string("JetCorrectionsRecord"),
				    tag = cms.string("JetCorrectorParametersCollection_"+jec_era+"_AK8PFchs"),
				    label= cms.untracked.string("AK8PFchs")
				    ),
			    )
		    )
    # FastSim JEC is not available for Puppi jets
    if "FastSim" not in options.DataProcessing:
      process.jec.toGet.append(cms.PSet(
        			    record = cms.string("JetCorrectionsRecord"),
        			    tag = cms.string("JetCorrectorParametersCollection_"+jec_era+"_AK4PFPuppi"),
        			    label= cms.untracked.string("AK4PFPuppi")
        			    )
                              )
      process.jec.toGet.append(cms.PSet(
        			    record = cms.string("JetCorrectionsRecord"),
        			    tag = cms.string("JetCorrectorParametersCollection_"+jec_era+"_AK8PFPuppi"),
        			    label= cms.untracked.string("AK8PFPuppi")
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

### External JER =====================================================================================================
if "Data" in options.DataProcessing:
  jer_era = "Spring16_25nsV6_DATA"
elif "MC" in options.DataProcessing:
  jer_era = "Spring16_25nsV6_MC"
else: sys.exit("!!!!ERROR: Enter 'DataProcessing' period. Options are: Data_80X, MC_MiniAODv2_80X, MC_MiniAODv2_80X_reHLT or MC_MiniAODv2_80X_FastSim.\n")

if options.usePrivateSQLiteForJER:
    # JER
    #   https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyResolution#Accessing_factors_from_Global_Ta
    process.load("JetMETCorrections.Modules.JetResolutionESProducer_cfi")
    process.jer = cms.ESSource("PoolDBESSource",CondDBSetup,
            connect = cms.string('sqlite_file:'+jer_era+"_JER.db"), # '_JER' added to filename to distinguish from JEC file
            toGet = cms.VPSet(
                # Resolution
                cms.PSet(
                    record = cms.string('JetResolutionRcd'),
                    tag    = cms.string('JR_'+jer_era+'_PtResolution_AK4PF'),
                    label  = cms.untracked.string('AK4PF_pt')
                    ),
                cms.PSet(
                    record = cms.string('JetResolutionRcd'),
                    tag    = cms.string('JR_'+jer_era+'_PtResolution_AK4PFchs'),
                    label  = cms.untracked.string('AK4PFchs_pt')
                    ),
                cms.PSet(
                    record = cms.string('JetResolutionRcd'),
                    tag    = cms.string('JR_'+jer_era+'_PtResolution_AK4PFPuppi'),
                    label  = cms.untracked.string('AK4PFPuppi_pt')
                    ),
                cms.PSet(
                    record = cms.string('JetResolutionRcd'),
                    tag    = cms.string('JR_'+jer_era+'_PtResolution_AK8PF'),
                    label  = cms.untracked.string('AK8PF_pt')
                    ),
                cms.PSet(
                    record = cms.string('JetResolutionRcd'),
                    tag    = cms.string('JR_'+jer_era+'_PtResolution_AK8PFchs'),
                    label  = cms.untracked.string('AK8PFchs_pt')
                    ),
                cms.PSet(
                    record = cms.string('JetResolutionRcd'),
                    tag    = cms.string('JR_'+jer_era+'_PtResolution_AK8PFPuppi'),
                    label  = cms.untracked.string('AK8PFPuppi_pt')
                    ),
                # Scale factors
                cms.PSet(
                    record = cms.string('JetResolutionScaleFactorRcd'),
                    tag    = cms.string('JR_'+jer_era+'_SF_AK4PF'),
                    label  = cms.untracked.string('AK4PF')
                    ),
                cms.PSet(
                    record = cms.string('JetResolutionScaleFactorRcd'),
                    tag    = cms.string('JR_'+jer_era+'_SF_AK4PFchs'),
                    label  = cms.untracked.string('AK4PFchs')
                    ),
                cms.PSet(
                    record = cms.string('JetResolutionScaleFactorRcd'),
                    tag    = cms.string('JR_'+jer_era+'_SF_AK4PFPuppi'),
                    label  = cms.untracked.string('AK4PFPuppi')
                    ),
                cms.PSet(
                    record = cms.string('JetResolutionScaleFactorRcd'),
                    tag    = cms.string('JR_'+jer_era+'_SF_AK8PF'),
                    label  = cms.untracked.string('AK8PF')
                    ),
                cms.PSet(
                    record = cms.string('JetResolutionScaleFactorRcd'),
                    tag    = cms.string('JR_'+jer_era+'_SF_AK8PFchs'),
                    label  = cms.untracked.string('AK8PFchs')
                    ),
                cms.PSet(
                    record = cms.string('JetResolutionScaleFactorRcd'),
                    tag    = cms.string('JR_'+jer_era+'_SF_AK8PFPuppi'),
                    label  = cms.untracked.string('AK8PFPuppi')
                    ),
                ),
            )
    process.es_prefer_jer = cms.ESPrefer('PoolDBESSource', 'jer')


### ------------------------------------------------------------------
### Puppi
### (https://twiki.cern.ch/twiki/bin/viewauth/CMS/PUPPI)
### ------------------------------------------------------------------
process.load('CommonTools/PileupAlgos/Puppi_cff')
process.puppi.candName = cms.InputTag('packedPFCandidates')
process.puppi.vertexName = cms.InputTag('offlineSlimmedPrimaryVertices')
process.puppi.useExistingWeights = cms.bool(True)

#process.puppiOnTheFly = process.puppi.clone()

#process.puppiOnTheFly.useExistingWeights = True

### ------------------------------------------------------------------
### Recluster jets and adding subtructure tools from jetToolbox 
### (https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetToolbox)
### ------------------------------------------------------------------
from JMEAnalysis.JetToolbox.jetToolbox_cff import jetToolbox
listBTagInfos = [
     'pfInclusiveSecondaryVertexFinderTagInfos',
     ]
listBtagDiscriminatorsAK4 = [ 
		'pfJetProbabilityBJetTags',
		'pfCombinedInclusiveSecondaryVertexV2BJetTags',
		'pfCombinedMVAV2BJetTags',
		'pfCombinedCvsLJetTags',
		'pfCombinedCvsBJetTags',
		]
listBtagDiscriminatorsAK8 = [ 
		'pfJetProbabilityBJetTags',
		'pfCombinedInclusiveSecondaryVertexV2BJetTags',
		'pfCombinedMVAV2BJetTags',
		'pfCombinedCvsLJetTags',
		'pfCombinedCvsBJetTags',
    'pfBoostedDoubleSecondaryVertexAK8BJetTags',
    'pfBoostedDoubleSecondaryVertexCA15BJetTags',
		]

runMC = ("MC" in options.DataProcessing)

ak4Cut='pt > 25 && abs(eta) < 5.'
ak8Cut='pt > 200 && abs(eta) < 2.4'
jetToolbox( process, 'ak4', 'analysisPath', 'edmNtuplesOut', runOnMC=runMC, addQGTagger=True, bTagInfos=listBTagInfos, bTagDiscriminators=listBtagDiscriminatorsAK4, Cut=ak4Cut )
jetToolbox( process, 'ak4', 'analysisPath', 'edmNtuplesOut', runOnMC=runMC, PUMethod='Puppi', newPFCollection=True, nameNewPFCollection='puppi', bTagInfos=listBTagInfos, bTagDiscriminators=listBtagDiscriminatorsAK4, Cut=ak4Cut )
jetToolbox( process, 'ak8', 'analysisPath', 'edmNtuplesOut', runOnMC=runMC, addSoftDropSubjets=True, addTrimming=True, rFiltTrim=0.1, addPruning=True, addFiltering=True, addSoftDrop=True, addNsub=True, bTagInfos=listBTagInfos, bTagDiscriminators=listBtagDiscriminatorsAK8, addCMSTopTagger=True, Cut=ak8Cut , addNsubSubjets=True, subjetMaxTau=4 )
jetToolbox( process, 'ak8', 'analysisPath', 'edmNtuplesOut', runOnMC=runMC, PUMethod='Puppi', newPFCollection=True, nameNewPFCollection='puppi', addSoftDropSubjets=True, addTrimming=True, addPruning=True, addFiltering=True, addSoftDrop=True, addNsub=True, bTagInfos=listBTagInfos, bTagDiscriminators=listBtagDiscriminatorsAK8, addCMSTopTagger=True, Cut=ak8Cut, addNsubSubjets=True, subjetMaxTau=4 )
# Added , addNsubSubjets=True, subjetMaxTau=4  to both ak8 above

jLabel		= 'selectedPatJetsAK4PFCHS'
jLabelAK8	= 'selectedPatJetsAK8PFCHS'
jLabelPuppi	= 'selectedPatJetsAK4PFPuppi'
jLabelAK8Puppi 	= 'selectedPatJetsAK8PFPuppi'

# JER Twiki:
#   https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyResolution#Scale_factors
# Get Latest txt files from:
#   https://github.com/cms-jet/JRDatabase/tree/master/textFiles
jetAlgo         = 'AK4PFchs'
jetAlgoPuppi    = 'AK4PFPuppi'
jetAlgoAK8      = 'AK8PFchs'
jetAlgoAK8Puppi = 'AK8PFPuppi'


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

### -------------------------------------------------------------------
### Latest Run II MET Filter recommendations
### -------------------------------------------------------------------

# Using recipe:
# https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2?rev=99#How_to_run_the_Bad_Charged_Hadro

process.load('RecoMET.METFilters.BadPFMuonFilter_cfi')
process.BadPFMuonFilter.muons = cms.InputTag("slimmedMuons")
process.BadPFMuonFilter.PFCandidates = cms.InputTag("packedPFCandidates")

process.load('RecoMET.METFilters.BadChargedCandidateFilter_cfi')
process.BadChargedCandidateFilter.muons = cms.InputTag("slimmedMuons")
process.BadChargedCandidateFilter.PFCandidates = cms.InputTag("packedPFCandidates")


### ------------------------------------------------------------------
### Configure UserData
### ------------------------------------------------------------------

### Selected leptons and jets

### Check PackedCandidate vertex to make sure there is one
process.chs.cut = 'vertexRef.isNonnull() && fromPV()'

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
    triggerResults = cms.InputTag(triggerResultsLabel,"",hltProcess),
    triggerSummary = cms.InputTag(triggerSummaryLabel,"",hltProcess),
    hltMuonFilter  = cms.InputTag(hltMuonFilterLabel),
    hltPath            = cms.string("HLT_IsoMu40_eta2p1_v11"),
    hlt2reco_deltaRmax = cms.double(0.1),
    # mainROOTFILEdir    = cms.string("../data/")
    )

process.jetUserData = cms.EDProducer(
    'JetUserData',
    jetLabel          = cms.InputTag(jLabel),
    rho               = cms.InputTag(rhoLabel),
    coneSize          = cms.double(0.4),
    getJERFromTxt     = cms.bool(False),
    jetCorrLabel      = cms.string(jetAlgo),
    jerLabel          = cms.string(jetAlgo),
    resolutionsFile   = cms.string(jer_era+'_PtResolution_'+jetAlgo+'.txt'),
    scaleFactorsFile  = cms.string(jer_era+'_SF_'+jetAlgo+'.txt'),
    ### TTRIGGER ###
    triggerResults = cms.InputTag(triggerResultsLabel,"",hltProcess),
    triggerSummary = cms.InputTag(triggerSummaryLabel,"",hltProcess),
    hltJetFilter       = cms.InputTag("hltPFHT"),
    hltPath            = cms.string("HLT_PFHT800"),
    hlt2reco_deltaRmax = cms.double(0.2),
    candSVTagInfos         = cms.string("pfInclusiveSecondaryVertexFinder"), 
    )

process.jetUserDataPuppi = cms.EDProducer(
    'JetUserData',
    jetLabel          = cms.InputTag(jLabelPuppi),
    rho               = cms.InputTag(rhoLabel),
    coneSize          = cms.double(0.4),
    getJERFromTxt     = cms.bool(False),
    jetCorrLabel      = cms.string(jetAlgoPuppi),
    jerLabel          = cms.string(jetAlgoPuppi),
    resolutionsFile   = cms.string(jer_era+'_PtResolution_'+jetAlgoPuppi+'.txt'),
    scaleFactorsFile  = cms.string(jer_era+'_SF_'+jetAlgoPuppi+'.txt'),
    ### TTRIGGER ###
    triggerResults = cms.InputTag(triggerResultsLabel,"",hltProcess),
    triggerSummary = cms.InputTag(triggerSummaryLabel,"",hltProcess),
    hltJetFilter       = cms.InputTag("hltPFHT"),
    hltPath            = cms.string("HLT_PFHT800"),
    hlt2reco_deltaRmax = cms.double(0.2),
    candSVTagInfos         = cms.string("pfInclusiveSecondaryVertexFinder"), 
    )


process.jetUserDataAK8 = cms.EDProducer(
    'JetUserData',
    jetLabel          = cms.InputTag(jLabelAK8),
    rho               = cms.InputTag(rhoLabel),
    coneSize          = cms.double(0.8),
    getJERFromTxt     = cms.bool(False),
    jetCorrLabel      = cms.string(jetAlgoAK8),
    jerLabel          = cms.string(jetAlgoAK8),
    resolutionsFile   = cms.string(jer_era+'_PtResolution_'+jetAlgoAK8+'.txt'),
    scaleFactorsFile  = cms.string(jer_era+'_SF_'+jetAlgoAK8+'.txt'),
    ### TTRIGGER ###
    triggerResults = cms.InputTag(triggerResultsLabel,"",hltProcess),
    triggerSummary = cms.InputTag(triggerSummaryLabel,"",hltProcess),
    hltJetFilter       = cms.InputTag("hltAK8PFJetsTrimR0p1PT0p03"),
    hltPath            = cms.string("HLT_AK8PFHT650_TrimR0p1PT0p03Mass50"),
    hlt2reco_deltaRmax = cms.double(0.2), 
    candSVTagInfos         = cms.string("pfInclusiveSecondaryVertexFinder"), 
)

process.boostedJetUserDataAK8 = cms.EDProducer(
    'BoostedJetToolboxUserData',
    jetLabel  = cms.InputTag('jetUserDataAK8'),
    #topjetLabel = cms.InputTag('patJetsCMSTopTagCHSPacked'),
    vjetLabel = cms.InputTag('selectedPatJetsAK8PFCHSSoftDropPacked'),
    distMax = cms.double(0.8)
)


process.jetUserDataAK8Puppi = cms.EDProducer(
    'JetUserData',
    jetLabel          = cms.InputTag( jLabelAK8Puppi ),
    rho               = cms.InputTag(rhoLabel),
    coneSize          = cms.double(0.8),
    getJERFromTxt     = cms.bool(False),
    jetCorrLabel      = cms.string(jetAlgoAK8Puppi),
    jerLabel          = cms.string(jetAlgoAK8Puppi),
    resolutionsFile   = cms.string(jer_era+'_PtResolution_'+jetAlgoAK8Puppi+'.txt'),
    scaleFactorsFile  = cms.string(jer_era+'_SF_'+jetAlgoAK8Puppi+'.txt'),
    ### TTRIGGER ###
    triggerResults = cms.InputTag(triggerResultsLabel,"",hltProcess),
    triggerSummary = cms.InputTag(triggerSummaryLabel,"",hltProcess),
    hltJetFilter       = cms.InputTag("hltAK8PFJetsTrimR0p1PT0p03"),
    hltPath            = cms.string("HLT_AK8PFHT650_TrimR0p1PT0p03Mass50"),
    hlt2reco_deltaRmax = cms.double(0.2), 
    candSVTagInfos         = cms.string("pfInclusiveSecondaryVertexFinder"), 
    )

# FastSim JEC is not available for Puppi jets, use CHS instead
if "FastSim" in options.DataProcessing:
  process.jetUserDataPuppi.jetCorrLabel    = jetAlgo
  process.jetUserDataAK8Puppi.jetCorrLabel = jetAlgoAK8

process.boostedJetUserDataAK8Puppi = cms.EDProducer(
    'BoostedJetToolboxUserData',
    jetLabel  = cms.InputTag('jetUserDataAK8Puppi'),
    #topjetLabel = cms.InputTag('patJetsCMSTopTagPuppiPacked'),
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
    triggerResults = cms.InputTag(triggerResultsLabel,"",hltProcess),
    triggerSummary = cms.InputTag(triggerSummaryLabel,"",hltProcess),
    hltElectronFilter  = cms.InputTag(hltElectronFilterLabel),  ##trigger matching code to be fixed!
    hltPath             = cms.string("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL"),
    eleVetoIdFullInfoMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-veto"),
    eleLooseIdFullInfoMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-loose"),
    eleMediumIdFullInfoMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-medium"),
    eleTightIdFullInfoMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-tight"),
    eleHEEPIdIdFullInfoMap = cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV60"),
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

my_eid_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Summer16_80X_V1_cff',
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
    bits = cms.InputTag(triggerResultsLabel,"",hltProcess),
    prescales = cms.InputTag("patTrigger"),
    storePrescales = cms.untracked.bool(True), 
    hltProcName = cms.untracked.string(hltProcess), 
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

# Filter rarely used gen particle types
# Taken from: https://cmssdt.cern.ch/SDT/doxygen/CMSSW_8_0_9/doc/html/d1/d1f/prunedGenParticles__cfi_8py_source.html
process.filteredPrunedGenParticles = cms.EDProducer(
  "GenParticlePruner",
  src = cms.InputTag("prunedGenParticles"),
  select = cms.vstring(
    "drop  *", # this is the default
    # keep leptons, with history
    #OPT "++keep abs(pdgId) == 11 || abs(pdgId) == 13 || abs(pdgId) == 15",
    #OPT --> keep leptons
    "keep abs(pdgId) == 11 || abs(pdgId) == 13 || abs(pdgId) == 15",
    # keep neutrinos
    "keep abs(pdgId) == 12 || abs(pdgId) == 14 || abs(pdgId) == 16",
    # drop the shower part of the history
    #"drop   status == 2",
    # keep gamma above 10 GeV (or all prompt) and its first parent
    #"+keep pdgId == 22 && status == 1 && (pt > 10 || isPromptFinalState())",
    #OPT  # keep first parent of electrons above 3 GeV (or prompt)
    #OPT  "+keep abs(pdgId) == 11 && status == 1 && (pt > 3 || isPromptFinalState())",
    # but keep keep taus with decays
    #"keep++ abs(pdgId) == 15",
    #remove pythia8 garbage
    #"drop  status > 30 && status < 70 ",
    #remove pythia8 garbage
    #"drop  pdgId == 21 && pt < 5",
    # but remove again gluons in the inheritance chain
    #"drop   status == 2 && abs(pdgId) == 21",
    # keep VIP(articles)s
    #OPT "keep abs(pdgId) == 23 || abs(pdgId) == 24 || abs(pdgId) == 25 || abs(pdgId) == 6 || abs(pdgId) == 37 ",
    #OPT --> keep VIP(articles)s (Z,W,h,t,H+) and their (grand)children, except for Z, also parent of top
    #"keep++ abs(pdgId) == 24 || abs(pdgId) == 25 || abs(pdgId) == 37 ",
    #"keep abs(pdgId) == 23 ",
    #"+keep++ abs(pdgId) == 6",
    # keep K0
    #OPT "keep abs(pdgId) == 310 && abs(eta) < 2.5 && pt > 1 ",
    # keep heavy flavour quarks for parton-based jet flavour
    #"keep (4 <= abs(pdgId) <= 5) & (status = 2 || status = 11 || status = 71 || status = 72)",
    #OPT  # keep light-flavour quarks and gluons for parton-based jet flavour
    #OPT  "keep (1 <= abs(pdgId) <= 3 || pdgId = 21) & (status = 2 || status = 11 || status = 71 || status = 72) && pt>5", 
    #OPT  # keep b and c hadrons for hadron-based jet flavour
    #OPT  "keep (400 < abs(pdgId) < 600) || (4000 < abs(pdgId) < 6000)",
    #OPT  # additional c hadrons for jet fragmentation studies
    #OPT  "keep abs(pdgId) = 10411 || abs(pdgId) = 10421 || abs(pdgId) = 10413 || abs(pdgId) = 10423 || abs(pdgId) = 20413 || abs(pdgId) = 20423 || abs(pdgId) = 10431 || abs(pdgId) = 10433 || abs(pdgId) = 20433", 
    #OPT  # additional b hadrons for jet fragmentation studies
    #OPT  "keep abs(pdgId) = 10511 || abs(pdgId) = 10521 || abs(pdgId) = 10513 || abs(pdgId) = 10523 || abs(pdgId) = 20513 || abs(pdgId) = 20523 || abs(pdgId) = 10531 || abs(pdgId) = 10533 || abs(pdgId) = 20533 || abs(pdgId) = 10541 || abs(pdgId) = 10543 || abs(pdgId) = 20543", 
    #keep SUSY particles
    #"keep (1000001 <= abs(pdgId) <= 1000039 ) || ( 2000001 <= abs(pdgId) <= 2000015)",
    # keep protons 
    "keep pdgId = 2212",
    #keep event summary (status=3 for pythia6, 21 <= status <= 29 for pythia8)
    "keep status == 3 || ( 21 <= status <= 29) || ( 11 <= status <= 19)",
    #keep event summary based on status flags
    "keep isHardProcess() || fromHardProcessFinalState() || fromHardProcessDecayed() || fromHardProcessBeforeFSR() || (statusFlags().fromHardProcess() && statusFlags().isLastCopy())",
    )
    #select = cms.vstring(
    #"keep *",
    #"drop (1 <= abs(pdgId) <= 4)"
    #)
)

### Including ntuplizer 
process.options.allowUnscheduled = cms.untracked.bool(True)
process.load("Analysis.B2GAnaFW.b2gedmntuples_cff")

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
    "keep *_jetsAK4CHS_*_*",
    "keep *_jetsAK8CHS_*_*",
    "keep *_subjetsAK8CHS_*_*",
    "keep *_jetKeysAK4CHS_*_*",
    "keep *_jetKeysAK8CHS_*_*",
    "keep *_subjetKeysAK8CHS_*_*",
    "keep *_jetsAK4Puppi_*_*",
    "keep *_jetsAK8Puppi_*_*",
    "keep *_subjetsAK8Puppi_*_*",
    "keep *_jetKeysAK4Puppi_*_*",
    "keep *_jetKeysAK8Puppi_*_*",
    "keep *_subjetKeysAK8Puppi_*_*",
    #"keep *_eventShape*_*_*",
    "keep *_*_*centrality*_*",
    "keep *_metFull_*_*",
    "keep *_metNoHF_*_*",
    "keep *_METUserData*_trigger*_*",
    "keep *_eventInfo_*_*",
    "keep *_electronKeys_*_*",   
    "keep *_photonKeys_*_*",   
    "keep *_muonKeys_*_*",
    "keep *_TriggerUserData*_trigger*_*",
    "keep *_fixedGridRhoFastjetAll_*_*",
    "keep *_eventUserData_*_*",
    "keep *_BadPFMuonFilter_*_*",
    "keep *_BadChargedCandidateFilter_*_*",
    "keep *_eeBadScFilter_*_*"
    ),
    dropMetaData = cms.untracked.string('ALL'),
    )


  ### keep NoHF jets if needed:
if( options.useNoHFMET ):
  process.edmNtuplesOut.outputCommands+=('keep *_jetsAK4NoHF_*_*',)

if "MC" in options.DataProcessing: 
  process.edmNtuplesOut.outputCommands+=(
      'keep *_generator_*_*',
      #"keep *_genPart_*_*",
      "keep *_filteredPrunedGenParticles_*_*",
      "keep *_genJetsAK8*_*_*",
      "keep *_genJetsAK8SoftDrop*_*_*",      
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

open('B2GEntupleFileDump.py','w').write(process.dumpPython())
