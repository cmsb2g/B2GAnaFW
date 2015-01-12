### *****************************************************************************************
### Usage:
###
### cmsRun topplusdmanaEDMntuples_cfg.py maxEvts=N sample="mySample/sample.root" version="71" outputLabel="myoutput"
###
### Default values for the options are set:
### maxEvts     = -1
### sample      = 'file:/scratch/decosa/ttDM/testSample/tlbsm_53x_v3_mc_10_1_qPV.root'
### outputLabel = 'analysisTTDM.root'
### *****************************************************************************************
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as opts

options = opts.VarParsing ('analysis')

options.register('maxEvts',
                 100,# default value: process all events
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.int,
                 'Number of events to process')

options.register('sample',
                 #'file:/afs/cern.ch/work/d/decosa/public/DMtt/miniAOD_Phys14.root',
                 #/TprimeJetToTH_allHdecays_M1200GeV_Tune4C_13TeV-madgraph-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM
                 '/store/mc/Phys14DR/TprimeJetToTH_allHdecays_M1200GeV_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/20000/94117DA2-009A-E411-9DFB-002590494CB2.root',
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.string,
                 'Sample to analyze')


options.register('outputLabel',
                 'analysisTprime.root',
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.string,
                 'Output label')

options.register('isData',
                 False,
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.bool,
                 'Is data?')

options.register('LHE',
                 True,
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.bool,
                 'Keep LHEProducts')

options.parseArguments()

if(options.isData):options.LHE = False

    
###inputTag labels
muLabel  = 'slimmedMuons'
elLabel  = 'slimmedElectrons'
jLabel = 'slimmedJets'
jLabelAK8 = 'slimmedJetsAK8'
#ak8jetLabel = 'patJetsSlimmedJetsAK8BTagged'
ak8jetLabel = 'selectedPatJetsAK8PFCHS' 
ak8subjetLabel = 'selectedPatJetsAK8PFCHSPrunedSubjets' 
pvLabel  = 'offlineSlimmedPrimaryVertices'
convLabel = 'reducedEgamma:reducedConversions'
particleFlowLabel = 'packedPFCandidates'    
metLabel = 'slimmedMETs'

triggerResultsLabel = "TriggerResults"
triggerSummaryLabel = "hltTriggerSummaryAOD"
hltMuonFilterLabel       = "hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f40QL3crIsoRhoFiltered0p15"
hltPathLabel             = "HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL"
hltElectronFilterLabel  = "hltL1sL1Mu3p5EG12ORL1MuOpenEG12L3Filtered8"
lheLabel = "source"


process = cms.Process("b2gEDMNtuples")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.categories.append('HLTrigReport')
### Output Report
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
### Number of maximum events to process
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvts) )
### Source file
process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring(
        options.sample
        )
)

process.load("PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff")
process.load("Configuration.EventContent.EventContent_cff")
process.load('Configuration.StandardSequences.Geometry_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag as customiseGlobalTag
process.GlobalTag = customiseGlobalTag(process.GlobalTag, globaltag = 'auto:startup_GRun')
process.GlobalTag.connect   = 'frontier://FrontierProd/CMS_COND_31X_GLOBALTAG'
process.GlobalTag.pfnPrefix = cms.untracked.string('frontier://FrontierProd/')
for pset in process.GlobalTag.toGet.value():
    pset.connect = pset.connect.value().replace('frontier://FrontierProd/', 'frontier://FrontierProd/')
#   Fix for multi-run processing:
process.GlobalTag.RefreshEachRun = cms.untracked.bool( False )
process.GlobalTag.ReconnectEachRun = cms.untracked.bool( False )
    


###
### AK8 jets with subjet b-tagging
###

## Filter out neutrinos from packed GenParticles
process.packedGenParticlesForJetsNoNu = cms.EDFilter("CandPtrSelector", 
    src = cms.InputTag("packedGenParticles"), 
    cut = cms.string("abs(pdgId) != 12 && abs(pdgId) != 14 && abs(pdgId) != 16")
    )

## Fat GenJets
from RecoJets.JetProducers.ak4GenJets_cfi import ak4GenJets
process.ak8GenJetsNoNu = ak4GenJets.clone(
    rParam = cms.double(0.8),
    src = cms.InputTag("packedGenParticlesForJetsNoNu")
    )

## Pruned fat GenJets (two jet collections are produced, fat jets and subjets)
from RecoJets.JetProducers.SubJetParameters_cfi import SubJetParameters
process.ak8GenJetsNoNuPruned = ak4GenJets.clone(
   SubJetParameters,
   rParam = cms.double(0.8),
   src = cms.InputTag("packedGenParticlesForJetsNoNu"),
   usePruning = cms.bool(True),
   writeCompound = cms.bool(True),
   jetCollInstanceName=cms.string("SubJets")
   )

## Select charged hadron subtracted packed PF candidates
process.pfCHS = cms.EDFilter("CandPtrSelector", 
    src = cms.InputTag("packedPFCandidates"), 
    cut = cms.string("fromPV")
    )

## Fat PFJets
from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets
process.ak8PFJetsCHS = ak4PFJets.clone(
   rParam = cms.double(0.8),
   src = cms.InputTag("pfCHS"),
   doAreaFastjet = cms.bool(True),
   jetPtMin = cms.double(50.)
   )

## Pruned fat PFJets (two jet collections are produced, fat jets and subjets)
from RecoJets.JetProducers.ak5PFJetsPruned_cfi import ak5PFJetsPruned
process.ak8PFJetsCHSPruned = ak5PFJetsPruned.clone(
   rParam = cms.double(0.8),
   src = cms.InputTag("pfCHS"),
   doAreaFastjet = cms.bool(True),
   writeCompound = cms.bool(True),
   jetCollInstanceName=cms.string("SubJets"),
   jetPtMin = cms.double(50.)
   )

#################################################
## Make PAT jets
#################################################

#for Inclusive Vertex Finder
process.load("RecoBTag/Configuration/RecoBTag_cff")
process.load('RecoVertex/AdaptiveVertexFinder/inclusiveVertexing_cff')
process.inclusiveVertexFinder.tracks = cms.InputTag("unpackedTracksAndVertices")
process.inclusiveVertexFinder.primaryVertices = cms.InputTag("unpackedTracksAndVertices")
process.trackVertexArbitrator.tracks = cms.InputTag("unpackedTracksAndVertices")
process.trackVertexArbitrator.primaryVertices = cms.InputTag("unpackedTracksAndVertices")

#new input for impactParameterTagInfos, softleptons, IVF
process.impactParameterTagInfos.jetTracks = cms.InputTag("jetTracksAssociatorAtVertexSlimmedJetsAK8BTagged")
process.impactParameterTagInfos.primaryVertex = cms.InputTag("unpackedTracksAndVertices")
process.inclusiveVertexFinder.primaryVertices = cms.InputTag("unpackedTracksAndVertices")
process.trackVertexArbitrator.primaryVertices = cms.InputTag("unpackedTracksAndVertices")
process.softPFMuonsTagInfos.primaryVertex = cms.InputTag("unpackedTracksAndVertices")
process.softPFElectronsTagInfos.primaryVertex = cms.InputTag("unpackedTracksAndVertices")
process.softPFMuonsTagInfos.jets = cms.InputTag("patJetsSlimmedJetsAK8BTagged")
process.softPFElectronsTagInfos.jets = cms.InputTag("patJetsSlimmedJetsAK8BTagged") 
process.inclusiveSecondaryVertexFinderTagInfosV2 = process.inclusiveSecondaryVertexFinderTagInfos.clone()
process.inclusiveSecondaryVertexFinderTagInfosV2.trackSelection.qualityClass = cms.string('any')

## b-tag discriminators
bTagDiscriminators = [
    'trackCountingHighEffBJetTags',
    'trackCountingHighPurBJetTags',
    'jetProbabilityBJetTags',
    'jetBProbabilityBJetTags',
    'simpleSecondaryVertexHighEffBJetTags',
    'simpleSecondaryVertexHighPurBJetTags',
    'combinedSecondaryVertexBJetTags',
    'combinedInclusiveSecondaryVertexV2BJetTags'
    ]

from PhysicsTools.PatAlgos.tools.jetTools import addJetCollection

addJetCollection(
    process,
    labelName = 'AK8PFCHS',
    jetSource = cms.InputTag('ak8PFJetsCHS'),
    algo = 'ak',  # needed for jet flavor clustering
    rParam = 0.8, # needed for jet flavor clustering
    pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
    pfCandidates = cms.InputTag('packedPFCandidates'),
    #svSource = cms.InputTag('slimmedSecondaryVertices'),
    btagDiscriminators = bTagDiscriminators,
    jetCorrections = ('AK4PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None'),
    genJetCollection = cms.InputTag('ak8GenJetsNoNu')
    )

getattr(process,'patJetPartons').particles = cms.InputTag('prunedGenParticles')
getattr(process,'patJetPartonMatchAK8PFCHS').matched = cms.InputTag('prunedGenParticles')

#if hasattr(process,'pfInclusiveSecondaryVertexFinderTagInfosAK8PFCHS'):
#  getattr(process,'pfInclusiveSecondaryVertexFinderTagInfosAK8PFCHS').extSVCollection = cms.InputTag('slimmedSecondaryVertices')

getattr(process,'patJetsAK8PFCHS').addAssociatedTracks = cms.bool(False) # needs to be disabled since there is no track collection present in MiniAOD
getattr(process,'patJetsAK8PFCHS').addJetCharge = cms.bool(False)        # needs to be disabled since there is no track collection present in MiniAOD

process.jetTracksAssociatorAtVertexAK8PFCHS.tracks = cms.InputTag("unpackedTracksAndVertices")

### PATify pruned fat jets
addJetCollection(
    process,
    labelName = 'AK8PFCHSPruned',
    jetSource = cms.InputTag('ak8PFJetsCHSPruned'),
    btagDiscriminators = ['None'],
    jetCorrections = ('AK4PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None'),
    genJetCollection = cms.InputTag('ak8GenJetsNoNu'),
    getJetMCFlavour = False # jet flavor disabled
    )
getattr(process,'patJetPartonMatchAK8PFCHSPruned').matched = cms.InputTag('prunedGenParticles')

### PATify pruned subjets
addJetCollection(
    process,
    labelName = 'AK8PFCHSPrunedSubjets',
    jetSource = cms.InputTag('ak8PFJetsCHSPruned','SubJets'),
    algo = 'ak',  # needed for subjet flavor clustering
    rParam = 0.8, # needed for subjet flavor clustering
    pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
    pfCandidates = cms.InputTag('packedPFCandidates'),
    #svSource = cms.InputTag('slimmedSecondaryVertices'),
    btagDiscriminators = bTagDiscriminators,
    jetCorrections = ('AK4PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None'),
    genJetCollection = cms.InputTag('ak8GenJetsNoNuPruned','SubJets'),
    #Apparently only in CMSSW_7_3_X explicitJTA = True,  # needed for subjet b tagging
    #Apparently only in CMSSW_7_3_X svClustering = True, # needed for subjet b tagging
    #Apparently only in CMSSW_7_3_X fatJets=cms.InputTag('ak8PFJetsCHS'),             # needed for subjet flavor clustering
    #Apparently only in CMSSW_7_3_X groomedFatJets=cms.InputTag('ak8PFJetsCHSPruned') # needed for subjet flavor clustering
    )

if hasattr( process, 'jetTracksAssociatorAtVertex' + 'AK8PFCHSPrunedSubjets' ):
  process.jetTracksAssociatorAtVertexAK8PFCHSPrunedSubjets.tracks = cms.InputTag("unpackedTracksAndVertices")
#  from RecoJets.JetAssociationProducers.ak4aTA_cff import ak4JetTracksAssociatorExplicit
#  m = 'jetTracksAssociatorAtVertex' + 'AK8PFCHSPrunedSubjets'
#  print 'Switching ' + m + ' to explicit jet-track association'
#  setattr( process, m, ak4JetTracksAssociatorExplicit.clone(
#    jets = getattr(getattr(process,m),'jets'),
#    tracks = cms.InputTag("unpackedTracksAndVertices")
#    ) 
#    )

getattr(process,'patJetPartonMatchAK8PFCHSPrunedSubjets').matched = cms.InputTag('prunedGenParticles')
#if hasattr(process,'pfInclusiveSecondaryVertexFinderTagInfosAK8PFCHSPrunedSubjets'):
#  getattr(process,'pfInclusiveSecondaryVertexFinderTagInfosAK8PFCHSPrunedSubjets').extSVCollection = cms.InputTag('slimmedSecondaryVertices')
getattr(process,'patJetsAK8PFCHSPrunedSubjets').addAssociatedTracks = cms.bool(False) # needs to be disabled since there is no track collection present in MiniAOD
getattr(process,'patJetsAK8PFCHSPrunedSubjets').addJetCharge = cms.bool(False)        # needs to be disabled since there is no track collection present in MiniAOD

## Establish references between PATified fat jets and subjets using the BoostedJetMerger
process.selectedPatJetsAK8PFCHSPrunedPacked = cms.EDProducer("BoostedJetMerger",
    jetSrc=cms.InputTag("selectedPatJetsAK8PFCHSPruned"),
    subjetSrc=cms.InputTag("selectedPatJetsAK8PFCHSPrunedSubjets")
    )

from PhysicsTools.PatAlgos.tools.pfTools import *
## Adapt primary vertex collection
adaptPVs(process, pvCollection=cms.InputTag('offlineSlimmedPrimaryVertices'))

addJetCollection(
    process,
    postfix   = "",
    labelName = 'SlimmedJetsAK8BTagged',
    jetSource = cms.InputTag('slimmedJetsAK8'),
    trackSource = cms.InputTag('unpackedTracksAndVertices'),
    pfCandidates = cms.InputTag('packedPFCandidates'), 
    pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
    jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'Type-2'),
    btagDiscriminators = ['combinedSecondaryVertexBJetTags', 'combinedInclusiveSecondaryVertexV2BJetTags'],
    algo= 'AK', 
    rParam = 0.8
    )

process.patJetFlavourAssociation.jets = "slimmedJetsAK8"
process.patJetFlavourAssociation.rParam = 0.8 

#adjust MC matching
process.patJetGenJetMatchSlimmedJetsAK8BTagged.matched = "slimmedGenJets"
process.patJetPartonMatchSlimmedJetsAK8BTagged.matched = "prunedGenParticles"
process.patJetPartons.particles = "prunedGenParticles"

#adjust PV used for Jet Corrections
process.patJetCorrFactorsSlimmedJetsAK8BTagged.primaryVertices = "unpackedTracksAndVertices"

#adjust JTA cone size 
process.jetTracksAssociatorAtVertexSlimmedJetsAK8BTagged.coneSize = 0.8 

process.load('PhysicsTools.PatAlgos.slimming.unpackedTracksAndVertices_cfi')
process.combinedSecondaryVertex.trackMultiplicityMin = 1 #silly sv, uses un filtered tracks.. i.e. any pt

##
## Jet substructure variables as user floats
##
#process.patJetsSlimmedJetsAK8BTagged.userData.userFloats.src = []
#process.selectedPatJetsSlimmedJetsAK8BTagged.cut = cms.string("pt > 100")
#from RecoJets.Configuration.RecoPFJets_cff import ak8PFJetsCHSPruned, ak8PFJetsCHSFiltered, ak8PFJetsCHSTrimmed
#process.ak8PFJetsCHSPruned   = ak8PFJetsCHSPruned.clone()
#process.ak8PFJetsCHSTrimmed  = ak8PFJetsCHSTrimmed.clone()
#process.ak8PFJetsCHSFiltered = ak8PFJetsCHSFiltered.clone()
#
#process.load("RecoJets.JetProducers.ak8PFJetsCHS_groomingValueMaps_cfi")
#process.ak8PFJetsCHSPrunedLinks.src  = cms.InputTag("slimmedJetsAK8")
#process.ak8PFJetsCHSTrimmedLinks.src  = cms.InputTag("slimmedJetsAK8")
#process.ak8PFJetsCHSFilteredLinks.src = cms.InputTag("slimmedJetsAK8")
#
#process.patJetsSlimmedJetsAK8BTagged.userData.userFloats.src += ['ak8PFJetsCHSPrunedLinks','ak8PFJetsCHSTrimmedLinks','ak8PFJetsCHSFilteredLinks']
#
#process.load('RecoJets.JetProducers.nJettinessAdder_cfi')
#process.NjettinessAK8 = process.Njettiness.clone(src=cms.InputTag("slimmedJetsAK8"),)
#process.NjettinessAK8.cone = cms.double(0.8)
#process.patJetsSlimmedJetsAK8BTagged.userData.userFloats.src += ['NjettinessAK8:tau1','NjettinessAK8:tau2','NjettinessAK8:tau3']


### Selected leptons and jets
process.skimmedPatMuons = cms.EDFilter(
    "PATMuonSelector",
    src = cms.InputTag(muLabel),
    cut = cms.string("pt > 30 && abs(eta) < 2.4")
    )

process.skimmedPatElectrons = cms.EDFilter(
    "PATElectronSelector",
    src = cms.InputTag(elLabel),
    cut = cms.string("pt > 30 && abs(eta) < 2.5")
    )

process.skimmedPatMET = cms.EDFilter(
    "PATMETSelector",
    src = cms.InputTag(metLabel),
    cut = cms.string("")
    )


process.skimmedPatJets = cms.EDFilter(
    "PATJetSelector",
    src = cms.InputTag(jLabel),
    cut = cms.string(" pt > 25 && abs(eta) < 4.")
    )


process.skimmedPatJetsAK8 = cms.EDFilter(
    "CandViewSelector",
    src = cms.InputTag(ak8jetLabel),
    cut = cms.string("pt > 100 && abs(eta) < 4.")    
    )

process.skimmedPatSubJetsAK8 = cms.EDFilter(
    "CandViewSelector",
    src = cms.InputTag(ak8subjetLabel),
    cut = cms.string("")
    )


#process.skimmedPatJetsAK8 = cms.EDFilter(
#    "PATJetSelector",
#    src = cms.InputTag(jLabelAK8),
#    cut = cms.string("pt > 100 && abs(eta) < 4.")
#)



### Asking for at least 2 jets satisfying the selection above 
process.jetFilter = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("skimmedPatJets"),
    minNumber = cms.uint32(2),
    filter = cms.bool(True)
    )



### Electron ID

# Load tools and function definitions
#from PhysicsTools.SelectorUtils.tools.vid_id_tools import *

# # Turn on VID producer
### switchOnVIDElectronIdProducer(process)
#process.load("RecoEgamma.ElectronIdentification.egmGsfElectronIDs_cfi")
#process.egmGsfElectronIDs.physicsObjectSrc = cms.InputTag('slimmedElectrons')

#from PhysicsTools.SelectorUtils.centralIDRegistry import central_id_registry
#process.egmGsfElectronIDSequence = cms.Sequence(process.egmGsfElectronIDs)

# Define which IDs we want to produce
# Each of these two example IDs contains all four standard 
# cut-based ID working points (only two WP of the PU20bx25 are actually used here).

#my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_PHYS14_PU20bx25_V0_miniAOD_cff']

#Add them to the VID producer
#for idmod in my_id_modules:
#    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)


process.muonUserData = cms.EDProducer(
    'MuonUserData',
    muonLabel = cms.InputTag("skimmedPatMuons"),
    pv        = cms.InputTag(pvLabel),
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
    jetLabel  = cms.InputTag("skimmedPatJets"),
    pv        = cms.InputTag(pvLabel),
    ### TTRIGGER ###
    triggerResults = cms.InputTag(triggerResultsLabel,"","HLT"),
    triggerSummary = cms.InputTag(triggerSummaryLabel,"","HLT"),
    hltJetFilter       = cms.InputTag("hltSixCenJet20L1FastJet"),
    hltPath            = cms.string("HLT_QuadJet60_DiJet20_v6"),
    hlt2reco_deltaRmax = cms.double(0.2),
    )


process.jetUserDataAK8 = cms.EDProducer(
    'JetUserData',
    jetLabel  = cms.InputTag("skimmedPatJetsAK8"),
    pv        = cms.InputTag(pvLabel),
    ### TTRIGGER ###
    triggerResults = cms.InputTag(triggerResultsLabel,"","HLT"),
    triggerSummary = cms.InputTag(triggerSummaryLabel,"","HLT"),
    hltJetFilter       = cms.InputTag("hltSixCenJet20L1FastJet"),
    hltPath            = cms.string("HLT_QuadJet60_DiJet20_v6"),
    hlt2reco_deltaRmax = cms.double(0.2)
    )

process.subjetUserDataAK8 = cms.EDProducer(
    'JetUserData',
    jetLabel  = cms.InputTag("skimmedPatSubJetsAK8"),
    pv        = cms.InputTag(pvLabel),
    ### TTRIGGER ###
    triggerResults = cms.InputTag(triggerResultsLabel,"","HLT"),
    triggerSummary = cms.InputTag(triggerSummaryLabel,"","HLT"),
    hltJetFilter       = cms.InputTag("hltSixCenJet20L1FastJet"),
    hltPath            = cms.string("HLT_QuadJet60_DiJet20_v6"),
    hlt2reco_deltaRmax = cms.double(0.2)
    )

process.electronUserData = cms.EDProducer(
    'ElectronUserData',
    eleLabel = cms.InputTag("skimmedPatElectrons"),
    pv        = cms.InputTag(pvLabel),
    conversion        = cms.InputTag(convLabel),
    triggerResults = cms.InputTag(triggerResultsLabel),
    triggerSummary = cms.InputTag(triggerSummaryLabel),
    hltElectronFilter  = cms.InputTag(hltElectronFilterLabel),  ##trigger matching code to be fixed!
    hltPath             = cms.string("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL"),
    #electronVetoIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V0-miniAOD-standalone-veto"),
    #electronTightIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V0-miniAOD-standalone-tight"),
    )



from PhysicsTools.CandAlgos.EventShapeVars_cff import *
process.eventShapePFVars = pfEventShapeVars.clone()
process.eventShapePFVars.src = cms.InputTag(particleFlowLabel)

process.eventShapePFJetVars = pfEventShapeVars.clone()
process.eventShapePFJetVars.src = cms.InputTag("skimmedPatJets")

process.centrality = cms.EDProducer("CentralityUserData",
    src = cms.InputTag("skimmedPatJets")
    )                                    

### Including ntuplizer 
process.load("B2GAnaFW.B2GAnaFW.b2gedmntuples_cff")

process.options.allowUnscheduled = cms.untracked.bool(True)


### definition of Analysis sequence
process.analysisPath = cms.Path(
    #process.selectedPatJetsAK8PFCHS + 
    #process.selectedPatJetsAK8PFCHSPrunedSubjets + 
    #process.selectedPatJetsAK8PFCHSPrunedPacked + 
    process.skimmedPatElectrons +
    process.skimmedPatMuons +
    process.skimmedPatJets +
    process.skimmedPatJetsAK8 +
    process.skimmedPatSubJetsAK8 +
    process.skimmedPatMET +
    process.eventShapePFVars +
    process.eventShapePFJetVars +
    process.centrality
    )

#process.analysisPath+=process.jetFilter

#process.analysisPath+=process.egmGsfElectronIDSequence
process.analysisPath+=process.muonUserData
process.analysisPath+=process.jetUserData
#process.analysisPath+=process.jetUserDataAK8
#process.analysisPath+=process.subjetUserDataAK8
process.analysisPath+=process.electronUserData
process.analysisPath+=process.genPart
process.analysisPath+=process.muons
process.analysisPath+=process.electrons
process.analysisPath+=process.jets
process.analysisPath+=process.jetsAK8
process.analysisPath+=process.subjetsAK8
process.analysisPath+=process.met

### Creating the filter path to use in order to select events
process.filterPath = cms.Path(
    process.jetFilter
    )


### keep info from LHEProducts if they are stored in PatTuples
if(options.LHE):
  process.LHEUserData = cms.EDProducer("LHEUserData",
  lheLabel = cms.InputTag("source")
  )

process.analysisPath+=process.LHEUserData

process.edmNtuplesOut.outputCommands+=('keep *_*LHE*_*_*',)

process.edmNtuplesOut.outputCommands+=('keep LHEEventProduct_*_*_*',)


### end LHE products     

process.edmNtuplesOut.SelectEvents = cms.untracked.PSet(
    SelectEvents = cms.vstring('filterPath')
    )


process.fullPath = cms.Schedule(
    process.analysisPath,
    process.filterPath
    )

#process.edmNtuplesOut.SelectEvents ='filterPath'


process.endPath = cms.EndPath(process.edmNtuplesOut)

open('junk.py','w').write(process.dumpPython())
