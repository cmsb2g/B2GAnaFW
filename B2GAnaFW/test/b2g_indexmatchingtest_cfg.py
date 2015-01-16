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
                 1,# default value: process all events
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.int,
                 'Number of events to process')

options.register('sample',
                 'file:/afs/cern.ch/work/d/decosa/public/DMtt/miniAOD_Phys14.root',
                 #'/store/mc/Phys14DR/TprimeJetToTH_allHdecays_M1200GeV_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/20000/94117DA2-009A-E411-9DFB
-002590494CB2.root',
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.string,
                 'Sample to analyze')


options.register('outputLabel',
                 'analysisTTDM.root',
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
ak8jetLabel = 'patJetsSlimmedJetsAK8BTagged'
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


process = cms.Process("ttDManalysisEDMNtuples")

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

from PhysicsTools.PatAlgos.tools.jetTools import addJetCollection
addJetCollection(
    process,
    postfix   = "",
    labelName = 'SlimmedJetsAK8BTagged',
    jetSource = cms.InputTag('slimmedJetsAK8'),
    trackSource = cms.InputTag('unpackedTracksAndVertices'),
    pfCandidates = cms.InputTag('packedPFCandidates'),
    pvSource = cms.InputTag('unpackedTracksAndVertices'),
    jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'Type-2'),
    btagDiscriminators = [      'combinedSecondaryVertexBJetTags', 'combinedInclusiveSecondaryVertexV2BJetTags'     ]
    ,algo= 'AK', rParam = 0.8
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

### Subjet b-tagging
# Build jet collection with reco tools
from RecoJets.JetProducers.ak5PFJets_cfi import ak5PFJets
# Gen jets, with neutrinos subtracted
process.packedGenParticlesForJetsNoNu = cms.EDFilter("CandPtrSelector", src = cms.InputTag("packedGenParticles"), cut = cms.string("abs(pdgId) != 12 && abs(pdgId) !=
14 && abs(pdgId) != 16"))
from RecoJets.JetProducers.ak4GenJets_cfi import ak4GenJets
process.ak8GenJetsNoNuPruned = ak4GenJets.clone(
    rParam = cms.double(0.8),
    src = cms.InputTag("packedGenParticlesForJetsNoNu"),
    usePruning = cms.bool(True),
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets"),
    nFilt = cms.int32(3), # number of subjets, exclusive
    zcut = cms.double(0.1),
    rcut_factor = cms.double(0.5)
)
# Reco jets, with CHS
process.chs = cms.EDFilter("CandPtrSelector", src = cms.InputTag("packedPFCandidates"), cut = cms.string("fromPV()>0"))
process.ak8PFJetsCHS = ak5PFJets.clone(
                                      src = 'chs',
                                      rParam = cms.double(0.8),
                                      jetPtMin = cms.double(100.0),
                                      jetAlgorithm = cms.string("AntiKt")
                                      )

process.ak8PFJetsCHSPruned = ak5PFJets.clone(
                                        src = 'chs',
                                        rParam = cms.double(0.8),
                                        jetPtMin = cms.double(100.0),
                                        jetAlgorithm = cms.string("AntiKt"),
                                        usePruning = cms.bool(True),
                                        useExplicitGhosts = cms.bool(True),
                                        jetCollInstanceName = cms.string("SubJets"),
                                        writeCompound = cms.bool(True),
                                        nFilt = cms.int32(3), # number of subjets, exclusive
                                        zcut = cms.double(0.1),
                                        rcut_factor = cms.double(0.5)
                                        )
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
# Run Pat tools on AK8 jets
addJetCollection(
    process,
    labelName = 'AK8PFCHSPruned',
    jetSource = cms.InputTag('ak8PFJetsCHSPruned'),
    trackSource = cms.InputTag('unpackedTracksAndVertices'),
    btagDiscriminators = ['None'],
    jetCorrections = ('AK7PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None'),
    genJetCollection = cms.InputTag('ak8GenJetsNoNuPruned'),
)
# fix references
getattr(process,'patJetPartonMatchAK8PFCHSPruned').matched = cms.InputTag('prunedGenParticles')
process.patJetGenJetMatchAK8PFCHSPruned.matched = "slimmedGenJets"
process.patJetPartonMatchAK8PFCHSPruned.matched = "prunedGenParticles"
process.patJetCorrFactorsAK8PFCHSPruned.primaryVertices = "unpackedTracksAndVertices"
# Run Pat tools on AK8 subjets
addJetCollection(
    process,
    labelName = 'AK8PFCHSPrunedSubjets',
    jetSource = cms.InputTag('ak8PFJetsCHSPruned','SubJets'),
    trackSource = cms.InputTag('unpackedTracksAndVertices'),
    algo = 'ak',  # needed for subjet flavor clustering
    rParam = 0.8, # needed for subjet flavor clustering
    pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
    pfCandidates = cms.InputTag('packedPFCandidates'),
    btagDiscriminators = bTagDiscriminators,
    jetCorrections = ('AK7PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None'),
    genJetCollection = cms.InputTag('ak8GenJetsNoNuPruned','SubJets'),
)
# fix references
getattr(process,'patJetPartonMatchAK8PFCHSPrunedSubjets').matched = cms.InputTag('prunedGenParticles')
process.patJetGenJetMatchAK8PFCHSPrunedSubjets.matched = "slimmedGenJets"
process.patJetPartonMatchAK8PFCHSPrunedSubjets.matched = "prunedGenParticles"
process.patJetCorrFactorsAK8PFCHSPrunedSubjets.primaryVertices = "unpackedTracksAndVertices"
if hasattr(process,'pfInclusiveSecondaryVertexFinderTagInfosAK8PFCHSPrunedSubjets'):
    getattr(process,'pfInclusiveSecondaryVertexFinderTagInfosAK8PFCHSPrunedSubjets').extSVCollection = cms.InputTag('slimmedSecondaryVertices')
getattr(process,'patJetsAK8PFCHSPrunedSubjets').addAssociatedTracks = cms.bool(False) # no track collection present in MiniAOD
getattr(process,'patJetsAK8PFCHSPrunedSubjets').addJetCharge = cms.bool(False)        # no track collection present in MiniAOD
# Establish references between PATified fat jets and subjets using the BoostedJetMerger
process.selectedPatJetsAK8PFCHSPrunedPacked = cms.EDProducer("BoostedJetMerger",
    jetSrc=cms.InputTag("selectedPatJetsAK8PFCHSPruned"),
    subjetSrc=cms.InputTag("selectedPatJetsAK8PFCHSPrunedSubjets")
)

### Add Nsubjettiness
from RecoJets.JetProducers.nJettinessAdder_cfi import Njettiness
process.Njettiness = Njettiness.clone(
    src = cms.InputTag("ak8PFJetsCHS"),
    cone = cms.double(0.8)
    )

process.patJetsAK8PFCHS.userData.userFloats.src += ['Njettiness:tau1','Njettiness:tau2','Njettiness:tau3']

#$#$#$#$#$#$#$#$#$#
#   TOP TAG JETS  #
# patJetsCMSTopTagCHS
from RecoJets.JetProducers.PFJetParameters_cfi import *
from RecoJets.JetProducers.AnomalousCellParameters_cfi import *
from RecoJets.JetProducers.CATopJetParameters_cfi import *
process.cmsTopTagCHS = cms.EDProducer(
    "CATopJetProducer",
    PFJetParameters.clone( src = cms.InputTag('chs'),
                           doAreaFastjet = cms.bool(True),
                           doRhoFastjet = cms.bool(False),
                           jetPtMin = cms.double(200.0)
                           ),
    AnomalousCellParameters,
    CATopJetParameters.clone( jetCollInstanceName = cms.string("SubJets"),
                              verbose = cms.bool(False),
                              algorithm = cms.int32(1), # 0 = KT, 1 = CA, 2 = anti-KT
                              tagAlgo = cms.int32(0), #0=legacy top
                              useAdjacency = cms.int32(2), # modified adjacency
                              centralEtaCut = cms.double(2.5), # eta for defining "central" jets
                              sumEtBins = cms.vdouble(0,1600,2600), # sumEt bins over which cuts vary. vector={bin 0 lower bound, bin 1 lower bound, ...}
                              rBins = cms.vdouble(0.8,0.8,0.8), # Jet distance paramter R. R values depend on sumEt bins.
                              ptFracBins = cms.vdouble(0.05,0.05,0.05), # minimum fraction of central jet pt for subjets (deltap)
                              deltarBins = cms.vdouble(0.19,0.19,0.19), # Applicable only if useAdjacency=1. deltar adjacency values for each sumEtBin
                              nCellBins = cms.vdouble(1.9,1.9,1.9),
                            ),
    jetAlgorithm = cms.string("CambridgeAachen"),
    rParam = cms.double(0.8),
    writeCompound = cms.bool(True)
    )
process.CATopTagInfos = cms.EDProducer("CATopJetTagger",
                                    src = cms.InputTag("cmsTopTagCHS"),
                                    TopMass = cms.double(171),
                                    TopMassMin = cms.double(0.),
                                    TopMassMax = cms.double(250.),
                                    WMass = cms.double(80.4),
                                    WMassMin = cms.double(0.0),
                                    WMassMax = cms.double(200.0),
                                    MinMassMin = cms.double(0.0),
                                    MinMassMax = cms.double(200.0),
                                    verbose = cms.bool(False)
                                    )
addJetCollection(
    process,
    labelName = 'CMSTopTagCHS',
    jetSource = cms.InputTag('cmsTopTagCHS'),
    jetCorrections = ('AK7PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
    trackSource = cms.InputTag('unpackedTracksAndVertices'),
    pvSource = cms.InputTag("unpackedTracksAndVertices"),
    btagDiscriminators = ['combinedSecondaryVertexBJetTags'],
    getJetMCFlavour = False
    )
process.patJetPartonMatchCMSTopTagCHS.matched='prunedGenParticles'
process.patJetCorrFactorsCMSTopTagCHS.primaryVertices = "unpackedTracksAndVertices"
process.patJetGenJetMatchCMSTopTagCHS.matched = 'slimmedGenJets'
process.patJetPartonMatchCMSTopTagCHS.matched = 'prunedGenParticles'
#process.jetTracksAssociatorAtVertexCMSTopTagCHS=process.ak5JetTracksAssociatorAtVertexPF.clone(jets = cms.InputTag('cmsTopTagCHS'), coneSize = 0.8)
process.secondaryVertexTagInfosCMSTopTagCHS.trackSelection.jetDeltaRMax = cms.double(0.8) # default is 0.3
process.secondaryVertexTagInfosCMSTopTagCHS.vertexCuts.maxDeltaRToJetAxis = cms.double(0.8) # default is 0.5
process.combinedSecondaryVertexCMSTopTagCHS= process.combinedSecondaryVertex.clone()
process.combinedSecondaryVertexCMSTopTagCHS.trackSelection.jetDeltaRMax = cms.double(0.8)
process.combinedSecondaryVertexCMSTopTagCHS.trackPseudoSelection.jetDeltaRMax = cms.double(0.8)
process.combinedSecondaryVertexBJetTagsCMSTopTagCHS.jetTagComputer = cms.string('combinedSecondaryVertexCMSTopTagCHS')
process.patJetsCMSTopTagCHS.addTagInfos = True
process.patJetsCMSTopTagCHS.tagInfoSources = cms.VInputTag(
    cms.InputTag('CATopTagInfos')
    )

addJetCollection(
    process,
    labelName = 'CMSTopTagCHSSubjets',
    jetSource = cms.InputTag('cmsTopTagCHS','SubJets'),
    jetCorrections = ('AK7PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
    trackSource = cms.InputTag('unpackedTracksAndVertices'),
    pvSource = cms.InputTag("unpackedTracksAndVertices"),
    btagDiscriminators = ['combinedSecondaryVertexBJetTags'],
    getJetMCFlavour = False,
    )
process.patJetPartonMatchCMSTopTagCHSSubjets.matched='prunedGenParticles'
process.patJetCorrFactorsCMSTopTagCHSSubjets.primaryVertices = "unpackedTracksAndVertices"
process.patJetGenJetMatchCMSTopTagCHSSubjets.matched = 'slimmedGenJets'
process.patJetPartonMatchCMSTopTagCHSSubjets.matched = 'prunedGenParticles'

process.patJetsCMSTopTagCHSPacked = cms.EDProducer("BoostedJetMerger",
    jetSrc=cms.InputTag("patJetsCMSTopTagCHS" ),
    subjetSrc=cms.InputTag("patJetsCMSTopTagCHSSubjets")
      )

#$#$#$#$#$#$#$#$#$#

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
    cut = cms.string("pt > 1")
)

process.skimmedCMSTOPTAGSubJets = cms.EDFilter(
    "CandViewSelector",
    src = cms.InputTag("selectedPatJetsCMSTopTagCHSSubjets"),
    cut = cms.string("pt > 1")
)


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

process.patjetUserData = cms.EDProducer(
    'PatJetUserData',
    jetLabel  = cms.InputTag("selectedPatJetsAK8PFCHS"),
    packedjetLabel  = cms.InputTag("selectedPatJetsAK8PFCHSPrunedPacked"),
    subjetLabel  = cms.InputTag("selectedPatJetsAK8PFCHSPrunedSubjets"),
    pv        = cms.InputTag(pvLabel),
    ### TTRIGGER ###
    triggerResults = cms.InputTag(triggerResultsLabel,"","HLT"),
    triggerSummary = cms.InputTag(triggerSummaryLabel,"","HLT"),
    hltJetFilter       = cms.InputTag("hltSixCenJet20L1FastJet"),
    hltPath            = cms.string("HLT_QuadJet60_DiJet20_v6"),
    hlt2reco_deltaRmax = cms.double(0.2),
    doSubjets          = cms.bool(True)
)

process.cmstoptagjetUserData = cms.EDProducer(
    'PatJetUserData',
    jetLabel  = cms.InputTag("patJetsCMSTopTagCHS"),
    packedjetLabel  = cms.InputTag("patJetsCMSTopTagCHSPacked"),
    subjetLabel  = cms.InputTag("patJetsCMSTopTagCHSSubjets"),
    pv        = cms.InputTag(pvLabel),
    ### TTRIGGER ###
    triggerResults = cms.InputTag(triggerResultsLabel,"","HLT"),
    triggerSummary = cms.InputTag(triggerSummaryLabel,"","HLT"),
    hltJetFilter       = cms.InputTag("hltSixCenJet20L1FastJet"),
    hltPath            = cms.string("HLT_QuadJet60_DiJet20_v6"),
    hlt2reco_deltaRmax = cms.double(0.2),
    doSubjets          = cms.bool(True)
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
    process.skimmedPatElectrons +
    process.skimmedPatMuons +
    process.skimmedPatJets +
    process.skimmedPatJetsAK8 +
    process.skimmedPatSubJetsAK8+
    process.skimmedCMSTOPTAGSubJets+
    process.skimmedPatMET +
    process.eventShapePFVars +
    process.eventShapePFJetVars +
    process.centrality
)

#process.analysisPath+=process.jetFilter

#process.analysisPath+=process.egmGsfElectronIDSequence
process.analysisPath+=process.muonUserData
process.analysisPath+=process.jetUserData
process.analysisPath+=process.patjetUserData
process.analysisPath+=process.cmstoptagjetUserData
process.analysisPath+=process.electronUserData
process.analysisPath+=process.genPart
process.analysisPath+=process.muons
process.analysisPath+=process.electrons
process.analysisPath+=process.jetsAK4
process.analysisPath+=process.jetsAK8
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
