import FWCore.ParameterSet.Config as cms
import copy

#### basic set of variables which are commons to all the objects
basic =  cms.EDProducer(
    "CandViewNtpProducer",
    src=cms.InputTag("skimmedPatMuons"),
    lazyParser=cms.untracked.bool(True),
    prefix=cms.untracked.string("basic"),
    eventInfo=cms.untracked.bool(False),
    variables = cms.VPSet(
    cms.PSet(
    tag = cms.untracked.string("Mass"),
    quantity = cms.untracked.string("mass")
    ),
    cms.PSet(
    tag = cms.untracked.string("Pt"),
    quantity = cms.untracked.string("pt")
    ),
    cms.PSet(
    tag = cms.untracked.string("Eta"),
    quantity = cms.untracked.string("eta")
    ),
    cms.PSet(
    tag = cms.untracked.string("Y"),
    quantity = cms.untracked.string("rapidity")
    ),
    cms.PSet(
    tag = cms.untracked.string("Phi"),
    quantity = cms.untracked.string("phi")
    ),
    cms.PSet(
    tag = cms.untracked.string("E"),
    quantity = cms.untracked.string("energy")
    ),
    cms.PSet(
    tag = cms.untracked.string("Charge"),
    quantity = cms.untracked.string("charge")
    ),
    )
    )

met =  cms.EDProducer(
    "CandViewNtpProducer",
    src=cms.InputTag("skimmedPatMET"),
    lazyParser=cms.untracked.bool(True),
    prefix=cms.untracked.string("met"),
    eventInfo=cms.untracked.bool(False),
    variables = cms.VPSet(
    cms.PSet(
    tag = cms.untracked.string("Pt"),
    quantity = cms.untracked.string("pt")
    ),
    cms.PSet(
    tag = cms.untracked.string("Px"),
    quantity = cms.untracked.string("px")
    ),
    cms.PSet(
    tag = cms.untracked.string("Py"),
    quantity = cms.untracked.string("py")
    ),
    cms.PSet(
    tag = cms.untracked.string("Phi"),
    quantity = cms.untracked.string("phi")
    ),
    )
    )

### muon variables
muonVars = (
#   cms.PSet(
#        tag = cms.untracked.string("dB"),
#        quantity = cms.untracked.string("userFloat('dB')")
#   ),
#   cms.PSet(
#        tag = cms.untracked.string("dBPV2D"),
#        quantity = cms.untracked.string("userFloat('dBPV2D')")
#   ),
#   cms.PSet(
#        tag = cms.untracked.string("dBPV3D"),
#        quantity = cms.untracked.string("userFloat('dBPV3D')")
#   ),
#   cms.PSet(
#        tag = cms.untracked.string("dBBS2D"),
#        quantity = cms.untracked.string("userFloat('dBBS2D')")
#   ),
#   cms.PSet(
#        tag = cms.untracked.string("dBBS3D"),
#        quantity = cms.untracked.string("userFloat('dBBS3D')")
#   ),
   cms.PSet(
        tag = cms.untracked.string("Key"),
        quantity = cms.untracked.string("originalObjectRef().key()")
   ),
   cms.PSet(
        tag = cms.untracked.string("Iso04"),
        quantity = cms.untracked.string("userFloat('iso04')")
   ),
   cms.PSet(
        tag = cms.untracked.string("D0"),
        quantity = cms.untracked.string("dB")
   ),
   cms.PSet(
        tag = cms.untracked.string("D0err"),
        quantity = cms.untracked.string("edB")
   ),
   cms.PSet(
        tag = cms.untracked.string("Dxy"),
        quantity = cms.untracked.string("userFloat('dxy')")
   ),
   cms.PSet(
        tag = cms.untracked.string("Dxyerr"),
        quantity = cms.untracked.string("userFloat('dxyErr')")
   ),
   cms.PSet(
        tag = cms.untracked.string("Dz"),
        quantity = cms.untracked.string("userFloat('dz')")
        ),
   cms.PSet(
        tag = cms.untracked.string("Dzerr"),
        quantity = cms.untracked.string("userFloat('dzErr')")
        ),
### the following variables need have track embedded in the pat::muon
    cms.PSet(
    tag = cms.untracked.string("IsLooseMuon"),
    quantity = cms.untracked.string("isLooseMuon")
    ),
    cms.PSet(
    tag = cms.untracked.string("IsSoftMuon"),
    quantity = cms.untracked.string("userFloat('isSoftMuon')")
    ),
    cms.PSet(
    tag = cms.untracked.string("IsTightMuon"),
    quantity = cms.untracked.string("userFloat('isTightMuon')")
    ),
## variables used in ID
## https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId#Tight_Muon_selection
   ### LOOSE
   cms.PSet(
    tag = cms.untracked.string("IsPFMuon"),
    quantity = cms.untracked.string("isPFMuon")
   ),
   cms.PSet(
    tag = cms.untracked.string("IsGlobalMuon"),
    quantity = cms.untracked.string("isGlobalMuon")
   ),
   cms.PSet(
    tag = cms.untracked.string("IsTrackerMuon"),
    quantity = cms.untracked.string("isTrackerMuon")
   ),
   ### TIGHT
   cms.PSet(
    tag = cms.untracked.string("GlbTrkNormChi2"),
    quantity = cms.untracked.string("? globalTrack.isNonnull ? globalTrack.normalizedChi2 : -900")
   ),
   cms.PSet(
    tag = cms.untracked.string("NumberValidMuonHits"),
    quantity = cms.untracked.string("? globalTrack.isNonnull ? globalTrack.hitPattern.numberOfValidMuonHits : -900")
   ),
   cms.PSet(
    tag = cms.untracked.string("NumberMatchedStations"),
    quantity = cms.untracked.string("numberOfMatchedStations")
   ),
   cms.PSet(
    tag = cms.untracked.string("NumberValidPixelHits"),
    quantity = cms.untracked.string("? innerTrack.isNonnull ? innerTrack.hitPattern.numberOfValidPixelHits : -900")
   ),
   cms.PSet(
    tag = cms.untracked.string("NumberTrackerLayers"),
    quantity = cms.untracked.string("? track.isNonnull ? track.hitPattern.trackerLayersWithMeasurement : -900")
   ),
   ### SOFT
   cms.PSet(
    tag = cms.untracked.string("NumberOfValidTrackerHits"),
    quantity = cms.untracked.string("? innerTrack.isNonnull ? innerTrack.hitPattern.numberOfValidTrackerHits : -900")
   ),
   cms.PSet(
    tag = cms.untracked.string("NumberOfPixelLayers"),
    quantity = cms.untracked.string("? innerTrack.isNonnull ? innerTrack.hitPattern.pixelLayersWithMeasurement : -900")
   ),
   cms.PSet(
    tag = cms.untracked.string("InTrkNormChi2"),
    quantity = cms.untracked.string("? innerTrack.isNonnull ? innerTrack.normalizedChi2 : -900")
   ),
## variables used in isolation
## https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId#Accessing_PF_Isolation_from_reco
   cms.PSet(
    tag = cms.untracked.string("SumChargedHadronPt"),
    quantity = cms.untracked.string("pfIsolationR04().sumChargedHadronPt")
   ),
   cms.PSet(
    tag = cms.untracked.string("SumNeutralHadronPt"),
    quantity = cms.untracked.string("pfIsolationR04().sumNeutralHadronEt")
   ),
   cms.PSet(
    tag = cms.untracked.string("SumPhotonPt"),
    quantity = cms.untracked.string("pfIsolationR04().sumPhotonEt")
   ),
   cms.PSet(
    tag = cms.untracked.string("SumPUPt"),
    quantity = cms.untracked.string("pfIsolationR04().sumPUPt")
   ),
### genLepton
   cms.PSet(
    tag = cms.untracked.string("GenMuonY"),
    quantity = cms.untracked.string("? genParticleRef.isNonnull ? genLepton.rapidity : -900")
   ),
   cms.PSet(
    tag = cms.untracked.string("GenMuonEta"),
    quantity = cms.untracked.string("? genParticleRef.isNonnull ? genLepton.eta : -900")
   ),
   cms.PSet(
    tag = cms.untracked.string("GenMuonPhi"),
    quantity = cms.untracked.string("? genParticleRef.isNonnull ? genLepton.phi : -900")
   ),
   cms.PSet(
    tag = cms.untracked.string("GenMuonPt"),
    quantity = cms.untracked.string("? genParticleRef.isNonnull ? genLepton.pt : -900")
   ),
   cms.PSet(
    tag = cms.untracked.string("GenMuonE"),
    quantity = cms.untracked.string("? genParticleRef.isNonnull ? genLepton.energy : -900")
   ),
   cms.PSet(
    tag = cms.untracked.string("GenMuonCharge"),
    quantity = cms.untracked.string("? genParticleRef.isNonnull ? genLepton.charge : -900")
   ),
### trigger matching
 #  cms.PSet(
 #   tag = cms.untracked.string("HLTmuonDeltaR"),
 #   quantity = cms.untracked.string("userFloat('HLTmuonDeltaR')")
 #  ),
 #  cms.PSet(
 #   tag = cms.untracked.string("HLTmuonPt"),
 #   quantity = cms.untracked.string("userFloat('HLTmuonPt')")
 #  ),
 #  cms.PSet(
 #   tag = cms.untracked.string("HLTmuonEta"),
 #   quantity = cms.untracked.string("userFloat('HLTmuonEta')")
 #  ),
 #  cms.PSet(
 #   tag = cms.untracked.string("HLTmuonPhi"),
 #   quantity = cms.untracked.string("userFloat('HLTmuonPhi')")
 #  ),
 #  cms.PSet(
 #   tag = cms.untracked.string("HLTmuonE"),
 #   quantity = cms.untracked.string("userFloat('HLTmuonE')")
 #  ),
)

### jet variables
jetVars = (
### B-TAGGING
    cms.PSet(
     tag = cms.untracked.string("CSV"),
     quantity = cms.untracked.string("bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags')")
    ),
    cms.PSet(
     tag = cms.untracked.string("CSVV1"),
     quantity = cms.untracked.string("bDiscriminator('pfCombinedSecondaryVertexV1BJetTags')")
    ),
### GEN PARTON
    cms.PSet(
     tag = cms.untracked.string("GenPartonY"),
     quantity = cms.untracked.string("? genParticleRef.isNonnull ? genParton.rapidity : -900")
    ),
    cms.PSet(
     tag = cms.untracked.string("GenPartonEta"),
     quantity = cms.untracked.string("? genParticleRef.isNonnull ? genParton.eta : -900")
    ),
    cms.PSet(
     tag = cms.untracked.string("GenPartonPhi"),
     quantity = cms.untracked.string("? genParticleRef.isNonnull ? genParton.phi : -900")
    ),
    cms.PSet(
     tag = cms.untracked.string("GenPartonPt"),
     quantity = cms.untracked.string("? genParticleRef.isNonnull ? genParton.pt : -900")
    ),
    cms.PSet(
     tag = cms.untracked.string("GenPartonE"),
     quantity = cms.untracked.string("? genParticleRef.isNonnull ? genParton.energy : -900")
    ),
    cms.PSet(
     tag = cms.untracked.string("GenPartonCharge"),
     quantity = cms.untracked.string("? genParticleRef.isNonnull ? genParton.charge : -900")
    ),
###
    cms.PSet(
     tag = cms.untracked.string("PartonFlavour"),
     quantity = cms.untracked.string("partonFlavour")
    ),
    cms.PSet(
     tag = cms.untracked.string("HadronFlavour"),
     quantity = cms.untracked.string("hadronFlavour")
    ),
### GEN JET
    cms.PSet(
     tag = cms.untracked.string("GenJetY"),
     quantity = cms.untracked.string("? genJetFwdRef.isNonnull ? genJet.rapidity : -900")
    ),
    cms.PSet(
     tag = cms.untracked.string("GenJetEta"),
     quantity = cms.untracked.string("? genJetFwdRef.isNonnull ? genJet.eta : -900")
    ),
    cms.PSet(
     tag = cms.untracked.string("GenJetPhi"),
     quantity = cms.untracked.string("? genJetFwdRef.isNonnull ? genJet.phi : -900")
    ),
    cms.PSet(
     tag = cms.untracked.string("GenJetPt"),
     quantity = cms.untracked.string("? genJetFwdRef.isNonnull ? genJet.pt : -900")
    ),
    cms.PSet(
     tag = cms.untracked.string("GenJetE"),
     quantity = cms.untracked.string("? genJetFwdRef.isNonnull ? genJet.energy : -900")
    ),
    cms.PSet(
     tag = cms.untracked.string("GenJetCharge"),
     quantity = cms.untracked.string("? genJetFwdRef.isNonnull ? genJet.charge : -900")
    ),
### TRIGGER MATHING
   # cms.PSet(
   #  tag = cms.untracked.string("HLTjetEta"),
   #  quantity = cms.untracked.string("userFloat('HLTjetEta')")
   # ),
   # cms.PSet(
   #  tag = cms.untracked.string("HLTjetPhi"),
   #  quantity = cms.untracked.string("userFloat('HLTjetPhi')")
   # ),
   # cms.PSet(
   #  tag = cms.untracked.string("HLTjetPt"),
   #  quantity = cms.untracked.string("userFloat('HLTjetPt')")
   # ),
   # cms.PSet(
   #  tag = cms.untracked.string("HLTjetE"),
   #  quantity = cms.untracked.string("userFloat('HLTjetE')")
   # ),
   # cms.PSet(
   #  tag = cms.untracked.string("HLTjetDeltaR"),
   #  quantity = cms.untracked.string("userFloat('HLTjetDeltaR')")
   # ),
### CONSTITUENTS
    cms.PSet(
     tag = cms.untracked.string("muonMultiplicity"),
     quantity = cms.untracked.string("?isPFJet || isJPTJet ? muonMultiplicity : -1")
    ),
    cms.PSet(
     tag = cms.untracked.string("PhotonEnergy"),
     quantity = cms.untracked.string("? isPFJet ? photonEnergy : -1")
    ),
    cms.PSet(
     tag = cms.untracked.string("ElectronEnergy"),
     quantity = cms.untracked.string("? isPFJet ? electronEnergy : -1")
    ),
    cms.PSet(
     tag = cms.untracked.string("MuonEnergy"),
     quantity = cms.untracked.string("? isPFJet ? muonEnergy : -1")
    ),
    cms.PSet(
     tag = cms.untracked.string("HFHadronEnergy"),
     quantity = cms.untracked.string("? isPFJet ? HFHadronEnergy : -1")
    ),
    cms.PSet(
     tag = cms.untracked.string("HFEMEnergy"),
     quantity = cms.untracked.string("? isPFJet ? HFEMEnergy : -1")
    ),
    cms.PSet(
     tag = cms.untracked.string("ChargedHadronMultiplicity"),
     quantity = cms.untracked.string("? isPFJet ? chargedHadronMultiplicity : -1")
    ),
    cms.PSet(
      tag = cms.untracked.string("numberOfDaughters"),
      quantity = cms.untracked.string("? isPFJet ? numberOfDaughters : -1")
       ),
   cms.PSet(
      tag = cms.untracked.string("chargedMultiplicity"),
      quantity = cms.untracked.string("? isPFJet ? chargedMultiplicity : -1")
    ),
    cms.PSet(
     tag = cms.untracked.string("neutralHadronMultiplicity"),
     quantity = cms.untracked.string("? isPFJet ? neutralHadronMultiplicity : -1")
    ),
    cms.PSet(
     tag = cms.untracked.string("neutralHadronEnergy"),
     quantity = cms.untracked.string("? isPFJet ? neutralHadronEnergy : -1")
    ),
    cms.PSet(
    tag = cms.untracked.string("neutralEmEnergy"),
    quantity = cms.untracked.string("? isPFJet ? neutralEmEnergy : -1"),
    ),
    cms.PSet(
    tag = cms.untracked.string("chargedEmEnergy"),
    quantity = cms.untracked.string("? isPFJet ? chargedEmEnergy : -1"),
    ),
    cms.PSet(
    tag = cms.untracked.string("chargedHadronEnergy"),
    quantity = cms.untracked.string("? isPFJet ? chargedHadronEnergy : -1"),
    ),
     cms.PSet(
     tag = cms.untracked.string("photonMultiplicity"),
     quantity = cms.untracked.string("? isPFJet ? photonMultiplicity : -1")
    ),
    cms.PSet(
     tag = cms.untracked.string("electronMultiplicity"),
     quantity = cms.untracked.string("? isPFJet ? electronMultiplicity : -1")
    ),
    cms.PSet(
     tag = cms.untracked.string("HFHadronMultiplicity"),
     quantity = cms.untracked.string("? isPFJet ? HFHadronMultiplicity : -1")
    ),
    cms.PSet(
     tag = cms.untracked.string("HFEMMultiplicity"),
     quantity = cms.untracked.string("? isPFJet ? HFEMMultiplicity : -1")
    ),
    cms.PSet(
     tag = cms.untracked.string("ChargeMuEnergy"),
     quantity = cms.untracked.string("? isPFJet ? chargedMuEnergy : -1")
    ),
    cms.PSet(
     tag = cms.untracked.string("neutralMultiplicity"),
     quantity = cms.untracked.string("? isPFJet ? neutralMultiplicity : -1")
    ),
#### FOR LEPTON MATCHING 
#    cms.PSet(
#     tag = cms.untracked.string("pfKeys"),
#     quantity = cms.untracked.string("? hasUserData('pfKeys') ? userData('pfKeys') : 0")
#    ),
#    cms.PSet(
#     tag = cms.untracked.string("dummy"),
#     quantity = cms.untracked.string("? hasUserData('dummy') ? 1 : 0")
#    ),    
#### FOR JEC
    cms.PSet(
        tag = cms.untracked.string("jecFactor0"),
        quantity = cms.untracked.string("jecFactor(0)")
        ),
    cms.PSet(
        tag = cms.untracked.string("jetArea"),
        quantity = cms.untracked.string("jetArea")
        ),
#### FOR SYSTEMATICS
    cms.PSet(
     tag = cms.untracked.string("SmearedPt"),
     quantity = cms.untracked.string("userFloat('SmearedPt')")
    ),
    cms.PSet(
     tag = cms.untracked.string("SmearedPEta"),
     quantity = cms.untracked.string("userFloat('SmearedPEta')")
    ),
    cms.PSet(
     tag = cms.untracked.string("SmearedPhi"),
     quantity = cms.untracked.string("userFloat('SmearedPhi')")
    ),
    cms.PSet(
     tag = cms.untracked.string("SmearedE"),
     quantity = cms.untracked.string("userFloat('SmearedE')")
    ),
    cms.PSet(
     tag = cms.untracked.string("JERup"),
     quantity = cms.untracked.string("userFloat('JERup')")
    ),
    cms.PSet(
     tag = cms.untracked.string("JERdown"),
     quantity = cms.untracked.string("userFloat('JERdown')")
    ),
)

jetKeys = cms.EDProducer(
    "JetKeyProducer",
    jetLabel = cms.InputTag("ak4PFJetsCHS")
    )

genPartVars = (
    cms.PSet(
    tag = cms.untracked.string("ID"),
    quantity = cms.untracked.string("pdgId")
    ),
    cms.PSet(
    tag = cms.untracked.string("Status"),
    quantity = cms.untracked.string("status")
    ),
    cms.PSet(
    tag = cms.untracked.string("MomID"),
    quantity = cms.untracked.string("?numberOfMothers>0 ? mother(0).pdgId : -900")
    ),
    )



### jet variables
jetAK8Vars = (
#### SUBSTRUCTURE
     cms.PSet(
        tag = cms.untracked.string("vSubjetIndex0"),
        quantity = cms.untracked.string("? hasUserInt('VSubjet0') ? userInt('VSubjet0') : -1 ")
        ),
     cms.PSet(
        tag = cms.untracked.string("vSubjetIndex1"),
        quantity = cms.untracked.string("? hasUserInt('VSubjet1') ? userInt('VSubjet1') : -1 ")
        ),
     cms.PSet(
        tag = cms.untracked.string("topSubjetIndex0"),
        quantity = cms.untracked.string("? hasUserInt('TopSubjet0') ? userInt('TopSubjet0') : -1 ")
        ),
     cms.PSet(
        tag = cms.untracked.string("topSubjetIndex1"),
        quantity = cms.untracked.string("? hasUserInt('TopSubjet1') ? userInt('TopSubjet1') : -1 ")
        ),
     cms.PSet(
        tag = cms.untracked.string("topSubjetIndex2"),
        quantity = cms.untracked.string("? hasUserInt('TopSubjet2') ? userInt('TopSubjet2') : -1 ")
        ),
     cms.PSet(
        tag = cms.untracked.string("topSubjetIndex3"),
        quantity = cms.untracked.string("? hasUserInt('TopSubjet3') ? userInt('TopSubjet3') : -1 ")
        ),
     cms.PSet(
        tag = cms.untracked.string("tau1"),
        quantity = cms.untracked.string("userFloat('NjettinessAK8:tau1')")
        ),
     cms.PSet(
        tag = cms.untracked.string("tau2"),
        quantity = cms.untracked.string("userFloat('NjettinessAK8:tau2')")
        ),
     cms.PSet(
        tag = cms.untracked.string("tau3"),
        quantity = cms.untracked.string("userFloat('NjettinessAK8:tau3')")
        ),
     cms.PSet(
        tag = cms.untracked.string("softDropMass"),
        quantity = cms.untracked.string("userFloat('ak8PFJetsCHSSoftDropMass')")
        ),
     cms.PSet(
        tag = cms.untracked.string("trimmedMass"),
        quantity = cms.untracked.string("userFloat('ak8PFJetsCHSTrimmedMass')")
        ),
     cms.PSet(
        tag = cms.untracked.string("prunedMass"),
        quantity = cms.untracked.string("userFloat('ak8PFJetsCHSPrunedMass')")
        ),
     cms.PSet(
        tag = cms.untracked.string("filteredMass"),
        quantity = cms.untracked.string("userFloat('ak8PFJetsCHSFilteredMass')")
        ),
     cms.PSet(
        tag = cms.untracked.string("topMass"),
        quantity = cms.untracked.string("? hasTagInfo('caTop') ? tagInfo('caTop').properties().topMass : 0")
        ),
     cms.PSet(
        tag = cms.untracked.string("wMass"),
        quantity = cms.untracked.string("? hasTagInfo('caTop') ? tagInfo('caTop').properties().wMass : 0")
        ),
     cms.PSet(
        tag = cms.untracked.string("nSubJets"),
        quantity = cms.untracked.string("? hasTagInfo('caTop') ? tagInfo('caTop').properties().nSubJets : 0")
        ),
     cms.PSet(
        tag = cms.untracked.string("minmass"),
        quantity = cms.untracked.string("? hasTagInfo('caTop') ? tagInfo('caTop').properties().minMass : 0")
        ),
)


### copying the muon set of variables from basic,
### adding the set of variable which are related to muons only
muons = copy.deepcopy(basic)
muons.variables += muonVars
muons.prefix = cms.untracked.string("mu")
muons.src = cms.InputTag("muonUserData")
#muons.src = cms.InputTag("skimmedPatMuons")

###electrons
electronVars = (
    ###Cut-based ID variables
   cms.PSet(
        tag = cms.untracked.string("Key"),
        quantity = cms.untracked.string("originalObjectRef().key()")
   ),
  cms.PSet(
        tag = cms.untracked.string("Iso03"),
        quantity = cms.untracked.string("userFloat('iso03')")
   ),
   cms.PSet(
        tag = cms.untracked.string("D0"),
        quantity = cms.untracked.string("userFloat('d0')")
   ),
   cms.PSet(
        tag = cms.untracked.string("Dz"),
        quantity = cms.untracked.string("userFloat('dz')")
        ),
   cms.PSet(
        tag = cms.untracked.string("dEtaIn"),
        quantity = cms.untracked.string("deltaEtaSuperClusterTrackAtVtx")
        ),
   cms.PSet(
        tag = cms.untracked.string("dPhiIn"),
        quantity = cms.untracked.string("deltaPhiSuperClusterTrackAtVtx")
        ),
   cms.PSet(
        tag = cms.untracked.string("HoE"),
        quantity = cms.untracked.string("hcalOverEcal")
        ),
   cms.PSet(
        tag = cms.untracked.string("full5x5siee"),
        quantity = cms.untracked.string("full5x5_sigmaIetaIeta")
        ),
   cms.PSet(
        tag = cms.untracked.string("ooEmooP"),
        quantity = cms.untracked.string("userFloat('ooEmooP')")
        ),
   cms.PSet(
        tag = cms.untracked.string("missHits"),
        quantity = cms.untracked.string("userInt('missHits')")
        ),
   cms.PSet(
        tag = cms.untracked.string("hasMatchedConVeto"),
        quantity = cms.untracked.string("userInt('hasMatchConv')")
        ),
  cms.PSet(
        tag = cms.untracked.string("isEB"),
        quantity = cms.untracked.string("isEB")
        ),
  cms.PSet(
        tag = cms.untracked.string("isVeto"),
        quantity = cms.untracked.string("userInt('isVeto')")
        ),
  cms.PSet(
        tag = cms.untracked.string("isLoose"),
        quantity = cms.untracked.string("userInt('isLoose')")
        ),
  cms.PSet(
        tag = cms.untracked.string("isTight"),
        quantity = cms.untracked.string("userInt('isTight')")
        ),
  cms.PSet(
        tag = cms.untracked.string("isMedium"),
        quantity = cms.untracked.string("userInt('isMedium')")
        ),
  cms.PSet(
        tag = cms.untracked.string("scEta"),
        quantity = cms.untracked.string("superCluster().eta()")
        )

)


electrons = copy.deepcopy(basic)
electrons.variables += electronVars
electrons.prefix = cms.untracked.string("el")
electrons.src = cms.InputTag("electronUserData")


###jets
jetsAK4 = copy.deepcopy(basic)
jetsAK4.variables += jetVars
jetsAK4.prefix = cms.untracked.string("jetAK4")
jetsAK4.src = cms.InputTag("jetUserData")
jetKeysAK4 = copy.deepcopy( jetKeys )
jetKeysAK4.jetLabel = cms.InputTag("jetUserData")

###patjets
jetsAK8 = copy.deepcopy(basic)
jetsAK8.variables += jetVars
jetsAK8.variables += jetAK8Vars
jetsAK8.prefix = cms.untracked.string("jetAK8")
jetsAK8.src = cms.InputTag("jetUserDataAK8")
jetKeysAK8 = copy.deepcopy( jetKeys )
jetKeysAK8.jetLabel = cms.InputTag("jetUserDataAK8")

###subjetsAK8
subjetsAK8 = copy.deepcopy(basic)
subjetsAK8.variables += jetVars
subjetsAK8.prefix = cms.untracked.string("subjetAK8")
subjetsAK8.src = cms.InputTag("slimmedJetsAK8PFCHSSoftDropPacked", "SubJets")
subjetKeysAK8 = copy.deepcopy( jetKeys )
subjetKeysAK8.jetLabel = cms.InputTag("slimmedJetsAK8PFCHSSoftDropPacked", "SubJets")

###subjetsAK8
subjetsCmsTopTag = copy.deepcopy(basic)
subjetsCmsTopTag.variables += jetVars
subjetsCmsTopTag.prefix = cms.untracked.string("subjetCmsTopTag")
subjetsCmsTopTag.src = cms.InputTag("slimmedJetsCMSTopTagCHSPacked", "SubJets")
subjetsCmsTopTagKeys = copy.deepcopy( jetKeys )
subjetsCmsTopTagKeys.jetLabel = cms.InputTag("slimmedJetsCMSTopTagCHSPacked", "SubJets")


###genPart
genPart = copy.deepcopy(basic)
genPart.variables += genPartVars
genPart.prefix = cms.untracked.string("genPart")
genPart.src = cms.InputTag("prunedGenParticles")

###genJetsAK8
genJetsAK8 = copy.deepcopy(basic)
genJetsAK8.prefix = cms.untracked.string("genJetsAK8")
genJetsAK8.src = cms.InputTag("ak8GenJetsNoNu")

genJetsAK8SoftDrop = copy.deepcopy(basic)
genJetsAK8SoftDrop.prefix = cms.untracked.string("genJetsAK8SoftDrop")
genJetsAK8SoftDrop.src = cms.InputTag("ak8GenJetsNoNuSoftDrop")

###event variables
eventShapeVar = (
   cms.PSet(
        tag = cms.untracked.string("isotropy"),
        quantity = cms.untracked.string("isotropy")
   ),
   cms.PSet(
        tag = cms.untracked.string("circularity"),
        quantity = cms.untracked.string("circularity")
   ),
   cms.PSet(
        tag = cms.untracked.string("sphericity"),
        quantity = cms.untracked.string("sphericity")
   ),
   cms.PSet(
        tag = cms.untracked.string("aplanarity"),
        quantity = cms.untracked.string("aplanarity")
   ),
   cms.PSet(
        tag = cms.untracked.string("thrust"),
        quantity = cms.untracked.string("thrust")
   ),
)


### event info: EvtNumber, RunNumber, LumiBlock
eventInfo =  cms.EDProducer(
    "CandViewNtpProducer",
    src=cms.InputTag("skimmedPatMET"),
    lazyParser=cms.untracked.bool(True),
    prefix=cms.untracked.string("evtInfo"),
    eventInfo=cms.untracked.bool(True),
    variables = cms.VPSet()
    )

edmNtuplesOut = cms.OutputModule(
    "PoolOutputModule",
    fileName = cms.untracked.string('B2GEDMNtuple.root'),
    outputCommands = cms.untracked.vstring(
    "drop *",
    "keep *_genPart_*_*",
    "keep *_muons_*_*",
    "keep *_electrons_*_*",
    "keep *_jetsAK4_*_*",
    "keep *_jetsAK8_*_*",
    "keep *_eventShape*_*_*",
    "keep *_*_*centrality*_*",
    "keep *_met_*_*",
    "keep *_eventInfo_*_*",
    "keep *_subjetsAK8_*_*",
    "keep *_subjetsCmsTopTag*_*_*",
    "keep *_jetKeysAK4_*_*",
    "keep *_jetKeysAK8_*_*",
    "keep *_subjetKeysAK8_*_*",
    "keep *_subjetsCmsTopTagKeys_*_*",
    "keep *_genJetsAK8_*_*",
    "keep *_genJetsAK8SoftDrop_*_*",
    "keep *_TriggerUserData_trigger*_*",
    "keep *_fixedGridRhoFastjetAll_*_*",
    "keep *_eventUserData_*_*"
    ),
    dropMetaData = cms.untracked.string('ALL'),
    )
print "DONE STANDARD"
