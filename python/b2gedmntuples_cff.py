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

metFull =  cms.EDProducer(
    "CandViewNtpProducer",
    src=cms.InputTag("skimmedPatMET"),
    lazyParser=cms.untracked.bool(True),
    prefix=cms.untracked.string("metFull"),
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
        cms.PSet(
          tag = cms.untracked.string("uncorPt"),
          quantity = cms.untracked.string("uncorPt")
          ),
        cms.PSet(
          tag = cms.untracked.string("uncorPhi"),
          quantity = cms.untracked.string("uncorPhi")
          ),
        cms.PSet(
          tag = cms.untracked.string("uncorSumEt"),
          quantity = cms.untracked.string("uncorSumEt")
          ),
        ),
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
        tag = cms.untracked.string("MiniIso"),
        quantity = cms.untracked.string("userFloat('miniIso')")
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
        tag = cms.untracked.string("IsSoftMuon"),
        quantity = cms.untracked.string("userFloat('isSoftMuon')")
        ),
    cms.PSet(
        tag = cms.untracked.string("IsLooseMuon"),
        quantity = cms.untracked.string("userFloat('isLooseMuon')")
        ),
    cms.PSet(
        tag = cms.untracked.string("IsMediumMuon"),
        quantity = cms.untracked.string("userFloat('isMediumMuon')")
        ),
    cms.PSet(
        tag = cms.untracked.string("IsTightMuon"),
        quantity = cms.untracked.string("userFloat('isTightMuon')")
        ),
    cms.PSet(
        tag = cms.untracked.string("IsHighPtMuon"),
        quantity = cms.untracked.string("userFloat('isHighPtMuon')")
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
      tag = cms.untracked.string("CSVv2"),
      quantity = cms.untracked.string("bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags')")
      ),
    cms.PSet(
      tag = cms.untracked.string("DoubleBAK8"),
      quantity = cms.untracked.string("bDiscriminator('pfBoostedDoubleSecondaryVertexAK8BJetTags')")
      ),
    cms.PSet(
      tag = cms.untracked.string("DoubleBCA15"),
      quantity = cms.untracked.string("bDiscriminator('pfBoostedDoubleSecondaryVertexCA15BJetTags')")
      ),
    cms.PSet(
      tag = cms.untracked.string("CMVAv2"),
      quantity = cms.untracked.string("bDiscriminator('pfCombinedMVAV2BJetTags')")
      ),
    cms.PSet(
      tag = cms.untracked.string("CvsL"),
      quantity = cms.untracked.string("bDiscriminator('pfCombinedCvsLJetTags')")
      ),
    cms.PSet(
      tag = cms.untracked.string("CvsB"),
      quantity = cms.untracked.string("bDiscriminator('pfCombinedCvsBJetTags')")
      ),
    cms.PSet(
      tag = cms.untracked.string("CMVA"),
      quantity = cms.untracked.string("bDiscriminator('pfCombinedMVAV2BJetTags')")
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
   cms.PSet(
       tag = cms.untracked.string("neutralHadronEnergyFrac"),
       quantity = cms.untracked.string("? isPFJet ?neutralHadronEnergyFraction : -1")
       ),
   cms.PSet(
       tag = cms.untracked.string("neutralEmEnergyFrac"),
       quantity = cms.untracked.string("? isPFJet ?neutralEmEnergyFraction : -1")
       ),
   cms.PSet(
       tag = cms.untracked.string("chargedHadronEnergyFrac"),
       quantity = cms.untracked.string("? isPFJet ?chargedHadronEnergyFraction : -1")
       ),
   cms.PSet(
       tag = cms.untracked.string("muonEnergyFrac"),
       quantity = cms.untracked.string("? isPFJet ?muonEnergyFraction : -1")
       ),
   cms.PSet(
       tag = cms.untracked.string("chargedEmEnergyFrac"),
       quantity = cms.untracked.string("? isPFJet ?chargedEmEnergyFraction : -1")
       ),
   cms.PSet(
       tag = cms.untracked.string("chargedMultiplicity"),
       quantity = cms.untracked.string("? isPFJet ?chargedMultiplicity : -1")
       ),
   cms.PSet(
       tag = cms.untracked.string("NumConstituents"),
       quantity = cms.untracked.string("? isPFJet ? chargedMultiplicity + neutralMultiplicity : -1")
        ),

    #### FOR JEC
    cms.PSet(
        tag = cms.untracked.string("jecFactor0"),
        quantity = cms.untracked.string("jecFactor(0)")
        ),
#    cms.PSet(
#        tag = cms.untracked.string("jecFactorL1FastJet"),
#        quantity = cms.untracked.string("jecFactor('L1FastJet')")
#        ),
   cms.PSet(
        tag = cms.untracked.string("jecFactorL3Absolute"),
        quantity = cms.untracked.string("jecFactor('L3Absolute')")
        ),
#   cms.PSet(
#        tag = cms.untracked.string("jecFactorL3AbsoluteNum"),
#        quantity = cms.untracked.string("jecFactor(3)")
#        ),
#    cms.PSet(
#        tag = cms.untracked.string("jecFactorL2L3Residual"),
#        quantity = cms.untracked.string("jecFactor('L2L3Residual')")
#        ),
#    cms.PSet(
#        tag = cms.untracked.string("jecFactorL2L3ResidualNum"),
#        quantity = cms.untracked.string("jecFactor(4)")
#        ),
    cms.PSet(
        tag = cms.untracked.string("jetArea"),
        quantity = cms.untracked.string("jetArea")
        ),
    cms.PSet(
        tag = cms.untracked.string("nSV"),
        quantity = cms.untracked.string("? hasUserInt('nSV') ? userInt('nSV') : -999")
        ),
    cms.PSet(
        tag = cms.untracked.string("SV0mass"),
        quantity = cms.untracked.string("? hasUserFloat('SV0mass') ? userFloat('SV0mass') : -999")
        ),
    cms.PSet(
        tag = cms.untracked.string("SV1mass"),
        quantity = cms.untracked.string("? hasUserFloat('SV1mass') ? userFloat('SV1mass') : -999")
        ),
    )

#### FOR SYSTEMATICS
jetVarsForSys = (
    cms.PSet(
        tag = cms.untracked.string("jecUncertainty"),
        quantity = cms.untracked.string("userFloat('jecUncertainty')")
        ),
    cms.PSet(
        tag = cms.untracked.string("PtResolution"),
        quantity = cms.untracked.string("userFloat('PtResolution')")
        ),
    cms.PSet(
        tag = cms.untracked.string("JERSF"),
        quantity = cms.untracked.string("userFloat('JERSF')")
        ),
    cms.PSet(
        tag = cms.untracked.string("JERSFUp"),
        quantity = cms.untracked.string("userFloat('JERSFUp')")
        ),
    cms.PSet(
        tag = cms.untracked.string("JERSFDown"),
        quantity = cms.untracked.string("userFloat('JERSFDown')")
        ),
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
    )

jetVarsJEC = (
    ### FOR JEC
    cms.PSet(
        tag = cms.untracked.string("jecFactor0"),
        quantity = cms.untracked.string("jecFactor(0)")
        ),
    cms.PSet(
        tag = cms.untracked.string("jecFactorL1FastJet"),
        quantity = cms.untracked.string("jecFactor('L1FastJet')")
        ),
    cms.PSet(
        tag = cms.untracked.string("jecFactorL3Absolute"),
        quantity = cms.untracked.string("jecFactor('L3Absolute')")
        ),
    cms.PSet(
        tag = cms.untracked.string("jetArea"),
        quantity = cms.untracked.string("jetArea")
        ),
)

jetKeys = cms.EDProducer(
    "JetKeyProducer",
    jetLabel = cms.InputTag("jetUserData")
    )
electronKeys = cms.EDProducer(
    "SourceKeyProducer",
    srcLabel = cms.InputTag("electronUserData")
    )
photonKeys = cms.EDProducer(
    "SourceKeyProducer",
    srcLabel = cms.InputTag("photonUserData")
    )
muonKeys = cms.EDProducer(
    "SourceKeyProducer",
    srcLabel = cms.InputTag("muonUserData")
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
      tag = cms.untracked.string("Mom0ID"),
      quantity = cms.untracked.string("?numberOfMothers>0 ? mother(0).pdgId : -900")
      ),
    cms.PSet(
      tag = cms.untracked.string("Mom0Status"),
      quantity = cms.untracked.string("?numberOfMothers>0 ? mother(0).status : -900")
      ),
    cms.PSet(
      tag = cms.untracked.string("Mom1ID"),
      quantity = cms.untracked.string("?numberOfMothers>1 ? mother(1).pdgId : -900")
      ),
    cms.PSet(
      tag = cms.untracked.string("Mom1Status"),
      quantity = cms.untracked.string("?numberOfMothers>1 ? mother(1).status : -900")
      ),
    cms.PSet(
      tag = cms.untracked.string("Dau0ID"),
      quantity = cms.untracked.string("?numberOfDaughters>0 ? daughter(0).pdgId : -900")
      ),
    cms.PSet(
      tag = cms.untracked.string("Dau0Status"),
      quantity = cms.untracked.string("?numberOfDaughters>0 ? daughter(0).status : -900")
      ),
    cms.PSet(
      tag = cms.untracked.string("Dau1ID"),
      quantity = cms.untracked.string("?numberOfDaughters>1 ? daughter(1).pdgId : -900")
      ),
    cms.PSet(
      tag = cms.untracked.string("Dau1Status"),
      quantity = cms.untracked.string("?numberOfDaughters>1 ? daughter(1).status : -900")
      ),
    )



### jet variables
jetToolboxAK8Vars = (
#### SUBSTRUCTURE
     cms.PSet(
        tag = cms.untracked.string("vSubjetIndex0"),
        quantity = cms.untracked.string("? hasUserInt('VSubjet0') ? userInt('VSubjet0') : -1 ")
        ),
     cms.PSet(
        tag = cms.untracked.string("vSubjetIndex1"),
        quantity = cms.untracked.string("? hasUserInt('VSubjet1') ? userInt('VSubjet1') : -1 ")
        ),
#     cms.PSet(
#        tag = cms.untracked.string("topSubjetIndex0"),
#        quantity = cms.untracked.string("? hasUserInt('TopSubjet0') ? userInt('TopSubjet0') : -1 ")
#        ),
#     cms.PSet(
#        tag = cms.untracked.string("topSubjetIndex1"),
#        quantity = cms.untracked.string("? hasUserInt('TopSubjet1') ? userInt('TopSubjet1') : -1 ")
#        ),
#     cms.PSet(
#        tag = cms.untracked.string("topSubjetIndex2"),
#        quantity = cms.untracked.string("? hasUserInt('TopSubjet2') ? userInt('TopSubjet2') : -1 ")
#        ),
#     cms.PSet(
#        tag = cms.untracked.string("topSubjetIndex3"),
#        quantity = cms.untracked.string("? hasUserInt('TopSubjet3') ? userInt('TopSubjet3') : -1 ")
#        ),
     cms.PSet(
        tag = cms.untracked.string("tau1"),
        quantity = cms.untracked.string("userFloat('NjettinessAK8CHS:tau1')")
        ),
     cms.PSet(
        tag = cms.untracked.string("tau2"),
        quantity = cms.untracked.string("userFloat('NjettinessAK8CHS:tau2')")
        ),
     cms.PSet(
        tag = cms.untracked.string("tau3"),
        quantity = cms.untracked.string("userFloat('NjettinessAK8CHS:tau3')")
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
)


### jet variables
jetToolboxAK8SubjetVars = (
#### SUBSTRUCTURE
     cms.PSet(
        tag = cms.untracked.string("vSubjetIndex0"),
        quantity = cms.untracked.string("? hasUserInt('VSubjet0') ? userInt('VSubjet0') : -1 ")
        ),
     cms.PSet(
        tag = cms.untracked.string("vSubjetIndex1"),
        quantity = cms.untracked.string("? hasUserInt('VSubjet1') ? userInt('VSubjet1') : -1 ")
        ),
#     cms.PSet(
#        tag = cms.untracked.string("topSubjetIndex0"),
#        quantity = cms.untracked.string("? hasUserInt('TopSubjet0') ? userInt('TopSubjet0') : -1 ")
#        ),
#     cms.PSet(
#        tag = cms.untracked.string("topSubjetIndex1"),
#        quantity = cms.untracked.string("? hasUserInt('TopSubjet1') ? userInt('TopSubjet1') : -1 ")
#        ),
#     cms.PSet(
#        tag = cms.untracked.string("topSubjetIndex2"),
#        quantity = cms.untracked.string("? hasUserInt('TopSubjet2') ? userInt('TopSubjet2') : -1 ")
#        ),
#     cms.PSet(
#        tag = cms.untracked.string("topSubjetIndex3"),
#        quantity = cms.untracked.string("? hasUserInt('TopSubjet3') ? userInt('TopSubjet3') : -1 ")
#        ),
     cms.PSet(
        tag = cms.untracked.string("tau1"),
        quantity = cms.untracked.string("userFloat('NsubjettinessAK8PFCHSSoftDropSubjets:tau1')")
        ),
     cms.PSet(
        tag = cms.untracked.string("tau2"),
        quantity = cms.untracked.string("userFloat('NsubjettinessAK8PFCHSSoftDropSubjets:tau2')")
        ),
     cms.PSet(
        tag = cms.untracked.string("tau3"),
        quantity = cms.untracked.string("userFloat('NsubjettinessAK8PFCHSSoftDropSubjets:tau3')")
        )
)


### jet variables
jetToolboxAK8SubjetPuppiVars = (
#### SUBSTRUCTURE
     cms.PSet(
        tag = cms.untracked.string("vSubjetIndex0"),
        quantity = cms.untracked.string("? hasUserInt('VSubjet0') ? userInt('VSubjet0') : -1 ")
        ),
     cms.PSet(
        tag = cms.untracked.string("vSubjetIndex1"),
        quantity = cms.untracked.string("? hasUserInt('VSubjet1') ? userInt('VSubjet1') : -1 ")
        ),
#     cms.PSet(
#        tag = cms.untracked.string("topSubjetIndex0"),
#        quantity = cms.untracked.string("? hasUserInt('TopSubjet0') ? userInt('TopSubjet0') : -1 ")
#        ),
#     cms.PSet(
#        tag = cms.untracked.string("topSubjetIndex1"),
#        quantity = cms.untracked.string("? hasUserInt('TopSubjet1') ? userInt('TopSubjet1') : -1 ")
#        ),
#     cms.PSet(
#        tag = cms.untracked.string("topSubjetIndex2"),
#        quantity = cms.untracked.string("? hasUserInt('TopSubjet2') ? userInt('TopSubjet2') : -1 ")
#        ),
#     cms.PSet(
#        tag = cms.untracked.string("topSubjetIndex3"),
#        quantity = cms.untracked.string("? hasUserInt('TopSubjet3') ? userInt('TopSubjet3') : -1 ")
#        ),
     cms.PSet(
        tag = cms.untracked.string("tau1"),
        quantity = cms.untracked.string("userFloat('NsubjettinessAK8PFPuppiSoftDropSubjets:tau1')")
        ),
     cms.PSet(
        tag = cms.untracked.string("tau2"),
        quantity = cms.untracked.string("userFloat('NsubjettinessAK8PFPuppiSoftDropSubjets:tau2')")
        ),
     cms.PSet(
        tag = cms.untracked.string("tau3"),
        quantity = cms.untracked.string("userFloat('NsubjettinessAK8PFPuppiSoftDropSubjets:tau3')")
        )
)




jetToolboxAK8PuppiVars = (
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
        tag = cms.untracked.string("tau1"),
        quantity = cms.untracked.string("userFloat('NjettinessAK8Puppi:tau1')")
        ),
     cms.PSet(
        tag = cms.untracked.string("tau2"),
        quantity = cms.untracked.string("userFloat('NjettinessAK8Puppi:tau2')")
        ),
     cms.PSet(
        tag = cms.untracked.string("tau3"),
        quantity = cms.untracked.string("userFloat('NjettinessAK8Puppi:tau3')")
        ),
     cms.PSet(
        tag = cms.untracked.string("softDropMass"),
        quantity = cms.untracked.string("userFloat('ak8PFJetsPuppiSoftDropMass')")
        ),
     cms.PSet(
        tag = cms.untracked.string("trimmedMass"),
        quantity = cms.untracked.string("userFloat('ak8PFJetsPuppiTrimmedMass')")
        ),
     cms.PSet(
        tag = cms.untracked.string("prunedMass"),
        quantity = cms.untracked.string("userFloat('ak8PFJetsPuppiPrunedMass')")
        ),
     cms.PSet(
        tag = cms.untracked.string("filteredMass"),
        quantity = cms.untracked.string("userFloat('ak8PFJetsPuppiFilteredMass')")
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
      tag = cms.untracked.string("Iso03db"),
      quantity = cms.untracked.string("userFloat('iso03db')")
      ),
    cms.PSet(
      tag = cms.untracked.string("MiniIso"),
      quantity = cms.untracked.string("userFloat('miniIso')")
      ),
    cms.PSet(
      tag = cms.untracked.string("rho"),
      quantity = cms.untracked.string("userFloat('rho')")
      ),
    cms.PSet(
      tag = cms.untracked.string("EA"),
      quantity = cms.untracked.string("userFloat('EA')")
      ),
    cms.PSet(
      tag = cms.untracked.string("sumChargedHadronPt"),
      quantity = cms.untracked.string("userFloat('sumChargedHadronPt')")
      ),
    cms.PSet(
      tag = cms.untracked.string("sumNeutralHadronEt"),
      quantity = cms.untracked.string("userFloat('sumNeutralHadronEt')")
      ),
    cms.PSet(
      tag = cms.untracked.string("sumPhotonEt"),
      quantity = cms.untracked.string("userFloat('sumPhotonEt')")
      ),
    cms.PSet(
      tag = cms.untracked.string("sumPUPt"),
      quantity = cms.untracked.string("userFloat('sumPUPt')")
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
        quantity = cms.untracked.string("userFloat('missHits')")
        ),
    cms.PSet(
        tag = cms.untracked.string("hasMatchedConVeto"),
        quantity = cms.untracked.string("userFloat('hasMatchConv')")
        ),
    cms.PSet(
        tag = cms.untracked.string("vidVeto"),
        quantity = cms.untracked.string("userFloat('vidVeto')")
        ),
    cms.PSet(
        tag = cms.untracked.string("vidLoose"),
        quantity = cms.untracked.string("userFloat('vidLoose')")
        ),
    cms.PSet(
        tag = cms.untracked.string("vidTight"),
        quantity = cms.untracked.string("userFloat('vidTight')")
        ),
    cms.PSet(
        tag = cms.untracked.string("vidMedium"),
        quantity = cms.untracked.string("userFloat('vidMedium')")
        ),
    cms.PSet(
        tag = cms.untracked.string("vidHEEP"),
        quantity = cms.untracked.string("userFloat('vidHEEP')")
        ),
    cms.PSet(
        tag = cms.untracked.string("vidHEEPnoiso"),
        quantity = cms.untracked.string("userFloat('vidHEEPnoiso')")
        ),
    cms.PSet(
        tag = cms.untracked.string("SCEta"),
        quantity = cms.untracked.string("superCluster().eta()")
        ),
    cms.PSet(
        tag = cms.untracked.string("SCPhi"),
        quantity = cms.untracked.string("superCluster().phi()")
        )

    )


electrons  = copy.deepcopy(basic)
electrons.variables += electronVars
electrons.prefix = cms.untracked.string("el")
electrons.src = cms.InputTag("electronUserData")




###photons                                                           
photonVars = (
    cms.PSet(
        tag = cms.untracked.string("SCEta"),
        quantity = cms.untracked.string("superCluster().eta()")
        ),
    cms.PSet(
        tag = cms.untracked.string("SCPhi"),
        quantity = cms.untracked.string("superCluster.phi()")
        ),
    cms.PSet(
        tag = cms.untracked.string("HasPixelSeed"),
        quantity = cms.untracked.string("userInt('hasPixelSeed')")
        ),
    cms.PSet(
        tag = cms.untracked.string("SigmaIEtaIEta"),
        quantity = cms.untracked.string("userFloat('sigmaIetaIeta')")
        ),
    cms.PSet(
        tag = cms.untracked.string("HoverE"),
        quantity = cms.untracked.string("userFloat('hoe')")
        ),
    cms.PSet(
        tag = cms.untracked.string("R9"),
        quantity = cms.untracked.string("userFloat('r9')")
        ),
    cms.PSet(
        tag = cms.untracked.string("ChargedHadronIso"),
        quantity = cms.untracked.string("userFloat('isoC')")
        ),        
    cms.PSet(
        tag = cms.untracked.string("PhotonIso"),
        quantity = cms.untracked.string("userFloat('isoP')")
        ),
    cms.PSet(
        tag = cms.untracked.string("NeutralHadronIso"),
        quantity = cms.untracked.string("userFloat('isoN')")
        ),
    cms.PSet(
        tag = cms.untracked.string("ChargedHadronIsoEAcorrected"),
        quantity = cms.untracked.string("userFloat('isoC_EAcor')")
        ),        
    cms.PSet(
        tag = cms.untracked.string("PhotonIsoEAcorrected"),
        quantity = cms.untracked.string("userFloat('isoP_EAcor')")
        ),
    cms.PSet(
        tag = cms.untracked.string("NeutralHadronIsoEAcorrected"),
        quantity = cms.untracked.string("userFloat('isoN_EAcor')")
        ),
    cms.PSet(
        tag = cms.untracked.string("PassLooseID"),
        quantity = cms.untracked.string("userInt('isLoose')")
        ),
    cms.PSet(
        tag = cms.untracked.string("PassMediumID"),
        quantity = cms.untracked.string("userInt('isMedium')")
        ),
    cms.PSet(
        tag = cms.untracked.string("PassTightID"),
        quantity = cms.untracked.string("userInt('isTight')")
        )
    )
### copying the muon set of variables from basic,
### adding the set of variable which are related to muons only
photons  = copy.deepcopy(basic)
photons.variables += photonVars
photons.prefix = cms.untracked.string("pho")
photons.src = cms.InputTag("photonUserData")



###photonjets                                                                                                                              
photonjets =  cms.EDProducer(
    "CandViewNtpProducer",
    src=cms.InputTag("photonJets"),
    lazyParser=cms.untracked.bool(True),
    prefix=cms.untracked.string("basic"),
    eventInfo=cms.untracked.bool(False),
    variables = cms.VPSet(
    cms.PSet(
        tag = cms.untracked.string("JetIndex"),
        quantity = cms.untracked.string("userInt('jetIndex')")
      
        ),
    cms.PSet(
        tag = cms.untracked.string("PhotonIndex"),
        quantity = cms.untracked.string("userInt('phoIndex')")
            ),
    cms.PSet(
        tag = cms.untracked.string("SubwGammatIndex"),
        quantity = cms.untracked.string("userInt('subIndex')")
        ),
    cms.PSet(
        tag = cms.untracked.string("PhotonSubjetFrac"),
        quantity = cms.untracked.string("userFloat('phoSubjetPtFrac')")
        ),
    cms.PSet(
        tag = cms.untracked.string("SubjetPt0"),
        quantity = cms.untracked.string("userFloat('SubPt0')")
        ),
    cms.PSet(
        tag = cms.untracked.string("SubjetPt1"),
        quantity = cms.untracked.string("userFloat('SubPt1')")
        ),
    cms.PSet(
        tag = cms.untracked.string("SubjetPt2"),
        quantity = cms.untracked.string("userFloat('SubPt2')")
        ),
    cms.PSet(
        tag = cms.untracked.string("SubjetEta0"),
        quantity = cms.untracked.string("userFloat('SubEta0')")
        ),
    cms.PSet(
        tag = cms.untracked.string("SubjetEta1"),
        quantity = cms.untracked.string("userFloat('SubEta1')")
        ),
    cms.PSet(
        tag = cms.untracked.string("SubjetEta2"),
        quantity = cms.untracked.string("userFloat('SubEta2')")
        ),
    cms.PSet(
        tag = cms.untracked.string("SubjetPhi0"),
        quantity = cms.untracked.string("userFloat('SubPhi0')")
        ),
    cms.PSet(
        tag = cms.untracked.string("SubjetPhi1"),
        quantity = cms.untracked.string("userFloat('SubPhi1')")
        ),
    cms.PSet(
        tag = cms.untracked.string("SubjetPhi2"),
        quantity = cms.untracked.string("userFloat('SubPhi2')")
        ),
    cms.PSet(
        tag = cms.untracked.string("SubjetEne0"),
        quantity = cms.untracked.string("userFloat('SubEne0')")
        ),
    cms.PSet(
        tag = cms.untracked.string("SubjetEne1"),
        quantity = cms.untracked.string("userFloat('SubEne1')")
        ),
    cms.PSet(
        tag = cms.untracked.string("SubjetEne2"),
        quantity = cms.untracked.string("userFloat('SubEne2')")
        ),

    )
)
#photonjets = copy.deepcopy(PhotonVars)                                                                                                    
#photonjets.variables += photonJetVars                                                                                                     
#photonjets = photonJetVars                                                                                                                
#photonjets.prefix = cms.untracked.string("phoJets")


qglVars = (
    ### B-TAGGING
    cms.PSet(
      tag = cms.untracked.string("QGL"),
      #quantity = cms.untracked.string("userFloat('QGL')")
      quantity = cms.untracked.string("userFloat('QGTaggerAK4PFCHS:qgLikelihood')")
      ),
)

###jets
jetsAK4CHS = copy.deepcopy(basic)
jetsAK4CHS.variables += jetVars
jetsAK4CHS.variables += qglVars
jetsAK4CHS.variables += jetVarsForSys
jetsAK4CHS.prefix = cms.untracked.string("jetAK4CHS")
jetsAK4CHS.src = cms.InputTag("jetUserData")
jetKeysAK4CHS = copy.deepcopy( jetKeys )
jetKeysAK4CHS.jetLabel = cms.InputTag("jetUserData")


###jets with Puppi
jetsAK4Puppi = copy.deepcopy(basic)
jetsAK4Puppi.variables += jetVars
jetsAK4Puppi.variables += jetVarsForSys
jetsAK4Puppi.prefix = cms.untracked.string("jetAK4Puppi")
jetsAK4Puppi.src = cms.InputTag("jetUserDataPuppi")
jetKeysAK4Puppi = copy.deepcopy( jetKeys )
jetKeysAK4Puppi.jetLabel = cms.InputTag("jetUserDataPuppi")

###AK8 jets with CHS
jetsAK8CHS = copy.deepcopy(basic)
jetsAK8CHS.variables += jetVars
jetsAK8CHS.variables += jetVarsForSys
jetsAK8CHS.variables += jetToolboxAK8Vars
jetsAK8CHS.prefix = cms.untracked.string("jetAK8CHS")
jetsAK8CHS.src = cms.InputTag("boostedJetUserDataAK8")
jetKeysAK8CHS = copy.deepcopy( jetKeys )
jetKeysAK8CHS.jetLabel = cms.InputTag("boostedJetUserDataAK8")

### AK8 jets with Puppi
jetsAK8Puppi = copy.deepcopy(basic)
jetsAK8Puppi.variables += jetVars
jetsAK8Puppi.variables += jetVarsForSys
jetsAK8Puppi.variables += jetToolboxAK8PuppiVars 
jetsAK8Puppi.prefix = 'jetAK8Puppi'
jetsAK8Puppi.src = cms.InputTag( 'boostedJetUserDataAK8Puppi' )
jetKeysAK8Puppi = jetKeysAK8CHS.clone( jetLabel = 'boostedJetUserDataAK8Puppi' )


###subjetsAK8 with CHS		
subjetsAK8CHS = copy.deepcopy(basic)		
subjetsAK8CHS.variables += jetVars		
subjetsAK8CHS.variables += jetToolboxAK8SubjetVars		
subjetsAK8CHS.prefix = cms.untracked.string("subjetAK8CHS")		
subjetsAK8CHS.src = cms.InputTag("selectedPatJetsAK8PFCHSSoftDropPacked", "SubJets")		
subjetKeysAK8CHS = copy.deepcopy( jetKeys )		
subjetKeysAK8CHS.jetLabel = cms.InputTag("selectedPatJetsAK8PFCHSSoftDropPacked", "SubJets")

###subjetsAK8Puppi
subjetsAK8Puppi = copy.deepcopy(basic)
subjetsAK8Puppi.variables += jetVars
subjetsAK8Puppi.variables += jetToolboxAK8SubjetPuppiVars
subjetsAK8Puppi.prefix = cms.untracked.string("subjetAK8Puppi")
subjetsAK8Puppi.src = cms.InputTag("selectedPatJetsAK8PFPuppiSoftDropPacked", "SubJets")
subjetKeysAK8Puppi = copy.deepcopy( jetKeys )
subjetKeysAK8Puppi.jetLabel = cms.InputTag('selectedPatJetsAK8PFPuppiSoftDropPacked', "SubJets")


###genPart
genPart = copy.deepcopy(basic)
genPart.variables += genPartVars
genPart.prefix = cms.untracked.string("genPart")
genPart.src = cms.InputTag("prunedGenParticles")

###genJetsAK8
genJetsAK8 = copy.deepcopy(basic)
genJetsAK8.prefix = cms.untracked.string("genJetsAK8")
genJetsAK8.src = cms.InputTag("ak8GenJetsNoNu")

###genJetsAK8 soft drop
genJetsAK8SoftDrop = copy.deepcopy(basic)
genJetsAK8SoftDrop.prefix = cms.untracked.string("genJetsAK8SoftDrop")
genJetsAK8SoftDrop.src = cms.InputTag("ak8GenJetsNoNuSoftDrop")

###genJetsAK8 subjets
genJetsAK8SoftDropSubjets = copy.deepcopy(basic)
genJetsAK8SoftDropSubjets.prefix = cms.untracked.string("genJetsAK8SoftDropSubJets")
genJetsAK8SoftDropSubjets.src = cms.InputTag("ak8GenJetsNoNuSoftDrop", "SubJets")


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


### No HF MET used as default for the time being
'''
metNoHF = copy.deepcopy(metFull)
metNoHF.prefix = cms.untracked.string("metNoHF")
metNoHF.src = cms.InputTag("skimmedPatMETNoHF")
'''

print "DONE STANDARD"
