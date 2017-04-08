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
    #Calc cms.PSet(
    #Calc tag = cms.untracked.string("Mass"),
    #Calc quantity = cms.untracked.string("mass")
    #Calc ),
    cms.PSet(
    tag = cms.untracked.string("Pt"),
    quantity = cms.untracked.string("pt")
    ),
    cms.PSet(
    tag = cms.untracked.string("Eta"),
    quantity = cms.untracked.string("eta")
    ),
    #Calc cms.PSet(
    #Calc tag = cms.untracked.string("Y"),
    #Calc quantity = cms.untracked.string("rapidity")
    #Calc ),
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

puppimetFull = metFull.clone(
  src = cms.InputTag("skimmedPatPuppiMET"),
  prefix = cms.untracked.string("puppimetFull"),
  )

### muon variables
muonVars = (
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
   # Impact point
   cms.PSet(
        tag = cms.untracked.string("Dxy"),
        quantity = cms.untracked.string("userFloat('dxy')")
   ),
   # cms.PSet(
   #      tag = cms.untracked.string("Dxyerr"),
   #      quantity = cms.untracked.string("userFloat('dxyErr')")
   # ),
   cms.PSet(
        tag = cms.untracked.string("Dz"),
        quantity = cms.untracked.string("userFloat('dz')")
        ),
   # cms.PSet(
   #   tag = cms.untracked.string("Dzerr"),
   #   quantity = cms.untracked.string("userFloat('dzErr')")
   #   ),
   cms.PSet(
        tag = cms.untracked.string("DB"),
        quantity = cms.untracked.string("userFloat('dB')")
   ),
   cms.PSet(
        tag = cms.untracked.string("DBerr"),
        quantity = cms.untracked.string("userFloat('dBErr')")
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
        tag = cms.untracked.string("IsMediumMuon2016"),
        quantity = cms.untracked.string("userFloat('isMediumMuon2016')")
        ),
    cms.PSet(
        tag = cms.untracked.string("IsTightMuon"),
        quantity = cms.untracked.string("userFloat('isTightMuon')")
        ),
    cms.PSet(
        tag = cms.untracked.string("IsHighPtMuon"),
        quantity = cms.untracked.string("userFloat('isHighPtMuon')")
        ),
    ### High pT muon variables from https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideMuonIdRun2#High_pT_Muon_pT_assignment_detai
   cms.PSet(
       tag = cms.untracked.string("InnerTrackPt"),
       quantity = cms.untracked.string("? innerTrack.isNonnull ? innerTrack.pt : -900")
       ),
   cms.PSet(
       tag = cms.untracked.string("TunePMuonBestTrackPt"),
       quantity = cms.untracked.string("? tunePMuonBestTrack.isNonnull ? tunePMuonBestTrack.pt : -900")
       ),
   ### variables used in ID
   ### https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId#Tight_Muon_selection
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
   ### MEDIUM2016
   cms.PSet(
       tag = cms.untracked.string("CombQualChi2LocalPos"),
       quantity = cms.untracked.string("combinedQuality.chi2LocalPosition")
       ),
   cms.PSet(
       tag = cms.untracked.string("CombQualTrkKink"),
       quantity = cms.untracked.string("combinedQuality.trkKink")
       ),
   cms.PSet(
       tag = cms.untracked.string("InTrkValidFraction"),
       quantity = cms.untracked.string("? innerTrack.isNonnull ? innerTrack.validFraction : -900")
       ),
   cms.PSet(
       tag = cms.untracked.string("SegmentCompatibility"),
       quantity = cms.untracked.string("userFloat('segmentCompatibility')")
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
   ### variables used in isolation
   ### https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId#Accessing_PF_Isolation_from_reco
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
   cms.PSet(
       tag = cms.untracked.string("TrackerSumPt"),
       quantity = cms.untracked.string("isolationR03().sumPt")
       ),
   ### genLepton
   #Calc cms.PSet(
   #Calc     tag = cms.untracked.string("GenMuonY"),
   #Calc     quantity = cms.untracked.string("? genParticleRef.isNonnull ? genLepton.rapidity : -900")
   #Calc     ),
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
    #cms.PSet(
    #  tag = cms.untracked.string("CMVA"),
    #  quantity = cms.untracked.string("bDiscriminator('pfCombinedMVAV2BJetTags')")
    #  ),
    ### GEN PARTON
    #Calc cms.PSet(
    #Calc   tag = cms.untracked.string("GenPartonY"),
    #Calc   quantity = cms.untracked.string("? genParticleRef.isNonnull ? genParton.rapidity : -900")
    #Calc   ),
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
    #Calc cms.PSet(
    #Calc   tag = cms.untracked.string("GenJetY"),
    #Calc   quantity = cms.untracked.string("? genJetFwdRef.isNonnull ? genJet.rapidity : -900")
    #Calc   ),
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
    #NonID  cms.PSet(
    #NonID      tag = cms.untracked.string("HFHadronEnergy"),
    #NonID      quantity = cms.untracked.string("? isPFJet ? HFHadronEnergy : -1")
    #NonID      ),
    #NonID  cms.PSet(
    #NonID      tag = cms.untracked.string("HFEMEnergy"),
    #NonID      quantity = cms.untracked.string("? isPFJet ? HFEMEnergy : -1")
    #NonID      ),
    #NonID  cms.PSet(
    #NonID      tag = cms.untracked.string("ChargedHadronMultiplicity"),
    #NonID      quantity = cms.untracked.string("? isPFJet ? chargedHadronMultiplicity : -1")
    #NonID      ),
    #NonID  cms.PSet(
    #NonID      tag = cms.untracked.string("numberOfDaughters"),
    #NonID      quantity = cms.untracked.string("? isPFJet ? numberOfDaughters : -1")
    #NonID      ),
    #NonID  cms.PSet(
    #NonID      tag = cms.untracked.string("neutralHadronMultiplicity"),
    #NonID      quantity = cms.untracked.string("? isPFJet ? neutralHadronMultiplicity : -1")
    #NonID      ),
    #Calc cms.PSet(
    #Calc     tag = cms.untracked.string("neutralHadronEnergy"),
    #Calc     quantity = cms.untracked.string("? isPFJet ? neutralHadronEnergy : -1")
    #Calc     ),
    cms.PSet(
        tag = cms.untracked.string("neutralEmEnergy"),
        quantity = cms.untracked.string("? isPFJet ? neutralEmEnergy : -1"),
        ),
    #Calc cms.PSet(
    #Calc     tag = cms.untracked.string("chargedEmEnergy"),
    #Calc     quantity = cms.untracked.string("? isPFJet ? chargedEmEnergy : -1"),
    #Calc     ),
    #Calc cms.PSet(
    #Calc     tag = cms.untracked.string("chargedHadronEnergy"),
    #Calc     quantity = cms.untracked.string("? isPFJet ? chargedHadronEnergy : -1"),
    #Calc     ),
    cms.PSet(
        tag = cms.untracked.string("photonMultiplicity"),
        quantity = cms.untracked.string("? isPFJet ? photonMultiplicity : -1")
        ),
    cms.PSet(
        tag = cms.untracked.string("electronMultiplicity"),
        quantity = cms.untracked.string("? isPFJet ? electronMultiplicity : -1")
        ),
    #NonID  cms.PSet(
    #NonID      tag = cms.untracked.string("HFHadronMultiplicity"),
    #NonID      quantity = cms.untracked.string("? isPFJet ? HFHadronMultiplicity : -1")
    #NonID      ),
    #NonID  cms.PSet(
    #NonID      tag = cms.untracked.string("HFEMMultiplicity"),
    #NonID      quantity = cms.untracked.string("? isPFJet ? HFEMMultiplicity : -1")
    #NonID      ),
    cms.PSet(
        tag = cms.untracked.string("ChargedMuEnergy"),
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
    #Calc cms.PSet(
    #Calc     tag = cms.untracked.string("muonEnergyFrac"),
    #Calc     quantity = cms.untracked.string("? isPFJet ?muonEnergyFraction : -1")
    #Calc     ),
    cms.PSet(
        tag = cms.untracked.string("chargedEmEnergyFrac"),
        quantity = cms.untracked.string("? isPFJet ?chargedEmEnergyFraction : -1")
        ),
    cms.PSet(
        tag = cms.untracked.string("chargedMultiplicity"),
        quantity = cms.untracked.string("? isPFJet ?chargedMultiplicity : -1")
        ),
    #Calc cms.PSet(
    #Calc     tag = cms.untracked.string("NumConstituents"),
    #Calc     quantity = cms.untracked.string("? isPFJet ? chargedMultiplicity + neutralMultiplicity : -1")
    #Calc      ),

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
    #Calc cms.PSet(
    #Calc     tag = cms.untracked.string("SmearedE"),
    #Calc     quantity = cms.untracked.string("userFloat('SmearedE')")
    #Calc     ),
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
    cms.PSet(
      tag = cms.untracked.string("DoubleBAK8"),
      quantity = cms.untracked.string("bDiscriminator('pfBoostedDoubleSecondaryVertexAK8BJetTags')")
      ),
    cms.PSet(
      tag = cms.untracked.string("DoubleBCA15"),
      quantity = cms.untracked.string("bDiscriminator('pfBoostedDoubleSecondaryVertexCA15BJetTags')")
      ),
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
        tag = cms.untracked.string("vSubjetPuppiIndex0"),
        quantity = cms.untracked.string("? hasUserInt('VSubjetPuppi0') ? userInt('VSubjetPuppi0') : -1 ")
        ),
     cms.PSet(
        tag = cms.untracked.string("vSubjetPuppiIndex1"),
        quantity = cms.untracked.string("? hasUserInt('VSubjetPuppi1') ? userInt('VSubjetPuppi1') : -1 ")
        ),
     cms.PSet(
        tag = cms.untracked.string("tau1CHS"),
        quantity = cms.untracked.string("userFloat('NjettinessAK8CHS:tau1')")
        ),
     cms.PSet(
        tag = cms.untracked.string("tau2CHS"),
        quantity = cms.untracked.string("userFloat('NjettinessAK8CHS:tau2')")
        ),
     cms.PSet(
        tag = cms.untracked.string("tau3CHS"),
        quantity = cms.untracked.string("userFloat('NjettinessAK8CHS:tau3')")
        ),
     cms.PSet(
        tag = cms.untracked.string("softDropMassCHS"),
        quantity = cms.untracked.string("userFloat('ak8PFJetsCHSSoftDropMass')")
        ),
     cms.PSet(
        tag = cms.untracked.string("trimmedMassCHS"),
        quantity = cms.untracked.string("userFloat('ak8PFJetsCHSTrimmedMass')")
        ),
     cms.PSet(
        tag = cms.untracked.string("prunedMassCHS"),
        quantity = cms.untracked.string("userFloat('ak8PFJetsCHSPrunedMass')")
        ),
     cms.PSet(
        tag = cms.untracked.string("filteredMassCHS"),
        quantity = cms.untracked.string("userFloat('ak8PFJetsCHSFilteredMass')")
        ),
     cms.PSet(
        tag = cms.untracked.string("softDropMassPuppi"),
        quantity = cms.untracked.string("userFloat('ak8PFJetsPuppiValueMap:softDropMassPuppi')")
        ),
     cms.PSet(
        tag = cms.untracked.string("PtPuppi"),
        quantity = cms.untracked.string("userFloat('ak8PFJetsPuppiValueMap:pt')")
        ),
     cms.PSet(
        tag = cms.untracked.string("EtaPuppi"),
        quantity = cms.untracked.string("userFloat('ak8PFJetsPuppiValueMap:eta')")
        ),
     cms.PSet(
        tag = cms.untracked.string("PhiPuppi"),
        quantity = cms.untracked.string("userFloat('ak8PFJetsPuppiValueMap:phi')")
        ),
     cms.PSet(
        tag = cms.untracked.string("MassPuppi"),
        quantity = cms.untracked.string("userFloat('ak8PFJetsPuppiValueMap:mass')")
        ),
     cms.PSet(
        tag = cms.untracked.string("tau1Puppi"),
        quantity = cms.untracked.string("userFloat('ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau1')")
        ),
     cms.PSet(
        tag = cms.untracked.string("tau2Puppi"),
        quantity = cms.untracked.string("userFloat('ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau2')")
        ),
     cms.PSet(
        tag = cms.untracked.string("tau3Puppi"),
        quantity = cms.untracked.string("userFloat('ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau3')")
        ),
     cms.PSet(
        tag = cms.untracked.string("uncorrSDMassAK8Puppi"),
        quantity = cms.untracked.string("? hasUserFloat('subjetSumMassSoftDropPuppi') ? userFloat('subjetSumMassSoftDropPuppi') : -999 ")
        ),
)

### jet variables
jetToolboxAK8SubjetVars = (
#### SUBSTRUCTURE
#     cms.PSet(
#        tag = cms.untracked.string("vSubjetIndex0"),
#        quantity = cms.untracked.string("? hasUserInt('VSubjet0') ? userInt('VSubjet0') : -1 ")
#        ),
#     cms.PSet(
#        tag = cms.untracked.string("vSubjetIndex1"),
#        quantity = cms.untracked.string("? hasUserInt('VSubjet1') ? userInt('VSubjet1') : -1 ")
#        ),
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
        ),
)


### jet variables
jetToolboxAK8SubjetPuppiVars = (
#### SUBSTRUCTURE
#     cms.PSet(
#        tag = cms.untracked.string("vSubjetIndex0"),
#        quantity = cms.untracked.string("? hasUserInt('VSubjet0') ? userInt('VSubjet0') : -1 ")
#        ),
#     cms.PSet(
#        tag = cms.untracked.string("vSubjetIndex1"),
#        quantity = cms.untracked.string("? hasUserInt('VSubjet1') ? userInt('VSubjet1') : -1 ")
#        ),
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
    cms.PSet(
      tag = cms.untracked.string("DoubleBAK8"),
      quantity = cms.untracked.string("bDiscriminator('pfBoostedDoubleSecondaryVertexAK8BJetTags')")
      ),
    cms.PSet(
      tag = cms.untracked.string("DoubleBCA15"),
      quantity = cms.untracked.string("bDiscriminator('pfBoostedDoubleSecondaryVertexCA15BJetTags')")
      ),
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
        tag = cms.untracked.string("vSubjetPuppiIndex0"),
        quantity = cms.untracked.string("? hasUserInt('VSubjetPuppi0') ? userInt('VSubjetPuppi0') : -1 ")
        ),
     cms.PSet(
        tag = cms.untracked.string("vSubjetPuppiIndex1"),
        quantity = cms.untracked.string("? hasUserInt('VSubjetPuppi1') ? userInt('VSubjetPuppi1') : -1 ")
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
    # For Isolation
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
    # Impact point
    cms.PSet(
      tag = cms.untracked.string("Dxy"),
      quantity = cms.untracked.string("userFloat('dxy')")
      ),
    cms.PSet(
      tag = cms.untracked.string("Dz"),
      quantity = cms.untracked.string("userFloat('dz')")
      ),
    cms.PSet(
      tag = cms.untracked.string("DB"),
      quantity = cms.untracked.string("userFloat('dB')")
      ),
    cms.PSet(
      tag = cms.untracked.string("DBerr"),
      quantity = cms.untracked.string("userFloat('dBErr')")
      ),
    # Cut-based ID variables
    cms.PSet(
      tag = cms.untracked.string("dEtaIn"),
      quantity = cms.untracked.string("deltaEtaSuperClusterTrackAtVtx")
      ),
    cms.PSet(
      tag = cms.untracked.string("dEtaInSeed"),
      quantity = cms.untracked.string("userFloat('dEtaInSeed')")
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
        tag = cms.untracked.string("SCEta"),
        quantity = cms.untracked.string("superCluster().eta()")
        ),
    cms.PSet(
        tag = cms.untracked.string("SCPhi"),
        quantity = cms.untracked.string("superCluster().phi()")
        ),
    # IDs
    cms.PSet(
        tag = cms.untracked.string("vidVeto"),
        quantity = cms.untracked.string("userFloat('vidVeto')")
        ),
    cms.PSet(
        tag = cms.untracked.string("vidLoose"),
        quantity = cms.untracked.string("userFloat('vidLoose')")
        ),
    cms.PSet(
        tag = cms.untracked.string("vidMedium"),
        quantity = cms.untracked.string("userFloat('vidMedium')")
        ),
    cms.PSet(
        tag = cms.untracked.string("vidTight"),
        quantity = cms.untracked.string("userFloat('vidTight')")
        ),
    cms.PSet(
        tag = cms.untracked.string("vidHEEP"),
        quantity = cms.untracked.string("userFloat('vidHEEP')")
        ),
    # IDs sans iso
    cms.PSet(
        tag = cms.untracked.string("vidVetonoiso"),
        quantity = cms.untracked.string("userFloat('vidVetonoiso')")
        ),
    cms.PSet(
        tag = cms.untracked.string("vidLoosenoiso"),
        quantity = cms.untracked.string("userFloat('vidLoosenoiso')")
        ),
    cms.PSet(
        tag = cms.untracked.string("vidMediumnoiso"),
        quantity = cms.untracked.string("userFloat('vidMediumnoiso')")
        ),
    cms.PSet(
        tag = cms.untracked.string("vidTightnoiso"),
        quantity = cms.untracked.string("userFloat('vidTightnoiso')")
        ),
    cms.PSet(
        tag = cms.untracked.string("vidHEEPnoiso"),
        quantity = cms.untracked.string("userFloat('vidHEEPnoiso')")
        ),
    cms.PSet(
        tag = cms.untracked.string("vidMvaGPvalue"),
        quantity = cms.untracked.string("userFloat('vidMvaGPvalue')")
        ),
    cms.PSet(
        tag = cms.untracked.string("vidMvaGPcateg"),
        quantity = cms.untracked.string("userInt('vidMvaGPcateg')")
        ),
    cms.PSet(
        tag = cms.untracked.string("vidMvaHZZvalue"),
        quantity = cms.untracked.string("userFloat('vidMvaHZZvalue')")
        ),
    cms.PSet(
        tag = cms.untracked.string("vidMvaHZZcateg"),
        quantity = cms.untracked.string("userInt('vidMvaHZZcateg')")
        ),
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
        tag = cms.untracked.string("SCRawE"),
        quantity = cms.untracked.string("superCluster.rawEnergy()")
        ),
    cms.PSet(
        tag = cms.untracked.string("HasPixelSeed"),
        quantity = cms.untracked.string("userInt('hasPixelSeed')")
        ),
    cms.PSet(
        tag = cms.untracked.string("ElectronVeto"),
        quantity = cms.untracked.string("userInt('eleveto')")
        ),
    cms.PSet(
        tag = cms.untracked.string("SigmaIEtaIEta"),
        quantity = cms.untracked.string("userFloat('sigmaIetaIeta')")
        ),
    cms.PSet(
        tag = cms.untracked.string("SigmaIEtaIPhi"),
        quantity = cms.untracked.string("userFloat('sigmaIetaIphi')")
        ),
    cms.PSet(
        tag = cms.untracked.string("SigmaIPhiIPhi"),
        quantity = cms.untracked.string("userFloat('sigmaIphiIphi')")
        ),
    cms.PSet(
        tag = cms.untracked.string("E1x5"),
        quantity = cms.untracked.string("userFloat('e1x5')")
        ),
    cms.PSet(
        tag = cms.untracked.string("E5x5"),
        quantity = cms.untracked.string("userFloat('e5x5')")
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
        tag = cms.untracked.string("PhotonIsoEAcorrectedsp15"),
        quantity = cms.untracked.string("userFloat('isoP_EAcor')")
        ),
    cms.PSet(
        tag = cms.untracked.string("NeutralHadronIsoEAcorrectedsp15"),
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
photonjets =  (
    cms.PSet(
        tag = cms.untracked.string("PhotonIndex"),
        quantity = cms.untracked.string("? hasUserInt('phoIndex') ? userInt('phoIndex') : -999 ")
            ),
    cms.PSet(
        tag = cms.untracked.string("SubwGammatIndex"),
        quantity = cms.untracked.string("? hasUserInt('subIndex') ? userInt('subIndex') : -999 ")
        ),
    cms.PSet(
        tag = cms.untracked.string("PhotonSubjetFrac"),
        quantity = cms.untracked.string("? hasUserFloat('phoSubjetPtFrac') ? userFloat('phoSubjetPtFrac') : -999 ")
        ),
    cms.PSet(
        tag = cms.untracked.string("SubjetPt0"),
        quantity = cms.untracked.string("? hasUserFloat('SubPt0') ? userFloat('SubPt0') : -999 ")
        ),
    cms.PSet(
        tag = cms.untracked.string("SubjetPt1"),
        quantity = cms.untracked.string("? hasUserFloat('SubPt1') ? userFloat('SubPt1') : -999 ")
        ),
    cms.PSet(
        tag = cms.untracked.string("SubjetPt2"),
        quantity = cms.untracked.string("? hasUserFloat('SubPt2') ? userFloat('SubPt2') : -999 ")
        ),
    cms.PSet(
        tag = cms.untracked.string("SubjetEta0"),
        quantity = cms.untracked.string("? hasUserFloat('SubEta0') ? userFloat('SubEta0') : -999 ")
        ),
    cms.PSet(
        tag = cms.untracked.string("SubjetEta1"),
        quantity = cms.untracked.string("? hasUserFloat('SubEta1') ? userFloat('SubEta1') : -999 ")
        ),
    cms.PSet(
        tag = cms.untracked.string("SubjetEta2"),
        quantity = cms.untracked.string("? hasUserFloat('SubEta2') ? userFloat('SubEta2') : -999 ")
        ),
    cms.PSet(
        tag = cms.untracked.string("SubjetPhi0"),
        quantity = cms.untracked.string("? hasUserFloat('SubPhi0') ? userFloat('SubPhi0') : -999 ")
        ),
    cms.PSet(
        tag = cms.untracked.string("SubjetPhi1"),
        quantity = cms.untracked.string("? hasUserFloat('SubPhi1') ? userFloat('SubPhi1') : -999 ")
        ),
    cms.PSet(
        tag = cms.untracked.string("SubjetPhi2"),
        quantity = cms.untracked.string("? hasUserFloat('SubPhi2') ? userFloat('SubPhi2') : -999 ")
        ),
    cms.PSet(
        tag = cms.untracked.string("SubjetEne0"),
        quantity = cms.untracked.string("? hasUserFloat('SubEne0') ? userFloat('SubEne0') : -999 ")
        ),
    cms.PSet(
        tag = cms.untracked.string("SubjetEne1"),
        quantity = cms.untracked.string("? hasUserFloat('SubEne1') ? userFloat('SubEne1') : -999 ")
        ),
    cms.PSet(
        tag = cms.untracked.string("SubjetEne2"),
        quantity = cms.untracked.string("? hasUserFloat('SubEne2') ? userFloat('SubEne2') : -999 ")
        ),

)


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
jetsAK8CHS.variables += photonjets
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
#subjetsAK8CHS.src = cms.InputTag("slimmedJetsAK8PFCHSSoftDropPacked", "SubJets")		
subjetKeysAK8CHS = copy.deepcopy( jetKeys )		
subjetKeysAK8CHS.jetLabel = cms.InputTag("selectedPatJetsAK8PFCHSSoftDropPacked", "SubJets")
#subjetKeysAK8CHS.jetLabel = cms.InputTag("slimmedJetsAK8PFCHSSoftDropPacked", "SubJets")

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
genJetsAK8.src = cms.InputTag("slimmedGenJetsAK8")

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
