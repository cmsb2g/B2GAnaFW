#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDMException.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

// Fastjet (for creating subjets)
#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>
#include "fastjet/tools/Filter.hh"
#include <fastjet/ClusterSequence.hh>
#include <fastjet/ClusterSequenceArea.hh>

// N-subjettiness algos
#include "B2GAnaFW/B2GAnaFW/interface/Nsubjettiness.hh"
#include "B2GAnaFW/B2GAnaFW/interface/Njettiness.hh"

// Vertex
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

// trigger
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h" // gives access to the (release cycle dependent) trigger object codes

#include <TFile.h>
#include <TH1F.h>
#include <TGraphAsymmErrors.h>
#include <TLorentzVector.h>
#include <vector>

using namespace fastjet;
using namespace reco;
using namespace edm;
using namespace std;
using namespace trigger;

typedef std::vector<pat::Jet> PatJetCollection;

struct orderByPt {
  const std::string mCorrLevel;
  orderByPt(const std::string& fCorrLevel) : mCorrLevel(fCorrLevel) {}
  bool operator ()(PatJetCollection::const_iterator const& a, PatJetCollection::const_iterator const& b) {
    if( mCorrLevel=="Uncorrected" ) return a->correctedJet("Uncorrected").pt() > b->correctedJet("Uncorrected").pt();
    else return a->pt() > b->pt();
  }
};

class JetUserData : public edm::EDProducer {
  public:
    JetUserData( const edm::ParameterSet & );   

  private:
    void produce( edm::Event &, const edm::EventSetup & );
    bool isMatchedWithTrigger(const pat::Jet&, trigger::TriggerObjectCollection,int&,double&,double);
    double getResolutionRatio(double eta);
    double getJERup(double eta);
    double getJERdown(double eta);

    void matchPackedJets(const edm::Handle<PatJetCollection>& jets,
        const edm::Handle<PatJetCollection>& packedjets,
        std::vector<int>& matchedIndices) ;

    edm::EDGetTokenT<std::vector<pat::Jet> >     jetToken_;
    edm::EDGetTokenT<std::vector<reco::Vertex> > pvToken_;

    //InputTag jetLabel_;
    InputTag pvLabel_;
    InputTag jLabel_, packedjLabel_, sjLabel_ ; 
    InputTag elLabel_, muLabel_;
    InputTag triggerResultsLabel_, triggerSummaryLabel_;
    InputTag hltJetFilterLabel_;
    std::string hltPath_;
    double hlt2reco_deltaRmax_;
    HLTConfigProvider hltConfig;
    int triggerBit;
    bool doSubjets_ ;
};


JetUserData::JetUserData(const edm::ParameterSet& iConfig) :

  pvLabel_            (iConfig.getParameter<edm::InputTag>("pv")),   // "offlinePrimaryVertex"
  jLabel_             (iConfig.getParameter<edm::InputTag>("jetLabel")),
  packedjLabel_       (iConfig.getParameter<edm::InputTag>("packedjetLabel")),
  sjLabel_            (iConfig.getParameter<edm::InputTag>("subjetLabel")),
  elLabel_            (iConfig.getParameter<edm::InputTag>("elLabel")),
  muLabel_            (iConfig.getParameter<edm::InputTag>("muLabel")),

  triggerResultsLabel_(iConfig.getParameter<edm::InputTag>("triggerResults")),
  triggerSummaryLabel_(iConfig.getParameter<edm::InputTag>("triggerSummary")),
  hltJetFilterLabel_  (iConfig.getParameter<edm::InputTag>("hltJetFilter")),   //trigger objects we want to match
  hltPath_            (iConfig.getParameter<std::string>("hltPath")),
  hlt2reco_deltaRmax_ (iConfig.getParameter<double>("hlt2reco_deltaRmax")),
  doSubjets_          (iConfig.getParameter<bool>("doSubjets"))
{
  produces<vector<pat::Jet> >();
}

#include "DataFormats/JetReco/interface/Jet.h"

void JetUserData::produce( edm::Event& iEvent, const edm::EventSetup& iSetup) {

  bool isMC = (!iEvent.isRealData());

  edm::Handle<std::vector<reco::Vertex> > pvHandle;
  iEvent.getByLabel(pvLabel_, pvHandle);

  edm::Handle<std::vector<pat::Jet> > jetHandle, packedjetHandle, subjetHandle;
  iEvent.getByLabel(jLabel_, jetHandle);
  auto_ptr<vector<pat::Jet> > jetColl( new vector<pat::Jet> (*jetHandle) );

  edm::Handle<std::vector<pat::Electron> > elHandle ; 
  iEvent.getByLabel(elLabel_, elHandle);

  edm::Handle<std::vector<pat::Muon> > muHandle ; 
  iEvent.getByLabel(muLabel_, muHandle);

  //// TRIGGER (this is not really needed ...)
  bool changedConfig = false;
  bool pathFound = false;
  if (!hltConfig.init(iEvent.getRun(), iSetup, "HLT", changedConfig)) {
    edm::LogError("HLTConfigProvider") << "Initialization of HLTConfigProvider failed!!" << std::endl;
    return;
  }

  if (changedConfig){
    edm::LogInfo("HLTMenu") << "the current menu is " << hltConfig.tableName() << std::endl;
    triggerBit = -1;
    for (size_t j = 0; j < hltConfig.triggerNames().size(); j++) {
      if (TString(hltConfig.triggerNames()[j]).Contains(hltPath_)) {triggerBit = j;pathFound=true;}
    }
    if (triggerBit == -1) edm::LogError("NoHLTPath") << "HLT path not found" << std::endl;
  }

     edm::Handle<edm::TriggerResults> triggerResults;
     iEvent.getByLabel(triggerResultsLabel_, triggerResults);
  /* Why do we need this
     if (size_t(triggerBit) < triggerResults->size() && pathFound)
       if (triggerResults->accept(triggerBit))
         std::cout << "event pass : " << hltPath_ << std::endl;
  */ 

  //// TRIGGER MATCHING
  trigger::TriggerObjectCollection JetLegObjects;

  edm::Handle<trigger::TriggerEvent> triggerSummary;

  if ( triggerSummary.isValid() ) {
    iEvent.getByLabel(triggerSummaryLabel_, triggerSummary);

    // Results from TriggerEvent product - Attention: must look only for
    // modules actually run in this path for this event!
    if(pathFound){
      const unsigned int triggerIndex(hltConfig.triggerIndex(hltPath_));
      const vector<string>& moduleLabels(hltConfig.moduleLabels(triggerIndex));
      const unsigned int moduleIndex(triggerResults->index(triggerIndex));
      for (unsigned int j=0; j<=moduleIndex; ++j) {
        const string& moduleLabel(moduleLabels[j]);
        const string  moduleType(hltConfig.moduleType(moduleLabel));
        // check whether the module is packed up in TriggerEvent product
        const unsigned int filterIndex(triggerSummary->filterIndex(InputTag(moduleLabel,"","HLT")));
        if (filterIndex<triggerSummary->sizeFilters()) {
          TString lable = moduleLabel.c_str();
          if (lable.Contains(hltJetFilterLabel_.label())) {

            const trigger::Vids& VIDS (triggerSummary->filterIds(filterIndex));
            const trigger::Keys& KEYS(triggerSummary->filterKeys(filterIndex));
            const size_type nI(VIDS.size());
            const size_type nK(KEYS.size());
            assert(nI==nK);
            const size_type n(max(nI,nK));
            const trigger::TriggerObjectCollection& TOC(triggerSummary->getObjects());
            for (size_type i=0; i!=n; ++i) {
              const trigger::TriggerObject& TO(TOC[KEYS[i]]);
              JetLegObjects.push_back(TO);	  
            }
          }
        }
      }
    }
  }

  //// Do subjets 
  std::vector<int> matchedjetIndices;
  if (doSubjets_) {
    iEvent.getByLabel(packedjLabel_, packedjetHandle);
    iEvent.getByLabel(sjLabel_, subjetHandle);
    if( packedjetHandle->size() > jetHandle->size() )
      edm::LogError("TooManyGroomedJets") << "Groomed packed jet collection size " << packedjetHandle->size()
        << " bigger than original fat jet collection size " << jetHandle->size() ;
    matchPackedJets(jetHandle,packedjetHandle,matchedjetIndices);
  }

  for (size_t i = 0; i< jetColl->size(); i++){
    pat::Jet & jet = (*jetColl)[i];
    // make fastjets out of daughters:
    std::vector<fastjet::PseudoJet> FJparticles;
    for (unsigned int k = 0; k < jet.numberOfDaughters(); k++) // Grabs all the constituents of the jet:
    {
      const edm::Ptr<reco::Candidate> & this_constituent = jet.daughterPtr(k);
      FJparticles.push_back( fastjet::PseudoJet( this_constituent->px(), this_constituent->py(), this_constituent->pz(), this_constituent->energy() ) );
    }
    fastjet::PseudoJet combJet = fastjet::join(FJparticles);
    Nsubjettiness t1(1, Njettiness::AxesMode::onepass_kt_axes, 1.0, 0.8);
    Nsubjettiness t2(2, Njettiness::AxesMode::onepass_kt_axes, 1.0, 0.8);
    Nsubjettiness t3(3, Njettiness::AxesMode::onepass_kt_axes, 1.0, 0.8);
    Nsubjettiness t4(4, Njettiness::AxesMode::onepass_kt_axes, 1.0, 0.8);
    double T1 = t1.result(combJet);
    double T2 = t2.result(combJet);
    double T3 = t3.result(combJet);
    double T4 = t4.result(combJet);
    jet.addUserFloat("tau1",   T1);
    jet.addUserFloat("tau2",   T2);
    jet.addUserFloat("tau3",   T3);
    jet.addUserFloat("tau4",   T4);
    // SUBJET POPULATOR:
    // DiJet Case:
    fastjet::JetDefinition jet_def_2(fastjet::antikt_algorithm, 1.5708);
    fastjet::ClusterSequence clust_seq_2(FJparticles, jet_def_2);
    int nSubJetsTwo = 2;
    std::vector<fastjet::PseudoJet> subjets2 = sorted_by_pt(clust_seq_2.exclusive_jets_up_to(nSubJetsTwo));
    double sjc2_j0_pt = subjets2[0].pt();
    double sjc2_j1_pt = subjets2[1].pt();
    jet.addUserFloat("sjc2_j0_pt",   sjc2_j0_pt);
    jet.addUserFloat("sjc2_j1_pt",   sjc2_j1_pt);
    double sjc2_j0_mass = subjets2[0].m();
    double sjc2_j1_mass = subjets2[1].m();
    jet.addUserFloat("sjc2_j0_mass",   sjc2_j0_mass);
    jet.addUserFloat("sjc2_j1_mass",   sjc2_j1_mass);
    double sjc2_j0_eta = subjets2[0].eta();
    double sjc2_j1_eta = subjets2[1].eta();
    double sjc2_j0_phi = subjets2[0].phi();
    double sjc2_j1_phi = subjets2[1].phi();
    jet.addUserFloat("sjc2_j0_eta",   sjc2_j0_eta);
    jet.addUserFloat("sjc2_j1_eta",   sjc2_j1_eta);
    jet.addUserFloat("sjc2_j0_phi",   sjc2_j0_phi);
    jet.addUserFloat("sjc2_j1_phi",   sjc2_j1_phi);
    // TriJet Case:
    fastjet::JetDefinition jet_def_3(fastjet::antikt_algorithm, 1.5708);
    fastjet::ClusterSequence clust_seq_3(FJparticles, jet_def_3);
    int nSubJetsThree = 3;
    std::vector<fastjet::PseudoJet> subjets3 = sorted_by_pt(clust_seq_3.exclusive_jets_up_to(nSubJetsThree));
    double sjc3_j0_pt = subjets3[0].pt();
    double sjc3_j1_pt = subjets3[1].pt();
    double sjc3_j2_pt = subjets3[2].pt();
    jet.addUserFloat("sjc3_j0_pt",   sjc3_j0_pt);
    jet.addUserFloat("sjc3_j1_pt",   sjc3_j1_pt);
    jet.addUserFloat("sjc3_j2_pt",   sjc3_j2_pt);
    double sjc3_j0_mass = subjets3[0].m();
    double sjc3_j1_mass = subjets3[1].m();
    double sjc3_j2_mass = subjets3[2].m();
    jet.addUserFloat("sjc3_j0_mass",   sjc3_j0_mass);
    jet.addUserFloat("sjc3_j1_mass",   sjc3_j1_mass);
    jet.addUserFloat("sjc3_j2_mass",   sjc3_j2_mass);
    double sjc3_j0_eta = subjets3[0].eta();
    double sjc3_j1_eta = subjets3[1].eta();
    double sjc3_j2_eta = subjets3[2].eta();
    double sjc3_j0_phi = subjets3[0].phi();
    double sjc3_j1_phi = subjets3[1].phi();
    double sjc3_j2_phi = subjets3[2].phi();
    jet.addUserFloat("sjc3_j0_eta",   sjc3_j0_eta);
    jet.addUserFloat("sjc3_j1_eta",   sjc3_j1_eta);
    jet.addUserFloat("sjc3_j2_eta",   sjc3_j2_eta);
    jet.addUserFloat("sjc3_j0_phi",   sjc3_j0_phi);
    jet.addUserFloat("sjc3_j1_phi",   sjc3_j1_phi);
    jet.addUserFloat("sjc3_j2_phi",   sjc3_j2_phi);
    // QuadJet Case:
    fastjet::JetDefinition jet_def_4(fastjet::antikt_algorithm, 1.5708);
    fastjet::ClusterSequence clust_seq_4(FJparticles, jet_def_4);
    int nSubJetsFour = 4;
    std::vector<fastjet::PseudoJet> subjets4 = sorted_by_pt(clust_seq_4.exclusive_jets_up_to(nSubJetsFour));
    double sjc4_j0_pt = subjets4[0].pt();
    double sjc4_j1_pt = subjets4[1].pt();
    double sjc4_j2_pt = subjets4[2].pt();
    double sjc4_j3_pt = subjets4[3].pt();
    jet.addUserFloat("sjc4_j0_pt",   sjc4_j0_pt);
    jet.addUserFloat("sjc4_j1_pt",   sjc4_j1_pt);
    jet.addUserFloat("sjc4_j2_pt",   sjc4_j2_pt);
    jet.addUserFloat("sjc4_j3_pt",   sjc4_j3_pt);
    double sjc4_j0_mass = subjets4[0].m();
    double sjc4_j1_mass = subjets4[1].m();
    double sjc4_j2_mass = subjets4[2].m();
    double sjc4_j3_mass = subjets4[3].m();
    jet.addUserFloat("sjc4_j0_mass",   sjc4_j0_mass);
    jet.addUserFloat("sjc4_j1_mass",   sjc4_j1_mass);
    jet.addUserFloat("sjc4_j2_mass",   sjc4_j2_mass);
    jet.addUserFloat("sjc4_j3_mass",   sjc4_j3_mass);
    double sjc4_j0_eta = subjets4[0].eta();
    double sjc4_j1_eta = subjets4[1].eta();
    double sjc4_j2_eta = subjets4[2].eta();
    double sjc4_j3_eta = subjets4[3].eta();
    double sjc4_j0_phi = subjets4[0].phi();
    double sjc4_j1_phi = subjets4[1].phi();
    double sjc4_j2_phi = subjets4[2].phi();
    double sjc4_j3_phi = subjets4[3].phi();
    jet.addUserFloat("sjc4_j0_eta",   sjc4_j0_eta);
    jet.addUserFloat("sjc4_j1_eta",   sjc4_j1_eta);
    jet.addUserFloat("sjc4_j2_eta",   sjc4_j2_eta);
    jet.addUserFloat("sjc4_j3_eta",   sjc4_j3_eta);
    jet.addUserFloat("sjc4_j0_phi",   sjc4_j0_phi);
    jet.addUserFloat("sjc4_j1_phi",   sjc4_j1_phi);
    jet.addUserFloat("sjc4_j2_phi",   sjc4_j2_phi);
    jet.addUserFloat("sjc4_j3_phi",   sjc4_j3_phi);

    // BTAGGING
    // - working points : https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagPerformanceOP
    // - SF : https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagPOG#Recommendation_for_b_c_tagging_a
    //float isCSVL = (jet.bDiscriminator("combinedInclusiveSecondaryVertexV2BJetTags") >  0.244 ? 1. : 0.);
    //float isCSVM = (jet.bDiscriminator("combinedInclusiveSecondaryVertexV2BJetTags") >  0.679 ? 1. : 0.);
    //float isCSVT = (jet.bDiscriminator("combinedInclusiveSecondaryVertexV2BJetTags") >  0.898 ? 1. : 0.);

    //Individual jet operations to be added here
    // trigger matched 
    int idx       = -1;
    double deltaR = -1.;
    bool isMatched2trigger = isMatchedWithTrigger(jet, JetLegObjects, idx, deltaR, hlt2reco_deltaRmax_) ;
    double hltEta = ( isMatched2trigger ? JetLegObjects[0].eta()    : -999.);
    double hltPhi = ( isMatched2trigger ? JetLegObjects[0].phi()    : -999.);
    double hltPt  = ( isMatched2trigger ? JetLegObjects[0].pt()     : -999.);
    double hltE   = ( isMatched2trigger ? JetLegObjects[0].energy() : -999.);

    // SMEARING
    // http://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
    reco::Candidate::LorentzVector smearedP4;
    if(isMC) {
      const reco::GenJet* genJet=jet.genJet();
      if(genJet) {
        float smearFactor=getResolutionRatio(jet.eta());
        smearedP4=jet.p4()-genJet->p4();
        smearedP4*=smearFactor; // +- 3*smearFactorErr;
        smearedP4+=genJet->p4();
      }
    } else {
      smearedP4=jet.p4();
    }
    // JER
    double JERup   = getJERup  (jet.eta());
    double JERdown = getJERdown(jet.eta());

    //jet.addUserFloat("IsCSVL",isCSVL);
    //jet.addUserFloat("IsCSVM",isCSVM);
    //jet.addUserFloat("IsCSVT",isCSVT);

    jet.addUserFloat("HLTjetEta",   hltEta);
    jet.addUserFloat("HLTjetPhi",   hltPhi);
    jet.addUserFloat("HLTjetPt",    hltPt);
    jet.addUserFloat("HLTjetE",     hltE);
    jet.addUserFloat("HLTjetDeltaR",deltaR);

    jet.addUserFloat("SmearedPEta", smearedP4.eta());
    jet.addUserFloat("SmearedPhi",  smearedP4.phi());
    jet.addUserFloat("SmearedPt",   smearedP4.pt());
    jet.addUserFloat("SmearedE",    smearedP4.energy());

    jet.addUserFloat("JERup", JERup);
    jet.addUserFloat("JERup", JERdown);


    TLorentzVector jetp4 ; 
    jetp4.SetPtEtaPhiE(jet.pt(), jet.eta(), jet.phi(), jet.energy()) ; 

    //// Lepton matching
    int matchedEl(-1), matchedMu(-1) ; 
    edm::Ptr<reco::Candidate> elobj, muobj ; 
    int elIdx(-1), muIdx(-1) ; 

    double drmin = 999 ; 
    for (std::vector<pat::Electron>::const_iterator iel = elHandle->begin(); iel != elHandle->end(); ++iel) {
      TLorentzVector elp4 ; 
      elp4.SetPtEtaPhiE(iel->pt(), iel->eta(), iel->phi(), iel->energy()) ; 
      if (elp4.DeltaR(jetp4) < drmin) {
        drmin = elp4.DeltaR(jetp4) ; 
        elobj = iel->originalObjectRef() ; 
        elIdx = std::distance(elHandle->begin(),iel) ; 
      }
    }


    drmin = 999 ; 
    for (std::vector<pat::Muon>::const_iterator imu = muHandle->begin(); imu != muHandle->end(); ++imu) {
      TLorentzVector mup4 ; 
      mup4.SetPtEtaPhiE(imu->pt(), imu->eta(), imu->phi(), imu->energy()) ; 
      if (mup4.DeltaR(jetp4) < drmin) {
        drmin = mup4.DeltaR(jetp4) ; 
        muobj = imu->originalObjectRef() ; 
        muIdx = std::distance(muHandle->begin(),imu) ; 
      }
    }

    const std::vector<edm::Ptr<reco::Candidate> > daus = jet.daughterPtrVector() ; 
    for (std::vector<edm::Ptr<reco::Candidate> >::const_iterator idau = daus.begin(); idau != daus.end(); ++idau) {
      if (*idau == elobj) matchedEl = elIdx ; 
      if (*idau == muobj) matchedMu = muIdx ; 
    }

    jet.addUserFloat("matchedMuIdx", matchedMu);
    jet.addUserFloat("matchedElIdx", matchedEl);

    //// Do subjets
    if (doSubjets_) {
      int subjet0Id = -1, subjet1Id = -1 ;
      int pfjIdx = matchedjetIndices.at(i) ; //// Get index of matched packed jet
      int nSubjets = 0 ;
      std::vector<PatJetCollection::const_iterator> itsubjets ;
      if (pfjIdx >= 0) {
        nSubjets = packedjetHandle->at(pfjIdx).numberOfDaughters() ;
        const edm::Ptr<reco::Jet> originalObjRef = edm::Ptr<reco::Jet>( packedjetHandle->at(pfjIdx).originalObjectRef() );
        for (int isj = 0; isj < nSubjets; ++isj) {
          for (std::vector<pat::Jet>::const_iterator itsj = subjetHandle->begin(); itsj != subjetHandle->end(); ++itsj) {
            if ( originalObjRef->daughterPtr(isj) == itsj->originalObjectRef() ) {
              itsubjets.push_back(itsj) ;
            }
          }
        } //// Loop over subjets
      } //// If subjets 
      std::sort(itsubjets.begin(), itsubjets.end(), orderByPt("Uncorrected"));
      if ( itsubjets.size() > 1 ) {
        subjet0Id = itsubjets.at(0) - subjetHandle->begin() ;
        subjet1Id = itsubjets.at(1) - subjetHandle->begin() ;
        jet.addUserInt("subjetIndex0", subjet0Id) ;
        jet.addUserInt("subjetIndex1", subjet1Id) ;
        //std::cout << " 2 subjet found with indices " << subjet0Id << " and " << subjet1Id << "\n" ;
      }
    } //// Do subjets 

  } //// Loop over all jets 

  iEvent.put( jetColl );

}

// ------------ method called once each job just after ending the event loop  ------------
  bool
JetUserData::isMatchedWithTrigger(const pat::Jet& p, trigger::TriggerObjectCollection triggerObjects, int& index, double& deltaR, double deltaRmax = 0.2)
{
  for (size_t i = 0 ; i < triggerObjects.size() ; i++){
    float dR = sqrt(pow(triggerObjects[i].eta()-p.eta(),2)+ pow(acos(cos(triggerObjects[i].phi()-p.phi())),2)) ;
    if (dR<deltaRmax) {
      deltaR = dR;
      index  = i;
      return true;
    }
  }
  return false;
}

  double
JetUserData::getResolutionRatio(double eta)
{
  eta=fabs(eta);
  if(eta>=0.0 && eta<0.5) return 1.079; // +-0.005 +-0.026 
  if(eta>=0.5 && eta<1.1) return 1.099; // +-0.005 +-0.028 
  if(eta>=1.1 && eta<1.7) return 1.121; // +-0.005 +-0.029 
  if(eta>=1.7 && eta<2.3) return 1.208; // +-0.013 +-0.045 
  if(eta>=2.3 && eta<2.8) return 1.254; // +-0.026 +-0.056 
  if(eta>=2.8 && eta<3.2) return 1.395; // +-0.036 +-0.051 
  if(eta>=3.2 && eta<5.0) return 1.056; // +-0.048 +-0.185 
  return -1.;
}

  double
JetUserData::getJERup(double eta)
{
  eta=fabs(eta);
  if(eta>=0.0 && eta<0.5) return 1.053 ;
  if(eta>=0.5 && eta<1.1) return 1.071 ;
  if(eta>=1.1 && eta<1.7) return 1.092 ;
  if(eta>=1.7 && eta<2.3) return 1.162 ;
  if(eta>=2.3 && eta<2.8) return 1.192 ;
  if(eta>=2.8 && eta<3.2) return 1.332 ;
  if(eta>=3.2 && eta<5.0) return 0.865 ;
  return -1.;  
}

  double
JetUserData::getJERdown(double eta)
{
  eta=fabs(eta);
  if(eta>=0.0 && eta<0.5) return 1.105 ;
  if(eta>=0.5 && eta<1.1) return 1.127 ;
  if(eta>=1.1 && eta<1.7) return 1.150 ;
  if(eta>=1.7 && eta<2.3) return 1.254 ;
  if(eta>=2.3 && eta<2.8) return 1.316 ;
  if(eta>=2.8 && eta<3.2) return 1.458 ;
  if(eta>=3.2 && eta<5.0) return 1.247 ;
  return -1.;  
}

// ------------ method that matches packed and original jets based on minimum dR ------------
void JetUserData::matchPackedJets(const edm::Handle<PatJetCollection>& jets,
    const edm::Handle<PatJetCollection>& packedJets, std::vector<int>& matchedIndices) {

  std::vector<bool> jetLocks(jets->size(),false);
  std::vector<int>  jetIndices;

  for(size_t pj=0; pj<packedJets->size(); ++pj) {
    double matchedDR = 1e9;
    int matchedIdx = -1;

    for(size_t ij = 0; ij < jets->size(); ++ij) {
      if( jetLocks.at(ij) ) continue; // skip jets that have already been matched

      double tempDR = reco::deltaR( jets->at(ij).rapidity(), jets->at(ij).phi(), packedJets->at(pj).rapidity(), packedJets->at(pj).phi() );
      if( tempDR < matchedDR ) {
        matchedDR = tempDR;
        matchedIdx = ij;
      }
    }

    if( matchedIdx >= 0 ) jetLocks.at(matchedIdx) = true;
    jetIndices.push_back(matchedIdx);
  }

  if( std::find( jetIndices.begin(), jetIndices.end(), -1 ) != jetIndices.end() )
    edm::LogError("JetMatchingFailed") << "Matching groomed to original jets failed. Please check that the two jet collections belong to each other.";

  for(size_t ij = 0; ij < jets->size(); ++ij) {
    std::vector<int>::iterator matchedIndex = std::find( jetIndices.begin(), jetIndices.end(), ij );
    matchedIndices.push_back( matchedIndex != jetIndices.end() ? std::distance(jetIndices.begin(),matchedIndex) : -1 );
  }
}


#include "FWCore/Framework/interface/MakerMacros.h"


DEFINE_FWK_MODULE(JetUserData);
