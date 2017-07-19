#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/DependentRecordImplementation.h"
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

// Vertex
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

// trigger
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h" // gives access to the (release cycle dependent) trigger object codes
#include "DataFormats/JetReco/interface/Jet.h"

// SVTagInfo
#include "DataFormats/BTauReco/interface/SecondaryVertexTagInfo.h"

// JEC/JER
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"
#include <TRandom3.h>

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

#define DEBUG false

typedef std::vector<pat::Jet> PatJetCollection;

class JetUserData : public edm::EDProducer {
  public:
    JetUserData( const edm::ParameterSet & );   

  private:
    void produce( edm::Event &, const edm::EventSetup & );
    bool isMatchedWithTrigger(const pat::Jet&, trigger::TriggerObjectCollection,int&,double&,double);

    edm::EDGetTokenT<std::vector<pat::Jet> >     jetToken_;

    //InputTag jetLabel_;
    EDGetTokenT< std::vector< pat::Jet > > jLabel_;
    EDGetTokenT<double> rhoLabel_;
    double coneSize_;
    bool getJERFromTxt_;
    std::string jetCorrLabel_;
    std::string jerLabel_;
    std::string resolutionsFile_;
    std::string scaleFactorsFile_;
    EDGetTokenT< edm::TriggerResults > triggerResultsLabel_;
    EDGetTokenT< trigger::TriggerEvent > triggerSummaryLabel_;
    InputTag hltJetFilterLabel_;
    std::string hltPath_;
    double hlt2reco_deltaRmax_;
    std::string candSVTagInfos_;
    HLTConfigProvider hltConfig;
    int triggerBit;
    TRandom3 rnd_;
};


JetUserData::JetUserData(const edm::ParameterSet& iConfig) :
  jLabel_             (consumes<std::vector<pat::Jet>>(iConfig.getParameter<edm::InputTag>("jetLabel"))), 
  rhoLabel_           (consumes<double>(iConfig.getParameter<edm::InputTag>("rho"))),
  coneSize_           (iConfig.getParameter<double>("coneSize")),
  getJERFromTxt_      (iConfig.getParameter<bool>("getJERFromTxt")),
  jetCorrLabel_       (iConfig.getParameter<std::string>("jetCorrLabel")),
  triggerResultsLabel_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerResults"))),
  triggerSummaryLabel_(consumes<trigger::TriggerEvent>(iConfig.getParameter<edm::InputTag>("triggerSummary"))),
  hltJetFilterLabel_  (iConfig.getParameter<edm::InputTag>("hltJetFilter")),   //trigger objects we want to match
  hltPath_            (iConfig.getParameter<std::string>("hltPath")),
  hlt2reco_deltaRmax_ (iConfig.getParameter<double>("hlt2reco_deltaRmax")),
  candSVTagInfos_         (iConfig.getParameter<std::string>("candSVTagInfos"))
{
  if (getJERFromTxt_) {
    resolutionsFile_  = iConfig.getParameter<std::string>("resolutionsFile");
    scaleFactorsFile_ = iConfig.getParameter<std::string>("scaleFactorsFile");
  } else
    jerLabel_         = iConfig.getParameter<std::string>("jerLabel");
  produces<vector<pat::Jet> >();
}


void JetUserData::produce( edm::Event& iEvent, const edm::EventSetup& iSetup) {

  bool isMC = (!iEvent.isRealData());

  edm::Handle<std::vector<pat::Jet> > jetHandle, packedjetHandle;
  iEvent.getByToken(jLabel_, jetHandle);
  std::unique_ptr<vector<pat::Jet> > jetColl( new vector<pat::Jet> (*jetHandle) );


  //// TRIGGER (this is not really needed ...)
  //bool changedConfig = false;
  //bool pathFound = false;
  //if (!hltConfig.init(iEvent.getRun(), iSetup, "HLT", changedConfig)) {
  //  edm::LogError("HLTConfigProvider") << "Initialization of HLTConfigProvider failed!!" << std::endl;
  //  return;
  //}

  //if (changedConfig){
  //  edm::LogInfo("HLTMenu") << "the current menu is " << hltConfig.tableName() << std::endl;
  //  triggerBit = -1;
  //  for (size_t j = 0; j < hltConfig.triggerNames().size(); j++) {
  //    if (TString(hltConfig.triggerNames()[j]).Contains(hltPath_)) {triggerBit = j;pathFound=true;}
  //  }
  //  if (triggerBit == -1) edm::LogError("NoHLTPath") << "HLT path not found" << std::endl;
  //}

  //   edm::Handle<edm::TriggerResults> triggerResults;
  //   iEvent.getByToken(triggerResultsLabel_, triggerResults);
  /* Why do we need this
     if (size_t(triggerBit) < triggerResults->size() && pathFound)
       if (triggerResults->accept(triggerBit))
         std::cout << "event pass : " << hltPath_ << std::endl;
  */ 

  //// TRIGGER MATCHING
  //trigger::TriggerObjectCollection JetLegObjects;

  //edm::Handle<trigger::TriggerEvent> triggerSummary;

  //if ( triggerSummary.isValid() ) {
  //  iEvent.getByToken(triggerSummaryLabel_, triggerSummary);

  //  // Results from TriggerEvent product - Attention: must look only for
  //  // modules actually run in this path for this event!
  //  if(pathFound){
  //    const unsigned int triggerIndex(hltConfig.triggerIndex(hltPath_));
  //    const vector<string>& moduleLabels(hltConfig.moduleLabels(triggerIndex));
  //    const unsigned int moduleIndex(triggerResults->index(triggerIndex));
  //    for (unsigned int j=0; j<=moduleIndex; ++j) {
  //      const string& moduleLabel(moduleLabels[j]);
  //      const string  moduleType(hltConfig.moduleType(moduleLabel));
  //      // check whether the module is packed up in TriggerEvent product
  //      const unsigned int filterIndex(triggerSummary->filterIndex(InputTag(moduleLabel,"","HLT")));
  //      if (filterIndex<triggerSummary->sizeFilters()) {
  //        TString lable = moduleLabel.c_str();
  //        if (lable.Contains(hltJetFilterLabel_.label())) {

  //          const trigger::Vids& VIDS (triggerSummary->filterIds(filterIndex));
  //          const trigger::Keys& KEYS(triggerSummary->filterKeys(filterIndex));
  //          const size_type nI(VIDS.size());
  //          const size_type nK(KEYS.size());
  //          assert(nI==nK);
  //          const size_type n(max(nI,nK));
  //          const trigger::TriggerObjectCollection& TOC(triggerSummary->getObjects());
  //          for (size_type i=0; i!=n; ++i) {
  //            const trigger::TriggerObject& TO(TOC[KEYS[i]]);
  //            JetLegObjects.push_back(TO);	  
  //          }
  //        }
  //      }
  //    }
  //  }
  //}

  // JEC Uncertainty
  edm::ESHandle<JetCorrectorParametersCollection> JetCorrParColl;
  iSetup.get<JetCorrectionsRecord>().get(jetCorrLabel_, JetCorrParColl); 
  JetCorrectorParameters const & JetCorrPar = (*JetCorrParColl)["Uncertainty"];
  JetCorrectionUncertainty jecUnc(JetCorrPar);

  // JER
  // Twiki: https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyResolution#Scale_factors
  // Recipe taken from: https://github.com/blinkseb/cmssw/blob/jer_fix_76x/JetMETCorrections/Modules/plugins/JetResolutionDemo.cc
  edm::Handle<double> rho;
  iEvent.getByToken(rhoLabel_, rho);
  JME::JetParameters jetParam;
  JME::JetResolution resolution;
  JME::JetResolutionScaleFactor res_sf;
  if (getJERFromTxt_) {
    resolution = JME::JetResolution(resolutionsFile_);
    res_sf = JME::JetResolutionScaleFactor(scaleFactorsFile_);
  } else {
    resolution = JME::JetResolution::get(iSetup, jerLabel_+"_pt");
    res_sf = JME::JetResolutionScaleFactor::get(iSetup, jerLabel_);
  }

  for (size_t i = 0; i< jetColl->size(); i++){
    pat::Jet & jet = (*jetColl)[i];

    //Individual jet operations to be added here
    // trigger matched 
    //int idx       = -1;
    //double deltaR = -1.;
    //bool isMatched2trigger = isMatchedWithTrigger(jet, JetLegObjects, idx, deltaR, hlt2reco_deltaRmax_) ;
    //double hltEta = ( isMatched2trigger ? JetLegObjects[0].eta()    : -999.);
    //double hltPhi = ( isMatched2trigger ? JetLegObjects[0].phi()    : -999.);
    //double hltPt  = ( isMatched2trigger ? JetLegObjects[0].pt()     : -999.);
    //double hltE   = ( isMatched2trigger ? JetLegObjects[0].energy() : -999.);

    // JEC uncertainty
    jecUnc.setJetPt (jet.pt());
    jecUnc.setJetEta(jet.eta());
    double jecUncertainty = jecUnc.getUncertainty(true);

    // JER
    jetParam.setJetPt(jet.pt()).setJetEta(jet.eta()).setRho(*rho);
    float PtResolution = resolution.getResolution(jetParam);
    float JERSF        = res_sf.getScaleFactor(jetParam);
    float JERSFUp      = res_sf.getScaleFactor(jetParam, Variation::UP);
    float JERSFDown    = res_sf.getScaleFactor(jetParam, Variation::DOWN);

    // Hybrid scaling and smearing procedure applied:
    //   https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution#Smearing_procedures
    reco::Candidate::LorentzVector smearedP4 =jet.p4();
    if(isMC) {
      // Hybrid method: Scale jet four momentum for well matched jets ...
      bool isGenMatched = 0;
      const reco::GenJet* genJet=jet.genJet();
      if (genJet) {
	TLorentzVector jetp4, genjetp4;
	jetp4.SetPtEtaPhiE(jet.pt(), jet.eta(), jet.phi(), jet.energy());
	genjetp4.SetPtEtaPhiE(genJet->pt(), genJet->eta(), genJet->phi(), genJet->energy());
	float dR = jetp4.DeltaR(genjetp4);
	float dPt = jet.pt()-genJet->pt();
	if ((dR<coneSize_/2.0)&&(std::abs(dPt)<(3*PtResolution*jet.pt()))) {
	  isGenMatched = 1;
	  smearedP4     *= std::max(0., 1 + (JERSF     - 1) * dPt / jet.pt());
	}
      }
      // ... and gaussian smear the rest
      if (!isGenMatched && JERSF>1) {
	double sigma = std::sqrt(JERSF * JERSF - 1) * PtResolution ;
	smearedP4 *= 1 + rnd_.Gaus(0, sigma);
      }
    }

    jet.addUserFloat("jecUncertainty",   jecUncertainty);

    //jet.addUserFloat("HLTjetEta",   hltEta);
    //jet.addUserFloat("HLTjetPhi",   hltPhi);
    //jet.addUserFloat("HLTjetPt",    hltPt);
    //jet.addUserFloat("HLTjetE",     hltE);
    //jet.addUserFloat("HLTjetDeltaR",deltaR);

    jet.addUserFloat("PtResolution", PtResolution);
    jet.addUserFloat("JERSF",        JERSF);
    jet.addUserFloat("JERSFUp",      JERSFUp);
    jet.addUserFloat("JERSFDown",    JERSFDown);
    jet.addUserFloat("SmearedPt",    smearedP4.pt());
    jet.addUserFloat("SmearedE",     smearedP4.energy());

    unsigned int nSV(0);
    float SV0mass(-999), SV1mass(-999) ;

    std::vector<std::string>tagInfoLabels = jet.tagInfoLabels() ;
#if DEBUG
    bool hasCandSVTagInfo(jet.hasTagInfo(candSVTagInfos_)) ; 
    std::cout << " jetTagInfoLabels size = " << tagInfoLabels.size() << std::endl ; 
    for (std::string tagInfoLabel : tagInfoLabels ) {
      std::cout << ">>>> Jet has " << tagInfoLabel << std::endl ; 
    }
    std::cout << ">>>>>> candSVTagInfo label is " << candSVTagInfos_ << std::endl ; 
    std::cout << " hasCandSVTagInfo = " << hasCandSVTagInfo << std::endl ; 
#endif 

    if ( jet.hasTagInfo(candSVTagInfos_) ) {
      const reco::CandSecondaryVertexTagInfo *candSVTagInfo = jet.tagInfoCandSecondaryVertex("pfInclusiveSecondaryVertexFinder");
#if DEBUG
      if ( candSVTagInfo == nullptr ) std::cout << ">>>>>> candSVTagInfo ptr does not exist\n" ;
      else std::cout << ">>>>>> candSVTagInfo ptr exists\n" ;
#endif 
      nSV = candSVTagInfo->nVertices() ; 
      SV0mass = nSV > 0 ? ((candSVTagInfo->secondaryVertex(0)).p4()).mass() : -999 ;
      SV1mass = nSV > 1 ? ((candSVTagInfo->secondaryVertex(1)).p4()).mass() : -999 ;
      if ( nSV > 0 ) {
#if DEBUG
        std::cout << ">>>>> nSV = " << nSV 
          << " SV0 mass = " << SV0mass 
          << " SV1 mass = " << SV1mass 
          << std::endl ; 
#endif 
      }
    }

    jet.addUserInt("nSV"     , nSV     ); 
    jet.addUserFloat("SV0mass", SV0mass); 
    jet.addUserFloat("SV1mass", SV1mass); 

    //// Jet constituent indices for lepton matching
    std::vector<unsigned int> constituentIndices;
    auto constituents = jet.daughterPtrVector();
    for ( auto & constituent : constituents ) {
      constituentIndices.push_back( constituent.key() );
    }

    jet.addUserData("pfKeys", constituentIndices );


  } //// Loop over all jets 

  iEvent.put( std::move(jetColl) );

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

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(JetUserData);
