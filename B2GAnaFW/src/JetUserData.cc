#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDMException.h"


#include "DataFormats/PatCandidates/interface/Jet.h"

// dR and dPhi
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"

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
#include <vector>

using namespace reco;
using namespace edm;
using namespace std;
using namespace trigger;

class JetUserData : public edm::EDProducer {
public:
  JetUserData( const edm::ParameterSet & );   

private:
  void produce( edm::Event &, const edm::EventSetup & );
  bool isMatchedWithTrigger(const pat::Jet&, trigger::TriggerObjectCollection,int&,double&,double);
  double getResolutionRatio(double eta);
  double getJERup(double eta);
  double getJERdown(double eta);

  edm::EDGetTokenT<std::vector<pat::Jet> >     jetToken_;
  edm::EDGetTokenT<std::vector<reco::Vertex> > pvToken_;

  //InputTag jetLabel_;
  InputTag jLabel_, pvLabel_;
  InputTag triggerResultsLabel_, triggerSummaryLabel_;
  InputTag hltJetFilterLabel_;
  std::string hltPath_;
  double hlt2reco_deltaRmax_;
  HLTConfigProvider hltConfig;
  int triggerBit;
 };


JetUserData::JetUserData(const edm::ParameterSet& iConfig) :
 
   jLabel_(iConfig.getParameter<edm::InputTag>("jetLabel")),
   pvLabel_(iConfig.getParameter<edm::InputTag>("pv")),   // "offlinePrimaryVertex"

   triggerResultsLabel_(iConfig.getParameter<edm::InputTag>("triggerResults")),
   triggerSummaryLabel_(iConfig.getParameter<edm::InputTag>("triggerSummary")),
   hltJetFilterLabel_  (iConfig.getParameter<edm::InputTag>("hltJetFilter")),   //trigger objects we want to match
   hltPath_            (iConfig.getParameter<std::string>("hltPath")),
   hlt2reco_deltaRmax_ (iConfig.getParameter<double>("hlt2reco_deltaRmax"))
  {
    produces<vector<pat::Jet> >();
  }

#include "DataFormats/JetReco/interface/Jet.h"

void JetUserData::produce( edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  bool isMC = (!iEvent.isRealData());
  
  //PV
  edm::Handle<std::vector<reco::Vertex> > pvHandle;
  iEvent.getByLabel(pvLabel_, pvHandle);
  //const reco::Vertex& PV= pvHandle->front();
  
  //Jets
  edm::Handle<std::vector<pat::Jet> > jetHandle;
  iEvent.getByLabel(jLabel_, jetHandle);
  auto_ptr<vector<pat::Jet> > jetColl( new vector<pat::Jet> (*jetHandle) );

  /////////  /////////  /////////  /////////  /////////  /////////  /////////  /////////  /////////
  // TRIGGER (this is not really needed ...)
  bool changedConfig = false;
  bool pathFound = false;
  if (!hltConfig.init(iEvent.getRun(), iSetup, "HLT", changedConfig)) {
    std::cout << "Initialization of HLTConfigProvider failed!!" << std::endl;
    return;
  }

  if (changedConfig){
    std::cout << "the current menu is " << hltConfig.tableName() << std::endl;
    triggerBit = -1;
    for (size_t j = 0; j < hltConfig.triggerNames().size(); j++) {
      //std::cout << "hltConfig.triggerNames()[" << j << "]: " << hltConfig.triggerNames()[j] << std::endl;
      if (TString(hltConfig.triggerNames()[j]).Contains(hltPath_)) {triggerBit = j;pathFound=true;}
    }
    if (triggerBit == -1) std::cout << "HLT path not found" << std::endl;
  }

  edm::Handle<edm::TriggerResults> triggerResults;
  iEvent.getByLabel(triggerResultsLabel_, triggerResults);
  if (size_t(triggerBit) < triggerResults->size() && pathFound)
    if (triggerResults->accept(triggerBit))
      std::cout << "event pass : " << hltPath_ << std::endl;

  /////////  /////////  /////////  /////////  /////////  /////////  /////////  /////////  /////////
  // TRIGGER MATCHING
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
	  //      cout << " 'L3' filter in slot " << j << " - label/type " << moduleLabel << "/" << moduleType << endl;
	  TString lable = moduleLabel.c_str();
	  if (lable.Contains(hltJetFilterLabel_.label())) {
	    
	    const trigger::Vids& VIDS (triggerSummary->filterIds(filterIndex));
	    const trigger::Keys& KEYS(triggerSummary->filterKeys(filterIndex));
	    const size_type nI(VIDS.size());
	    const size_type nK(KEYS.size());
	    assert(nI==nK);
	    const size_type n(max(nI,nK));
	    //	cout << "   " << n  << " accepted TRIGGER objects found: " << endl;
	    const trigger::TriggerObjectCollection& TOC(triggerSummary->getObjects());
	    for (size_type i=0; i!=n; ++i) {
	      const trigger::TriggerObject& TO(TOC[KEYS[i]]);
	      JetLegObjects.push_back(TO);	  
	      //	  cout << "   " << i << " " << VIDS[i] << "/" << KEYS[i] << ": "
	      //	       << TO.id() << " " << TO.pt() << " " << TO.eta() << " " << TO.phi() << " " << TO.mass()
	      //	       << endl;
	    }
	  }
	}
      }
    }
  }
  /////////  /////////  /////////  /////////  /////////  /////////  /////////  /////////  /////////  /////////

  //std::cout << " JetUserData 2 " << std::endl;
  //std::cout << " jetColl " << jetColl->size() << std::endl;

  for (size_t i = 0; i< jetColl->size(); i++){
    pat::Jet & jet = (*jetColl)[i];
    
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
    
  }

  iEvent.put( jetColl );

}

// ------------ method called once each job just after ending the event loop  ------------
bool
JetUserData::isMatchedWithTrigger(const pat::Jet& p, trigger::TriggerObjectCollection triggerObjects, int& index, double& deltaR, double deltaRmax = 0.2)
{
  for (size_t i = 0 ; i < triggerObjects.size() ; i++){
    float dR = sqrt(pow(triggerObjects[i].eta()-p.eta(),2)+ pow(acos(cos(triggerObjects[i].phi()-p.phi())),2)) ;
    //    std::cout << "dR: " << dR << std::endl;
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


#include "FWCore/Framework/interface/MakerMacros.h"


DEFINE_FWK_MODULE(JetUserData);
