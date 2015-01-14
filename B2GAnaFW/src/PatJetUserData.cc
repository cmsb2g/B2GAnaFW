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

class PatJetUserData : public edm::EDProducer {
public:
  PatJetUserData( const edm::ParameterSet & );   

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


PatJetUserData::PatJetUserData(const edm::ParameterSet& iConfig) :
 
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

void PatJetUserData::produce( edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  //Jets
  edm::Handle<std::vector<pat::Jet> > jetHandle;
  iEvent.getByLabel(jLabel_, jetHandle);

  auto_ptr<vector<pat::Jet> > jetColl( new vector<pat::Jet> (*jetHandle) );

  for ( std::vector<pat::Jet>::const_iterator jetBegin = jetColl->begin(),
	  jetEnd = jetColl->end(), jet = jetBegin; jet != jetEnd; ++jet ) {

    pat::Jet* subjet1 = dynamic_cast<pat::Jet*>((jet->daughter(0)));
    double subjet1Bdisc = subjet1->bDiscriminator("combinedInclusiveSecondaryVertexV2BJetTags");
    std::cout << "success..." << subjet1Bdisc << std::endl;
    pat::Jet* subjet2 = dynamic_cast<pat::Jet*>((jet->daughter(1)));
    double subjet2Bdisc = subjet2->bDiscriminator("combinedInclusiveSecondaryVertexV2BJetTags");
    std::cout << "success..." << subjet2Bdisc << std::endl;
    pat::Jet* subjet3 = dynamic_cast<pat::Jet*>((jet->daughter(2)));
    double subjet3Bdisc = subjet3->bDiscriminator("combinedInclusiveSecondaryVertexV2BJetTags");
    std::cout << "success..." << subjet3Bdisc << std::endl;

    jet->addUserFloat("subjet1csv", subjet1Bdisc);
    jet->addUserFloat("subjet2csv", subjet2Bdisc);
    jet->addUserFloat("subjet3csv", subjet3Bdisc);

  }

  iEvent.put( jetColl );

}

// ------------ method called once each job just after ending the event loop  ------------
bool
PatJetUserData::isMatchedWithTrigger(const pat::Jet& p, trigger::TriggerObjectCollection triggerObjects, int& index, double& deltaR, double deltaRmax = 0.2)
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
PatJetUserData::getResolutionRatio(double eta)
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
PatJetUserData::getJERup(double eta)
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
PatJetUserData::getJERdown(double eta)
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


DEFINE_FWK_MODULE(PatJetUserData);
