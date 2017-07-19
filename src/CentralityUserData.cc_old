#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDMException.h"

#include "DataFormats/Math/interface/Vector3D.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"

#include<vector>

using namespace reco;
using namespace edm;
using namespace std;


class  CentralityUserData : public edm::EDProducer {
public:
  CentralityUserData( const edm::ParameterSet & );   

private:
  void produce( edm::Event &, const edm::EventSetup & );
  void put( edm::Event& evt, double value, const char* instanceName);

  edm::EDGetTokenT<edm::View<reco::Candidate> > srcToken_;

 };


CentralityUserData::CentralityUserData(const edm::ParameterSet& iConfig)
 {
   srcToken_ = consumes<edm::View<reco::Candidate> >(iConfig.getParameter<edm::InputTag>("src"));
   produces<double>("centrality");
 }

void CentralityUserData::put(edm::Event& evt, double value, const char* instanceName)
{
  std::auto_ptr<double> varPtr(new double(value));
  evt.put(varPtr, instanceName);
}

#include "TLorentzVector.h"
void CentralityUserData::produce( edm::Event& iEvent, const edm::EventSetup& iSetup) {
    
  edm::Handle<edm::View<reco::Candidate> > objects;
  iEvent.getByToken(srcToken_, objects);

  std::vector<TLorentzVector> boosted_objects;
  TLorentzVector all(0.,0.,0.,0.);
  for ( edm::View<reco::Candidate>::const_iterator vec = objects->begin(); vec != objects->end(); ++vec) {
    TLorentzVector tmp = TLorentzVector(vec->px(), vec->py(), vec->pz(), vec->energy());
    all+=(tmp);
    boosted_objects.push_back( tmp );
  }
  double sum_pt = 0.;
  double sum_E = 0.;
  for ( size_t i=0; i< boosted_objects.size(); i++ ) {
    boosted_objects[i].Boost( -all.BoostVector() );
    sum_pt+=boosted_objects[i].Pt();
    sum_E +=boosted_objects[i].Energy();
  }

  double value = sum_pt/sum_E;
  double centrality = ( value!=value || std::isinf(value) || value < 0) ? 0. : value;
  
  put(iEvent, centrality, "centrality");
  
}

// ------------ method called once each job just after ending the event loop  ------------


#include "FWCore/Framework/interface/MakerMacros.h"


DEFINE_FWK_MODULE(CentralityUserData);
