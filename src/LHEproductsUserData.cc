#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDMException.h"

#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include <Math/VectorUtil.h>


using namespace edm;
using namespace std;
using namespace lhef;

class  LHEUserData : public edm::EDProducer {
public:
  LHEUserData( const edm::ParameterSet & );   

private:
  void produce( edm::Event &, const edm::EventSetup & );
  InputTag LHELabel_;
 };


LHEUserData::LHEUserData(const edm::ParameterSet& iConfig):
   LHELabel_(iConfig.getParameter<edm::InputTag>("lheLabel"))
 {
   // produces<std::vector<pat::Muon> >();
  produces<int>( "lheNup" ).setBranchAlias( "lheNup" );
  produces<vector<int> >( "lheIDup" ).setBranchAlias( "lheIDup" );
  //  produces<vector<double> >( "lhePxup" ).setBranchAlias( "lhePxup" );
  //  produces<vector<double> >( "lhePyup" ).setBranchAlias( "lhePyup" );
  //  produces<vector<double> >( "lhePzup" ).setBranchAlias( "lhePzup" );
  //  produces<vector<double> >( "lheEup" ).setBranchAlias( "lheEup" );
  produces<std::vector< math::XYZTLorentzVector > >( "lhePup" ).setBranchAlias( "lhePup" );
 }


void LHEUserData::produce( edm::Event& iEvent, const edm::EventSetup& iSetup) {
    
  edm::Handle<LHEEventProduct> lheHandle; //to be replaced with PAT
  iEvent.getByLabel(LHELabel_, lheHandle);

  auto_ptr<int> lheNup( new int );
  auto_ptr<vector<int> > lheIDup( new vector<int> );
  //  auto_ptr<vector<double> > lhePx( new vector<double> );
  //  auto_ptr<vector<double> > lhePy( new vector<double> );
  //  auto_ptr<vector<double> > lhePz( new vector<double> );
  //  auto_ptr<vector<double> > lheE( new vector<double> );
  auto_ptr< std::vector< math::XYZTLorentzVector > > lhePup( new std::vector< math::XYZTLorentzVector > );
  *lheNup = -1;
  
  const lhef::HEPEUP hepeup_ = lheHandle->hepeup();
  *lheNup = hepeup_.NUP; 

  for(int i = 0; i<*lheNup ; ++i){
    lheIDup->push_back(hepeup_.IDUP[i]);
    // lhePx->push_back(hepeup_.PUP[i][0]);
    // lhePy->push_back(hepeup_.PUP[i][1]);
    // lhePz->push_back(hepeup_.PUP[i][2]);
    // lheE->push_back(hepeup_.PUP[i][3]);
 math::XYZTLorentzVector p4;
    p4.SetPxPyPzE(hepeup_.PUP[i][0], hepeup_.PUP[i][1], hepeup_.PUP[i][2], hepeup_.PUP[i][3]);
    lhePup->push_back(p4);
  }
  

  iEvent.put( lheNup, "lheNup" );
  iEvent.put( lheIDup, "lheIDup" );
  //  iEvent.put( lhePx, "lhePx" );
  //  iEvent.put( lhePy, "lhePy" );
  //  iEvent.put( lhePz, "lhePz" );
  //  iEvent.put( lheE, "lheE" );
  iEvent.put( lhePup, "lhePup" );


}

// ------------ method called once each job just after ending the event loop  ------------


#include "FWCore/Framework/interface/MakerMacros.h"


DEFINE_FWK_MODULE(LHEUserData);
