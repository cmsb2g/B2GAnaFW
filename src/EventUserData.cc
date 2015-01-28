#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDMException.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include <Math/VectorUtil.h>


using namespace edm;
using namespace std;

class  EventUserData : public edm::EDProducer {
public:
  EventUserData( const edm::ParameterSet & );   

private:
  void produce( edm::Event &, const edm::EventSetup & );
  //InputTag LHELabel_;
  edm::InputTag m_PileupSrc;

 };


EventUserData::EventUserData(const edm::ParameterSet& iConfig)
 {
   m_PileupSrc = iConfig.getUntrackedParameter<edm::InputTag>("pileup",edm::InputTag("addPileupInfo"));
   
   produces<std::vector<int> >("puBX");
   produces<std::vector<int> >("puNInt");
   produces<int>("puNtrueInt");
   
 
 }


void EventUserData::produce( edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  auto_ptr<int> puNtrueInt (new int);
  auto_ptr<vector<int> > puBX(  new vector<int> );
  auto_ptr<vector<int> > puNInt( new vector<int> );
    
  ///pileup 
  Handle<std::vector< PileupSummaryInfo > >  PupInfo;
  iEvent.getByLabel(m_PileupSrc, PupInfo);
  std::vector<PileupSummaryInfo>::const_iterator PVI;
  for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
    //std::cout << " Pileup Information: bunchXing, nInt, TrueNInt " << PVI->getBunchCrossing() << " " << PVI->getPU_NumInteractions() << " "<< PVI->getTrueNumInteractions() <<endl;
    *puNtrueInt = PVI->getTrueNumInteractions(); 
    puBX->push_back(  PVI->getBunchCrossing() ); 
    puNInt->push_back( PVI->getPU_NumInteractions() );
  }
  
  iEvent.put(puBX,"puBX");
  iEvent.put(puNInt,"puNInt");
  iEvent.put(puNtrueInt,"puNtrueInt");


}

// ------------ method called once each job just after ending the event loop  ------------


#include "FWCore/Framework/interface/MakerMacros.h"


DEFINE_FWK_MODULE(EventUserData);
