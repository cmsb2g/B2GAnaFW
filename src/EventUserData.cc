#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDMException.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include <Math/VectorUtil.h>


class  EventUserData : public edm::EDProducer {
public:
  EventUserData( const edm::ParameterSet & );   

private:
  void produce( edm::Event &, const edm::EventSetup & );
  //InputTag LHELabel_;
  edm::EDGetTokenT< std::vector< PileupSummaryInfo > > m_PileupSrc;
  edm::EDGetTokenT< std::vector< reco::Vertex > > m_pvSrc;
 };


EventUserData::EventUserData(const edm::ParameterSet& iConfig):
   m_PileupSrc(consumes<std::vector<PileupSummaryInfo>>(iConfig.getParameter<edm::InputTag>("pileup"))),
   m_pvSrc(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("pvSrc")))
 {
   
   produces<std::vector<int> >("puBX");
   produces<std::vector<int> >("puNInt");
   produces<int>("puNtrueInt");
 }


void EventUserData::produce( edm::Event& iEvent, const edm::EventSetup& iSetup) {

  ///pileup   
  std::unique_ptr<int> puNtrueInt (new int);
  std::unique_ptr<std::vector<int> > puBX(  new std::vector<int> );
  std::unique_ptr<std::vector<int> > puNInt( new std::vector<int> );

  if ( ! iEvent.eventAuxiliary().isRealData() ) {
    edm::Handle<std::vector< PileupSummaryInfo > >  PupInfo;
    iEvent.getByToken(m_PileupSrc, PupInfo);
    std::vector<PileupSummaryInfo>::const_iterator PVI;
    for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
      //std::cout << " Pileup Information: bunchXing, nInt, TrueNInt " << PVI->getBunchCrossing() << " " << PVI->getPU_NumInteractions() << " "<< PVI->getTrueNumInteractions() <<endl;
      *puNtrueInt = PVI->getTrueNumInteractions(); 
      puBX->push_back(  PVI->getBunchCrossing() ); 
      puNInt->push_back( PVI->getPU_NumInteractions() );
    }
  }

  iEvent.put(std::move(puBX),"puBX");
  iEvent.put(std::move(puNInt),"puNInt");
  iEvent.put(std::move(puNtrueInt),"puNtrueInt");
}

// ------------ method called once each job just after ending the event loop  ------------


#include "FWCore/Framework/interface/MakerMacros.h"


DEFINE_FWK_MODULE(EventUserData);
