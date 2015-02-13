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
  edm::InputTag m_PileupSrc;
  edm::InputTag m_pvSrc;
 };


EventUserData::EventUserData(const edm::ParameterSet& iConfig)
 {
   m_PileupSrc = iConfig.getParameter<edm::InputTag>("pileup");
   m_pvSrc = iConfig.getParameter<edm::InputTag>("pvSrc");
   
   produces<std::vector<int> >("puBX");
   produces<std::vector<int> >("puNInt");
   produces<int>("puNtrueInt");
   produces<int>("npv");
   produces<double>("vx");
   produces<double>("vy");
   produces<double>("vz"); 
 }


void EventUserData::produce( edm::Event& iEvent, const edm::EventSetup& iSetup) {

  ///pileup   
  std::auto_ptr<int> puNtrueInt (new int);
  std::auto_ptr<std::vector<int> > puBX(  new std::vector<int> );
  std::auto_ptr<std::vector<int> > puNInt( new std::vector<int> );
    

  if ( ! iEvent.eventAuxiliary().isRealData() ) {
    edm::Handle<std::vector< PileupSummaryInfo > >  PupInfo;
    iEvent.getByLabel(m_PileupSrc, PupInfo);
    std::vector<PileupSummaryInfo>::const_iterator PVI;
    for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
      //std::cout << " Pileup Information: bunchXing, nInt, TrueNInt " << PVI->getBunchCrossing() << " " << PVI->getPU_NumInteractions() << " "<< PVI->getTrueNumInteractions() <<endl;
      *puNtrueInt = PVI->getTrueNumInteractions(); 
      puBX->push_back(  PVI->getBunchCrossing() ); 
      puNInt->push_back( PVI->getPU_NumInteractions() );
    }
  }

 

  // primary vertices

  std::auto_ptr<int> npv (new int() );
  std::auto_ptr<double> vx (new double() );
  std::auto_ptr<double> vy (new double() );
  std::auto_ptr<double> vz (new double() );

  edm::Handle<std::vector<reco::Vertex> > h_vtx;
  iEvent.getByLabel( m_pvSrc, h_vtx );
  *npv = h_vtx->size();
  if ( h_vtx->size() > 0 ) {
    *vx = h_vtx->front().x();
    *vy = h_vtx->front().y();
    *vz = h_vtx->front().z();
  }

 
  iEvent.put(puBX,"puBX");
  iEvent.put(puNInt,"puNInt");
  iEvent.put(puNtrueInt,"puNtrueInt");
  iEvent.put( npv, "npv"); 
  iEvent.put( vx, "vx"); 
  iEvent.put( vy, "vy"); 
  iEvent.put( vz, "vz"); 

}

// ------------ method called once each job just after ending the event loop  ------------


#include "FWCore/Framework/interface/MakerMacros.h"


DEFINE_FWK_MODULE(EventUserData);
