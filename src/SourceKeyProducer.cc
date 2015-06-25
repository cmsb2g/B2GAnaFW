#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDMException.h"

#include "DataFormats/Candidate/interface/Candidate.h"

// trigger
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h" // gives access to the (release cycle dependent) trigger object codes


#include <vector>

class SourceKeyProducer : public edm::EDProducer {
  public:
  typedef std::vector< std::vector<int> >     index_collection; 
    SourceKeyProducer( const edm::ParameterSet & );   

  private:
    void produce( edm::Event &, const edm::EventSetup & );
  edm::InputTag jLabel_; 
};


SourceKeyProducer::SourceKeyProducer(const edm::ParameterSet& iConfig) :
  jLabel_             (iConfig.getParameter<edm::InputTag>("srcLabel"))
{
  produces< index_collection >();
}


void SourceKeyProducer::produce( edm::Event& iEvent, const edm::EventSetup& iSetup) {

  edm::Handle<edm::View<reco::Candidate> > candHandle;
  iEvent.getByLabel(jLabel_, candHandle);
  std::auto_ptr< index_collection > keys( new index_collection () );

  for ( auto const & cand : *candHandle ){

    //// Cand constituent indices for lepton matching
    std::vector<int> constituentIndices;
    for ( unsigned int isrc = 0; isrc < cand.numberOfSourceCandidatePtrs(); ++isrc ){
      constituentIndices.push_back( cand.sourceCandidatePtr(isrc).key() );
    }

    keys->push_back( constituentIndices );


  } //// Loop over all cands 

  iEvent.put( keys );

}


#include "FWCore/Framework/interface/MakerMacros.h"


DEFINE_FWK_MODULE(SourceKeyProducer);
