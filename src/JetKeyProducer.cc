#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDMException.h"

#include "DataFormats/JetReco/interface/Jet.h"

// trigger
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h" // gives access to the (release cycle dependent) trigger object codes


#include <vector>

class JetKeyProducer : public edm::EDProducer {
  public:
  typedef std::vector< std::vector<int> >     index_collection; 
    JetKeyProducer( const edm::ParameterSet & );   

  private:
    void produce( edm::Event &, const edm::EventSetup & );
  edm::InputTag jLabel_; 
};


JetKeyProducer::JetKeyProducer(const edm::ParameterSet& iConfig) :
  jLabel_             (iConfig.getParameter<edm::InputTag>("jetLabel"))
{
  produces< index_collection >();
}


void JetKeyProducer::produce( edm::Event& iEvent, const edm::EventSetup& iSetup) {

  edm::Handle<edm::View<reco::Jet> > jetHandle;
  iEvent.getByLabel(jLabel_, jetHandle);
  std::auto_ptr< index_collection > keys( new index_collection () );

  for ( auto const & jet : *jetHandle ){

    //// Jet constituent indices for lepton matching
    std::vector<int> constituentIndices;
    auto constituents = jet.daughterPtrVector();
    for ( auto & constituent : constituents ) {
      constituentIndices.push_back( constituent.key() );
    }

    keys->push_back( constituentIndices );


  } //// Loop over all jets 

  iEvent.put( keys );

}


#include "FWCore/Framework/interface/MakerMacros.h"


DEFINE_FWK_MODULE(JetKeyProducer);
