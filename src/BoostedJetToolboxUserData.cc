#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
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

typedef std::vector<pat::Jet> PatJetCollection;

class BoostedJetToolboxUserData : public edm::EDProducer {
  public:
    BoostedJetToolboxUserData( const edm::ParameterSet & );   

  private:
    void produce( edm::Event &, const edm::EventSetup & );


    edm::EDGetTokenT<std::vector<pat::Jet> >     jToken_;
    edm::EDGetTokenT<std::vector<pat::Jet> >     tToken_;
    edm::EDGetTokenT<std::vector<pat::Jet> >     vToken_;

    double distMax_;

};


BoostedJetToolboxUserData::BoostedJetToolboxUserData(const edm::ParameterSet& iConfig) :
  jToken_             (consumes<std::vector<pat::Jet> > ( iConfig.getParameter<edm::InputTag>("jetLabel") )),
  tToken_             (consumes<std::vector<pat::Jet> > ( iConfig.getParameter<edm::InputTag>("topjetLabel")) ),
  vToken_             (consumes<std::vector<pat::Jet> > ( iConfig.getParameter<edm::InputTag>("vjetLabel")) ),
  distMax_            ( iConfig.getParameter<double>( "distMax" ) )
{
  produces<vector<pat::Jet> >();
}


void BoostedJetToolboxUserData::produce( edm::Event& iEvent, const edm::EventSetup& iSetup) {

  edm::Handle<std::vector<pat::Jet> > jetHandle, topjetHandle, vjetHandle;
  iEvent.getByToken(jToken_, jetHandle);
  iEvent.getByToken(tToken_, topjetHandle);
  iEvent.getByToken(vToken_, vjetHandle);


  auto_ptr<vector<pat::Jet> > jetColl( new vector<pat::Jet> (*jetHandle) );


  for (size_t i = 0; i< jetHandle->size(); i++){
    pat::Jet & jet = (*jetColl)[i];

    for ( auto const & topJet : *topjetHandle ) {
      float temp_dR2 = reco::deltaR2(jet.eta(),jet.phi(),topJet.eta(),topJet.phi());
      if ( temp_dR2 < distMax_ ) {
	int topSubjet0=-1, topSubjet1=-1, topSubjet2=-1, topSubjet3=-1;
	if ( topJet.numberOfDaughters() > 0 )
	  topSubjet0 = topJet.daughterPtr(0).key();
	if ( topJet.numberOfDaughters() > 1 )
	  topSubjet1 = topJet.daughterPtr(1).key();
	if ( topJet.numberOfDaughters() > 2 )
	  topSubjet2 = topJet.daughterPtr(2).key();
	if ( topJet.numberOfDaughters() > 3 )
	  topSubjet3 = topJet.daughterPtr(3).key();
	jet.addUserInt("TopSubjet0", topSubjet0  );
	jet.addUserInt("TopSubjet1", topSubjet1  );
	jet.addUserInt("TopSubjet2", topSubjet2  );
	jet.addUserInt("TopSubjet3", topSubjet3  );	
	break;
      }

    }



    for ( auto const & vJet : *vjetHandle ) {
      float temp_dR2 = reco::deltaR2(jet.eta(),jet.phi(),vJet.eta(),vJet.phi());
      if ( temp_dR2 < distMax_ ) {

	int vSubjet0=-1, vSubjet1=-1;
	if ( vJet.numberOfDaughters() > 0 )
	  vSubjet0 = vJet.daughterPtr(0).key();
	if ( vJet.numberOfDaughters() > 1 )
	  vSubjet1 = vJet.daughterPtr(1).key();
	jet.addUserInt("VSubjet0", vSubjet0  );
	jet.addUserInt("VSubjet1", vSubjet1  );
	break;
      }

    }

  } //// Loop over all jets 

  iEvent.put( jetColl );

}


#include "FWCore/Framework/interface/MakerMacros.h"


DEFINE_FWK_MODULE(BoostedJetToolboxUserData);

