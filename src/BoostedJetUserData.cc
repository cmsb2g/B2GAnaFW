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

class BoostedJetUserData : public edm::EDProducer {
  public:
    BoostedJetUserData( const edm::ParameterSet & );   

  private:
    void produce( edm::Event &, const edm::EventSetup & );


    edm::EDGetTokenT<std::vector<pat::Jet> >     jToken_;
};


BoostedJetUserData::BoostedJetUserData(const edm::ParameterSet& iConfig) :
  jToken_             (consumes<std::vector<pat::Jet> > ( iConfig.getParameter<edm::InputTag>("jetLabel") ))
{
  produces<vector<pat::Jet> >();
}


void BoostedJetUserData::produce( edm::Event& iEvent, const edm::EventSetup& iSetup) {

  edm::Handle<std::vector<pat::Jet> > jetHandle;
  iEvent.getByToken(jToken_, jetHandle);


  auto_ptr<vector<pat::Jet> > jetColl( new vector<pat::Jet> (*jetHandle) );


  for (size_t i = 0; i< jetHandle->size(); i++){
    pat::Jet & jet = (*jetColl)[i];

    // Get top tagged subjets and their keys
    auto const & topSubjets = jet.subjets("CMSTopTag");
    int topSubjet0=-1, topSubjet1=-1, topSubjet2=-1, topSubjet3=-1;
    if ( topSubjets.size() > 0 )
      topSubjet0 = topSubjets[0].key();
    if ( topSubjets.size() > 1 )
      topSubjet1 = topSubjets[1].key();
    if ( topSubjets.size() > 2 )
      topSubjet2 = topSubjets[2].key();
    if ( topSubjets.size() > 3 )
      topSubjet3 = topSubjets[3].key();
    jet.addUserInt("TopSubjet0", topSubjet0  );
    jet.addUserInt("TopSubjet1", topSubjet1  );
    jet.addUserInt("TopSubjet2", topSubjet2  );
    jet.addUserInt("TopSubjet3", topSubjet3  );	

    // Get V-tagged subjets and their keys
    auto const & vSubjets = jet.subjets("SoftDrop");
    
    int vSubjet0=-1, vSubjet1=-1;
    if ( vSubjets.size() > 0 )
      vSubjet0 = vSubjets[0].key();
    if ( vSubjets.size() > 1 )
      vSubjet1 = vSubjets[1].key();
    jet.addUserInt("VSubjet0", vSubjet0  );
    jet.addUserInt("VSubjet1", vSubjet1  );

  } //// Loop over all jets 

  iEvent.put( jetColl );

}


#include "FWCore/Framework/interface/MakerMacros.h"


DEFINE_FWK_MODULE(BoostedJetUserData);
