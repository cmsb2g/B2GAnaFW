#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDMException.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/JetReco/interface/Jet.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/Ptr.h"


#include<vector>
#include <TMath.h>


using namespace reco;
using namespace edm;
using namespace std;

typedef std::vector<pat::Jet> PatJetCollection;

class QGLUserData : public edm::EDProducer {
public:
  QGLUserData( const edm::ParameterSet & );   
  
private:
  void produce( edm::Event &, const edm::EventSetup & );
  edm::EDGetTokenT<std::vector<pat::Jet> >     jetToken_;
  //  edm::EDGetTokenT<edm::ValueMap<float> >     qgToken_;
  
  InputTag jLabel_, qgToken_; 
  
};

QGLUserData::QGLUserData(const edm::ParameterSet& iConfig) :
  jLabel_             (iConfig.getParameter<edm::InputTag>("jetLabel")),
  qgToken_             (iConfig.getParameter<edm::InputTag>("qgtagger"))
{
   produces<vector<pat::Jet> >();
}
  
  
void QGLUserData::produce( edm::Event& iEvent, const edm::EventSetup& iSetup) {


  edm::Handle<std::vector<pat::Jet> > jetHandle;
  edm::Handle<edm::View<pat::Jet> > jets;
  iEvent.getByLabel(jLabel_, jetHandle);
  iEvent.getByLabel(jLabel_, jets);
  auto_ptr<std::vector<pat::Jet> > jetColl( new std::vector<pat::Jet> (*jetHandle) );

  edm::Handle<edm::ValueMap<float>> qgHandle; 
  iEvent.getByLabel(qgToken_, qgHandle);

  for (size_t i = 0; i< jetColl->size(); i++){
    pat::Jet & jet = (*jetColl)[i];
    edm::RefToBase<pat::Jet> jetRef(edm::Ref<std::vector<pat::Jet> >(jetHandle,i));
    float qgLikelihood = 0.;
    if(qgHandle.isValid()) qgLikelihood = (*qgHandle)[jetRef];
    //std::cout<<"QGL: "<<qgLikelihood<<std::endl; 

    jet.addUserFloat("QGL", qgLikelihood);

  }


  iEvent.put( jetColl );

}



#include "FWCore/Framework/interface/MakerMacros.h"


DEFINE_FWK_MODULE(QGLUserData);
