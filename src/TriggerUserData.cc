#include <memory>
#include <cmath>
#include <typeinfo>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDMException.h"


// dR and dPhi
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"

// Muons
//#include "ttbarDM/TopPlusDMAna/interface/Muons.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonPFIsolation.h"

// Vertex
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

// trigger
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h" // gives access to the (release cycle dependent) trigger object codes

#include <TFile.h>
#include <TH1F.h>
#include <TGraphAsymmErrors.h>
#include<vector>

using namespace reco;
using namespace edm;
using namespace std;
using namespace trigger;

class  TriggerUserData : public edm::EDProducer {
public:
  TriggerUserData( const edm::ParameterSet & );   

private:
  void produce( edm::Event &, const edm::EventSetup & );

  edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
  edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;

  bool storePrescales_ ; 
  TString mainROOTFILEdir_;
  std::string hltProcName_ ; 

  HLTConfigProvider hltConfig;
  int triggerBit;

 };


TriggerUserData::TriggerUserData(const edm::ParameterSet& iConfig):
  triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
  triggerPrescales_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescales"))),
  storePrescales_(iConfig.getUntrackedParameter<bool>("storePrescales")), 
  mainROOTFILEdir_(iConfig.getUntrackedParameter<std::string>("mainROOTFILEdir","")),
  hltProcName_(iConfig.getUntrackedParameter<std::string>("hltProcName"))
{
  produces<std::vector<float>>("triggerBitTree");
  if ( storePrescales_ ) produces<std::vector<int>>("triggerPrescaleTree");  
  produces<std::vector<std::string>>("triggerNameTree");  
 }

void TriggerUserData::produce( edm::Event& iEvent, const edm::EventSetup& iSetup) {
  bool changedConfig = false;

  if (!hltConfig.init(iEvent.getRun(), iSetup, hltProcName_, changedConfig)) {
    std::cout << "Initialization of HLTConfigProvider failed!!" << std::endl;
    return;
  }

  auto_ptr<vector<float> > triggerBitTree( new vector<float> );
  auto_ptr<vector<int> > triggerPrescaleTree( new vector<int> );
  auto_ptr<vector<std::string> > triggerNameTree( new vector<std::string> );
 
  edm::Handle<edm::TriggerResults> triggerBits;
  edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;

  iEvent.getByToken(triggerBits_, triggerBits);
  if ( storePrescales_ ) iEvent.getByToken(triggerPrescales_, triggerPrescales);

  for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {

    triggerBitTree->push_back(triggerBits->accept(i));
    if ( storePrescales_ ) triggerPrescaleTree->push_back(triggerPrescales->getPrescaleForIndex(i));
    //std::cout << "TEST " <<  (hltConfig.triggerNames()[i]) << std::endl;
    //const char * name = TString(hltConfig.triggerNames()[i]);
    triggerNameTree->push_back((hltConfig.triggerNames()[i]));
    //std::cout << typeid((hltConfig.triggerNames()[i])).name() << endl;
  }
  
  iEvent.put( triggerBitTree, "triggerBitTree" );
  if ( storePrescales_ ) iEvent.put( triggerPrescaleTree, "triggerPrescaleTree" );
  iEvent.put( triggerNameTree, "triggerNameTree" );

}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TriggerUserData);
