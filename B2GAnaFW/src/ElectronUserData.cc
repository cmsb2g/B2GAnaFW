#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDMException.h"
// dR and dPhi                                                                                                                                                              
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
// Vertex
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

// trigger
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h" // gives access to the (release cycle dependent) trigger object codes

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/Ptr.h"

#include<vector>

using namespace reco;
using namespace edm;
using namespace std;


class  ElectronUserData : public edm::EDProducer {
public:
  ElectronUserData( const edm::ParameterSet & );   

private:
  void produce( edm::Event &, const edm::EventSetup & );
  bool isMatchedWithTrigger(const pat::Electron, trigger::TriggerObjectCollection,int&);
  bool passIDWP(string, bool, float, float, float, float, float, float, float, float, bool, int);


  InputTag eleLabel_, pvLabel_, convLabel_;
  InputTag triggerResultsLabel_, triggerSummaryLabel_;
  InputTag hltElectronFilterLabel_;
  TString hltPath_;
  HLTConfigProvider hltConfig;
  int triggerBit;
  int debug_; 
  // ----------member data ---------------------------

  edm::EDGetTokenT<edm::ValueMap<bool> > electronVetoIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > electronTightIdMapToken_;


  std::vector<Float_t> pt_;
  std::vector<Float_t> etaSC_;
  std::vector<Float_t> phiSC_;

  std::vector<Int_t>   passVetoId_;     
  std::vector<Int_t>   passTightId_;     

 };


ElectronUserData::ElectronUserData(const edm::ParameterSet& iConfig):
   eleLabel_(iConfig.getParameter<edm::InputTag>("eleLabel")),
   pvLabel_(iConfig.getParameter<edm::InputTag>("pv")),   // "offlinePrimaryVertex"
   convLabel_(iConfig.getParameter<edm::InputTag>("conversion")),   // "offlinePrimaryVertex"
   triggerResultsLabel_(iConfig.getParameter<edm::InputTag>("triggerResults")),
   triggerSummaryLabel_(iConfig.getParameter<edm::InputTag>("triggerSummary")),
   hltElectronFilterLabel_ (iConfig.getParameter<edm::InputTag>("hltElectronFilter")),   //trigger objects we want to match
   hltPath_ (iConfig.getParameter<std::string>("hltPath")),
   electronVetoIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronVetoIdMap"))),
   electronTightIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronTightIdMap")))
{
  debug_ = iConfig.getUntrackedParameter<int>("debugLevel",int(0));
  
  produces<std::vector<pat::Electron> >();
}

void ElectronUserData::produce( edm::Event& iEvent, const edm::EventSetup& iSetup) {

    
  //PV
  edm::Handle<std::vector<reco::Vertex> > vertices;
  iEvent.getByLabel(pvLabel_, vertices);



  if(debug_>=1) cout<<"vtx size " << vertices->size()<<endl; 

  reco::TrackBase::Point vtxPoint(0,0, 0);
  if(  vertices->size() >= 1 ) {
    vtxPoint = vertices->at(0).position();
  }

  if(debug_>=1) cout<<"vtxPoint " <<vtxPoint.X()<<" "<< vtxPoint.Y()<<" "<< vtxPoint.Z()<<endl; 
    
  //Electrons
  edm::Handle<std::vector<pat::Electron> > eleHandle;
  iEvent.getByLabel(eleLabel_, eleHandle);
  
  auto_ptr<vector<pat::Electron> > eleColl( new vector<pat::Electron> (*eleHandle) );



  Handle<reco::BeamSpot> bsHandle;
  iEvent.getByLabel("offlineBeamSpot", bsHandle);
  const reco::BeamSpot &beamspot = *bsHandle.product();
    
  Handle<reco::ConversionCollection> conversions;
  iEvent.getByLabel(convLabel_, conversions);
  //  iEvent.getByLabel("reducedEgamma:reducedConversions", conversions);


  // Electron ID
  // Get the electron ID data from the event stream.
  edm::Handle<edm::ValueMap<bool> > veto_id_decisions;
  edm::Handle<edm::ValueMap<bool> > tight_id_decisions;
  iEvent.getByToken(electronVetoIdMapToken_,veto_id_decisions);
  iEvent.getByToken(electronTightIdMapToken_,tight_id_decisions);
  passVetoId_.clear();     
  passTightId_.clear();   
  
  // TRIGGER
  bool changedConfig = false;
  if (!hltConfig.init(iEvent.getRun(), iSetup, "HLT", changedConfig)) {
    cout << "Initialization of HLTConfigProvider failed!!" << endl;
    return;
  }

  
  if (changedConfig){
    std::cout << "the curent menu is " << hltConfig.tableName() << " "<< hltConfig.triggerNames().size() <<endl; 
    triggerBit = -1;
    for (size_t j = 0; j < hltConfig.triggerNames().size(); j++) {
      if(debug_ ==99) cout<<"trigName " <<j<<" "<< hltConfig.triggerNames()[j] <<endl; 
      if (TString(hltConfig.triggerNames()[j]).Contains(hltPath_)) triggerBit = j;
    }
    if (triggerBit == -1) cout << "HLT path not found" << endl;
  }
  edm::Handle<edm::TriggerResults> triggerResults;
  iEvent.getByLabel(triggerResultsLabel_, triggerResults);
  if (size_t(triggerBit) < triggerResults->size() )
    if (triggerResults->accept(triggerBit))
      std::cout << "event pass : " << hltPath_ << std::endl;
  // if (triggerResults->accept(triggerBit)){ //crashing here to be understood!
  //   std::cout << "event pass : " << hltPath_ << std::endl;
  // }

  edm::Handle<trigger::TriggerEvent> triggerSummary;
  iEvent.getByLabel(triggerSummaryLabel_, triggerSummary);

  //find the index corresponding to the event
  if(false){ ///to be done after understanding which trigger(s) to use.
  size_t ElectronFilterIndex = (*triggerSummary).filterIndex(hltElectronFilterLabel_);
  trigger::TriggerObjectCollection allTriggerObjects = triggerSummary->getObjects();
  trigger::TriggerObjectCollection ElectronLegObjects;
  if (ElectronFilterIndex < (*triggerSummary).sizeFilters()) { //check if the trigger object is present
    //save the trigger objects corresponding to ele leg
    const trigger::Keys &keysElectrons = (*triggerSummary).filterKeys(ElectronFilterIndex);
    for (size_t j = 0; j < keysElectrons.size(); j++) {
      trigger::TriggerObject foundObject = (allTriggerObjects)[keysElectrons[j]];
      ElectronLegObjects.push_back(foundObject);
    }
  }
  if( debug_ >=1)  std::cout << "ElectronLegObjects: " << ElectronLegObjects.size() << std::endl;
  /////////  /////////  /////////  /////////  /////////  /////////  /////////  /////////  /////////  /////////
  }


  for (size_t i = 0; i< eleColl->size(); i++){
    pat::Electron & el = (*eleColl)[i];

    // Isolation
    GsfElectron::PflowIsolationVariables pfIso = el.pfIsolationVariables();
    // Compute isolation with delta beta correction for PU

    bool isEB = el.isEB() ? true : false;
    //float ptEle  = el.pt();
    //float etaEle = el.superCluster()->eta();
    float dEtaIn = el.deltaEtaSuperClusterTrackAtVtx();
    float dPhiIn = el.deltaPhiSuperClusterTrackAtVtx();
    float full5x5 = el.full5x5_sigmaIetaIeta();
    float hoe = el.hadronicOverEm();
    float absiso = pfIso.sumChargedHadronPt + max(0.0 , pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - 0.5 * pfIso.sumPUPt );
    float relIsoWithDBeta_ = absiso/el.pt();
    float ooEmooP_; 
    if( el.ecalEnergy() == 0 ){
      printf("Electron energy is zero!\n");
      ooEmooP_ = 999;
    }else if( !std::isfinite(el.ecalEnergy())){
      printf("Electron energy is not finite!\n");
      ooEmooP_ = 998;
    }else{
      ooEmooP_ = fabs(1.0/el.ecalEnergy() - el.eSuperClusterOverP()/el.ecalEnergy() );
    }
    // Impact parameter
    float d0 = (-1) * el.gsfTrack()->dxy(vtxPoint);
    float dz = el.gsfTrack()->dz(vtxPoint);


    if(debug_>=1) cout<<" ele " << i <<" pt "<< el.pt()<<" eta "<<el.eta()<<"fabs(1/E-1/P) "<< ooEmooP_ <<" d0 "<< d0 <<" dz " << dz <<" iso " << relIsoWithDBeta_<<endl; 		    
    
    float missHits = el.gsfTrack()->hitPattern().numberOfLostTrackerHits(HitPattern::MISSING_INNER_HITS);
    // conversion rejection match
    bool hasMatchConv = ConversionTools::hasMatchedConversion(el, conversions, beamspot.position());

    bool isVeto = passIDWP("VETO",isEB, dEtaIn, dPhiIn, full5x5, hoe, d0, dz, ooEmooP_, relIsoWithDBeta_, hasMatchConv, missHits);
    bool isLoose = passIDWP("LOOSE",isEB, dEtaIn, dPhiIn, full5x5, hoe, d0, dz, ooEmooP_, relIsoWithDBeta_, hasMatchConv, missHits);
    bool isMedium = passIDWP("MEDIUM",isEB, dEtaIn, dPhiIn, full5x5, hoe, d0, dz, ooEmooP_, relIsoWithDBeta_, hasMatchConv, missHits);
    bool isTight = passIDWP("TIGHT",isEB, dEtaIn, dPhiIn, full5x5, hoe, d0, dz, ooEmooP_, relIsoWithDBeta_, hasMatchConv, missHits);
    // Look up the ID decision for this electron in 
    // the ValueMap object and store it. We need a Ptr object as the key.
    //    const Ptr<pat::Electron> elPtr(eleHandle, i);
    //const Ptr<pat::Electron> elPtr(eleHandle, i);
    //std::cout<<"ValueMap: "<<veto_id_decisions->contains(elPtr.id())<<std::endl;
    //bool isPassVeto  = (*veto_id_decisions)[ elPtr ];
    //bool isPassTight = (*tight_id_decisions)[ elPtr ];
    //passVetoId_.push_back( isPassVeto );
    //passTightId_.push_back( isPassTight );

    el.addUserFloat("d0",          d0);
    el.addUserFloat("dz",          dz);
    el.addUserFloat("iso03",       relIsoWithDBeta_);
    el.addUserFloat("ooEmooP",     ooEmooP_);
    el.addUserFloat("missHits",     missHits);
    el.addUserFloat("hasMatchConv",     hasMatchConv);
    el.addUserFloat("isVeto",     isVeto);
    el.addUserFloat("isLoose",    isLoose);
    el.addUserFloat("isMedium",   isMedium);
    el.addUserFloat("isTight",    isTight);
    //el.addUserFloat("isPassVeto",     isPassVeto);
    //    el.addUserFloat("isPassTight",     isPassTight);

        

  }

 iEvent.put( eleColl );

}

// ------------ method called once each job just after ending the event loop  ------------

bool
ElectronUserData::isMatchedWithTrigger(const pat::Electron p, trigger::TriggerObjectCollection triggerObjects, int& index)
{
  for (size_t i = 0 ; i < triggerObjects.size() ; i++){
    float deltaR = sqrt(pow(triggerObjects[i].eta()-p.eta(),2)+ pow(acos(cos(triggerObjects[i].phi()-p.phi())),2)) ;
    std::cout << "deltaR: " << deltaR << std::endl;
    if (deltaR<0.1) {
      index = i;
      return true;
    }
  }
  return false;
}

bool ElectronUserData::passIDWP(string WP, bool isEB, float dEtaIn, float dPhiIn, float full5x5, float hoe, float d0, float dz, float ooemoop, float reliso, bool conv, int missHits){
  bool pass = false;

  if(WP == "VETO"){
    if(isEB){
      pass = (fabs(dEtaIn) < 0.02) && (fabs(dPhiIn) < 0.2579) && (full5x5 < 0.0125) && (hoe < 0.2564) && (fabs(d0) < 0.025) && (fabs(dz) < 0.5863) && (fabs(ooemoop) < 0.1508) && (reliso < 0.3313) && conv && (missHits < 3);
    }
    else{
      pass = (fabs(dEtaIn) < 0.0141) && (fabs(dPhiIn) < 0.2591) && (full5x5 < 0.0371) && (hoe < 0.1335) && (fabs(d0) < 0.2232) && (fabs(dz) < 0.9513) && (fabs(ooemoop) < 0.1542) && (reliso < 0.3816) && conv && (missHits < 4);
    }
  }
  if(WP == "LOOSE"){
    if(isEB){
      pass = (fabs(dEtaIn) < 0.0181) && (fabs(dPhiIn) < 0.0936) && (full5x5 < 0.0123) && (hoe < 0.141) && (fabs(d0) < 0.0166) && (fabs(dz) < 0.54342) && (fabs(ooemoop) < 0.1353) && (reliso < 0.24) && conv && (missHits < 2);
    }
    else{
      pass = (fabs(dEtaIn) < 0.0124) && (fabs(dPhiIn) < 0.0642) && (full5x5 < 0.035) && (hoe < 0.1115) && (fabs(d0) < 0.098) && (fabs(dz) < 0.9187) && (fabs(ooemoop) < 0.1443) && (reliso < 0.3529) && conv && (missHits < 2);
    }
      }

  if(WP == "MEDIUM"){
    if(isEB){
      pass = (fabs(dEtaIn) < 0.0106) && (fabs(dPhiIn) < 0.0323) && (full5x5 < 0.0107) && (hoe < 0.067) && (fabs(d0) < 0.0131) && (fabs(dz) < 0.22310) && (fabs(ooemoop) < 0.1043) && (reliso < 0.2179) && conv && (missHits < 2);
    }
    else{
      pass = (fabs(dEtaIn) < 0.0108) && (fabs(dPhiIn) < 0.0455) && (full5x5 < 0.0318) && (hoe < 0.097) && (fabs(d0) < 0.0845) && (fabs(dz) < 0.7523) && (fabs(ooemoop) < 0.1201) && (reliso < 0.254) && (missHits < 2);
    }
      }

  if(WP == "TIGHT"){
    if(isEB){
      pass = (fabs(dEtaIn) < 0.0091) && (fabs(dPhiIn) < 0.031) && (full5x5 < 0.0106) && (hoe < 0.0532) && (fabs(d0) < 0.0126) && (fabs(dz) < 0.0116) && (fabs(ooemoop) < 0.0609) && (reliso < 0.1649) && conv && (missHits < 2);
    }
    else{
      pass = (fabs(dEtaIn) < 0.0106) && (fabs(dPhiIn) < 0.0359) && (full5x5 < 0.0305) && (hoe < 0.0835) && (fabs(d0) < 0.0163) && (fabs(dz) < 0.5999) && (fabs(ooemoop) < 0.1126) && (reliso < 0.2075) && conv && (missHits < 2);
    }
      }
  return pass;
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(ElectronUserData);
