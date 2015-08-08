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

// MiniIsolation
#include "Isolations.h"

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
#include <TMath.h>

using namespace reco;
using namespace edm;
using namespace std;


class  ElectronUserData : public edm::EDProducer {
public:
  ElectronUserData( const edm::ParameterSet & );   

private:
  void produce( edm::Event &, const edm::EventSetup & );
  float getEA(float);
  bool isMatchedWithTrigger(const pat::Electron, trigger::TriggerObjectCollection,int&);
  bool passIDWP(string, bool, float, float, float, float, float, float, float, float, bool, int);


  InputTag eleLabel_, pvLabel_, convLabel_, rho_;
  InputTag triggerResultsLabel_, triggerSummaryLabel_;
  InputTag hltElectronFilterLabel_;
  TString hltPath_;
  HLTConfigProvider hltConfig;
  int triggerBit;
  int debug_; 
  // ----------member data ---------------------------

  // edm::EDGetTokenT<edm::ValueMap<bool> > electronVetoIdMapToken_;
  // edm::EDGetTokenT<edm::ValueMap<bool> > electronTightIdMapToken_;


  //  std::vector<Float_t> pt_;
  //  std::vector<Float_t> etaSC_;
  //  std::vector<Float_t> phiSC_;

  //  std::vector<Int_t>   passVetoId_;     
  //  std::vector<Int_t>   passTightId_;     

 };


ElectronUserData::ElectronUserData(const edm::ParameterSet& iConfig):
   eleLabel_(iConfig.getParameter<edm::InputTag>("eleLabel")),
   pvLabel_(iConfig.getParameter<edm::InputTag>("pv")),   // "offlinePrimaryVertex"
   convLabel_(iConfig.getParameter<edm::InputTag>("conversion")),   // "offlinePrimaryVertex"
   rho_( iConfig.getParameter<edm::InputTag>("rho")),
   triggerResultsLabel_(iConfig.getParameter<edm::InputTag>("triggerResults")),
   triggerSummaryLabel_(iConfig.getParameter<edm::InputTag>("triggerSummary")),
   hltElectronFilterLabel_ (iConfig.getParameter<edm::InputTag>("hltElectronFilter")),   //trigger objects we want to match
   hltPath_ (iConfig.getParameter<std::string>("hltPath"))
   //electronVetoIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronVetoIdMap"))),
   //electronTightIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronTightIdMap")))
{
  debug_ = iConfig.getUntrackedParameter<int>("debugLevel",int(0));
  
  produces<std::vector<pat::Electron> >();
}

void ElectronUserData::produce( edm::Event& iEvent, const edm::EventSetup& iSetup) {

    
  //PV
  edm::Handle<std::vector<reco::Vertex> > vertices;
  iEvent.getByLabel(pvLabel_, vertices);

  //RHO
  edm::Handle<double> rhoHandle;
  iEvent.getByLabel(rho_,rhoHandle);
  double rho = *rhoHandle; 

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

  //PackedPFCands for Mini-iso
  edm::Handle<pat::PackedCandidateCollection> packedPFCands;
  iEvent.getByLabel("packedPFCandidates", packedPFCands);


  Handle<reco::BeamSpot> bsHandle;
  iEvent.getByLabel("offlineBeamSpot", bsHandle);
  const reco::BeamSpot &beamspot = *bsHandle.product();
    
  Handle<reco::ConversionCollection> conversions;
  iEvent.getByLabel(convLabel_, conversions);
  //  iEvent.getByLabel("reducedEgamma:reducedConversions", conversions);


  // Electron ID
  // Get the electron ID data from the event stream.
  //edm::Handle<edm::ValueMap<bool> > veto_id_decisions;
  //edm::Handle<edm::ValueMap<bool> > tight_id_decisions;
  //iEvent.getByToken(electronVetoIdMapToken_,veto_id_decisions);
  //iEvent.getByToken(electronTightIdMapToken_,tight_id_decisions);
  //passVetoId_.clear();     
  //passTightId_.clear();   
  
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
    double miniIso = getPFMiniIsolation(packedPFCands, dynamic_cast<const reco::Candidate *>(&el), 0.05, 0.2, 10., false);
    double EA = getEA(el.eta());
    //double rho = 1.0;
    float absiso_EA = pfIso.sumChargedHadronPt + max(0.0 , pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - rho * EA );
    float relIsoWithEA_ = absiso_EA/el.pt();
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
    //const Ptr<pat::Electron> elPtr(eleHandle, i);
    //std::cout<<"ValueMap: "<<veto_id_decisions->contains(elPtr.id())<<std::endl;
    //bool isPassVeto  = (*veto_id_decisions)[ elPtr ];
    //bool isPassTight = (*tight_id_decisions)[ elPtr ];
    //passVetoId_.push_back( isPassVeto );
    //passTightId_.push_back( isPassTight );

    el.addUserFloat("dEtaIn",      dEtaIn);
    el.addUserFloat("dPhiIn",      dPhiIn);
    el.addUserFloat("full5x5",     full5x5);
    el.addUserFloat("hoe",         hoe);
    el.addUserFloat("d0",          d0);
    el.addUserFloat("dz",          dz);
    el.addUserFloat("ooEmooP",     ooEmooP_);
    el.addUserFloat("missHits",     missHits);
    el.addUserFloat("hasMatchConv",     hasMatchConv);   
    el.addUserFloat("iso03db",    relIsoWithDBeta_);
    el.addUserFloat("iso03",       relIsoWithEA_);
    el.addUserFloat("miniIso",     miniIso);
    el.addUserFloat("sumChargedHadronPt",   pfIso.sumChargedHadronPt);
    el.addUserFloat("sumNeutralHadronEt", pfIso.sumNeutralHadronEt  );
    el.addUserFloat("sumPhotonEt",  pfIso.sumPhotonEt );
    el.addUserFloat("sumPUPt", pfIso.sumPUPt  );


    el.addUserFloat("rho", rho );
    el.addUserFloat("EA", EA );
    el.addUserFloat("isVeto",     isVeto);
    el.addUserFloat("isLoose",    isLoose);
    el.addUserFloat("isMedium",   isMedium);
    el.addUserFloat("isTight",    isTight);




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


float ElectronUserData::getEA(float eta)
{
  // The following values refer to EA for cone 0.3 and fixedGridRhoFastjetAll. 
  // They are valid for electrons only, different EA are available for muons.
  float effArea = 0.;
  if(abs(eta)>0.0 && abs(eta)<=0.8) effArea = 0.1013;
  if(abs(eta)>0.8 && abs(eta)<=1.3) effArea = 0.0988;
  if(abs(eta)>1.3 && abs(eta)<=2.0) effArea = 0.0572;
  if(abs(eta)>2.0 && abs(eta)<=2.2) effArea = 0.0842;
  if(abs(eta)>2.2 && abs(eta)<=2.5) effArea = 0.1530;
  return effArea;
}

bool ElectronUserData::passIDWP(string WP, bool isEB, float dEtaIn, float dPhiIn, float full5x5, float hoe, float d0, float dz, float ooemoop, float reliso, bool conv, int missHits){
  bool pass = false;

  if(WP == "VETO"){
    if(isEB){
      pass = (fabs(dEtaIn) <  0.013625 ) && (fabs(dPhiIn) <  0.230374 ) && (full5x5 < 0.011586 ) && (hoe <  0.181130 ) && (fabs(d0) < 0.094095 ) && (fabs(dz) <  0.713070 ) && (fabs(ooemoop) <  0.295751 ) && (reliso < 0.158721 ) && !conv && (missHits <= 2);
    }
    else{
      pass = (fabs(dEtaIn) <  0.011932 ) && (fabs(dPhiIn) <  0.255450 ) && (full5x5 < 0.031849 ) && (hoe <  0.223870 ) && (fabs(d0) < 0.342293 ) && (fabs(dz) < 0.953461 ) && (fabs(ooemoop) < 0.155501 ) && (reliso < 0.177032 ) && !conv && (missHits <= 3);
    }
  }
  if(WP == "LOOSE"){
    if(isEB){
      pass = (fabs(dEtaIn) < 0.009277 ) && (fabs(dPhiIn) < 0.094739 ) && (full5x5 <  0.010331 ) && (hoe < 0.093068 ) && (fabs(d0) < 0.035904 ) && (fabs(dz) < 0.075496 ) && (fabs(ooemoop) <  0.189968 ) && (reliso < 0.130136 ) && !conv && (missHits <= 1);
    }
    else{
      pass = (fabs(dEtaIn) < 0.009833 ) && (fabs(dPhiIn) < 0.149934 ) && (full5x5 < 0.031838 ) && (hoe < 0.115754 ) && (fabs(d0) < 0.099266 ) && (fabs(dz) < 0.197897 ) && (fabs(ooemoop) < 0.140662 ) && (reliso < 0.163368 ) && !conv && (missHits <= 1);
    }
  }

  if(WP == "MEDIUM"){
    if(isEB){
      pass = (fabs(dEtaIn) <  0.008925 ) && (fabs(dPhiIn) <  0.035973 ) && (full5x5 <  0.009996 ) && (hoe <  0.050537 ) && (fabs(d0) <  0.012235 ) && (fabs(dz) <  0.042020 ) && (fabs(ooemoop) <  0.091942 ) && (reliso <  0.107587 ) && !conv && (missHits <= 1);
    }
    else{
      pass = (fabs(dEtaIn) <  0.007429 ) && (fabs(dPhiIn) <  0.067879 ) && (full5x5 <  0.030135 ) && (hoe <  0.086782 ) && (fabs(d0) <  0.036719 ) && (fabs(dz) <  0.138142 ) && (fabs(ooemoop) <  0.100683 ) && (reliso <  0.113254 ) && !conv && (missHits <= 1);
    }
  }

  if(WP == "TIGHT"){
    if(isEB){
      pass = (fabs(dEtaIn) <  0.006046 ) && (fabs(dPhiIn) <  0.028092 ) && (full5x5 <  0.009947 ) && (hoe <  0.045772 ) && (fabs(d0) <  0.008790 ) && (fabs(dz) <  0.021226 ) && (fabs(ooemoop) <  0.020118 ) && (reliso <  0.069537 ) && !conv && (missHits <= 1);
    }
    else{
      pass = (fabs(dEtaIn) <  0.007057 ) && (fabs(dPhiIn) <  0.030159 ) && (full5x5 <  0.028237 ) && (hoe <  0.067778 ) && (fabs(d0) <  0.027984 ) && (fabs(dz) <  0.133431 ) && (fabs(ooemoop) <  0.098919 ) && (reliso <  0.078265 ) && !conv && (missHits <= 1);
    }
  }
  return pass;
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(ElectronUserData);
