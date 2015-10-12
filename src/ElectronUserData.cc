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

// VID Debugging 
#include "DataFormats/PatCandidates/interface/VIDCutFlowResult.h"

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
  float getEAOld(float);
  float getEAWithGenWeightOld(float);
  bool isMatchedWithTrigger(const pat::Electron, trigger::TriggerObjectCollection,int&);
  bool passIDWP(string, bool, float, float, float, float, float, float, float, bool, int);

  void printCutFlowResult(vid::CutFlowResult &cutflow);

  InputTag eleLabel_, pvLabel_, packedPFCandsLabel_, convLabel_, rho_;
  InputTag triggerResultsLabel_, triggerSummaryLabel_;
  InputTag hltElectronFilterLabel_;
  TString hltPath_;
  HLTConfigProvider hltConfig;
  int triggerBit;
  int debug_; 
  // ----------member data ---------------------------

  edm::EDGetTokenT<edm::ValueMap<bool> > electronVetoIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > electronLooseIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > electronMediumIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > electronTightIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<vid::CutFlowResult> > electronHEEPIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<vid::CutFlowResult> > eleMediumIdFullInfoMapToken_;

  bool verboseIdFlag_;

  //  std::vector<Float_t> pt_;
  //  std::vector<Float_t> etaSC_;
  //  std::vector<Float_t> phiSC_;

  //  std::vector<Int_t>   passVetoId_;     
  //  std::vector<Int_t>   passTightId_;     

 };


ElectronUserData::ElectronUserData(const edm::ParameterSet& iConfig):
   eleLabel_(iConfig.getParameter<edm::InputTag>("eleLabel")),
   pvLabel_(iConfig.getParameter<edm::InputTag>("pv")),   // "offlinePrimaryVertex"
   packedPFCandsLabel_(iConfig.getParameter<edm::InputTag>("packedPFCands")),
   convLabel_(iConfig.getParameter<edm::InputTag>("conversion")),   // "offlinePrimaryVertex"
   rho_( iConfig.getParameter<edm::InputTag>("rho")),
   triggerResultsLabel_(iConfig.getParameter<edm::InputTag>("triggerResults")),
   triggerSummaryLabel_(iConfig.getParameter<edm::InputTag>("triggerSummary")),
   hltElectronFilterLabel_ (iConfig.getParameter<edm::InputTag>("hltElectronFilter")),   //trigger objects we want to match
   hltPath_ (iConfig.getParameter<std::string>("hltPath")),
   electronVetoIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronVetoIdMap"))),
   electronLooseIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronLooseIdMap"))),
   electronMediumIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronMediumIdMap"))),
   electronTightIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronTightIdMap"))),
   electronHEEPIdMapToken_(consumes<edm::ValueMap<vid::CutFlowResult> >(iConfig.getParameter<edm::InputTag>("electronHEEPIdMap"))), 
   eleMediumIdFullInfoMapToken_(consumes<edm::ValueMap<vid::CutFlowResult> >(iConfig.getParameter<edm::InputTag>("eleMediumIdFullInfoMap"))),
   verboseIdFlag_(iConfig.getParameter<bool>("eleIdVerbose"))
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
  iEvent.getByLabel(packedPFCandsLabel_, packedPFCands);

  Handle<reco::BeamSpot> bsHandle;
  iEvent.getByLabel("offlineBeamSpot", bsHandle);
  const reco::BeamSpot &beamspot = *bsHandle.product();
    
  Handle<reco::ConversionCollection> conversions;
  iEvent.getByLabel(convLabel_, conversions);
  //  iEvent.getByLabel("reducedEgamma:reducedConversions", conversions);


  // Electron ID
  // Get the electron ID data from the event stream.
  edm::Handle<edm::ValueMap<bool> > veto_id_decisions;
  edm::Handle<edm::ValueMap<bool> > loose_id_decisions;
  edm::Handle<edm::ValueMap<bool> > medium_id_decisions;
  edm::Handle<edm::ValueMap<bool> > tight_id_decisions;
  edm::Handle<edm::ValueMap<vid::CutFlowResult> > HEEP_id_cutflow_data;
  edm::Handle<edm::ValueMap<vid::CutFlowResult> > medium_id_cutflow_data;
  iEvent.getByToken(electronVetoIdMapToken_,veto_id_decisions);
  iEvent.getByToken(electronLooseIdMapToken_,loose_id_decisions);
  iEvent.getByToken(electronMediumIdMapToken_,medium_id_decisions);
  iEvent.getByToken(electronTightIdMapToken_,tight_id_decisions);
  iEvent.getByToken(electronHEEPIdMapToken_,HEEP_id_cutflow_data);
  iEvent.getByToken(eleMediumIdFullInfoMapToken_,medium_id_cutflow_data);
  //passVetoId_.clear();     
  //passTightId_.clear();

  // Cuts to be masked for non-iso version. 
  vector<string> maskCuts;
  maskCuts.push_back("GsfEleTrkPtIsoCut_0"); 
  maskCuts.push_back("GsfEleEmHadD1IsoRhoCut_0");   
  
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

    bool isVeto = passIDWP("VETO",isEB, dEtaIn, dPhiIn, full5x5, hoe, d0, dz, ooEmooP_, hasMatchConv, missHits);
    bool isLoose = passIDWP("LOOSE",isEB, dEtaIn, dPhiIn, full5x5, hoe, d0, dz, ooEmooP_, hasMatchConv, missHits);
    bool isMedium = passIDWP("MEDIUM",isEB, dEtaIn, dPhiIn, full5x5, hoe, d0, dz, ooEmooP_, hasMatchConv, missHits);
    bool isTight = passIDWP("TIGHT",isEB, dEtaIn, dPhiIn, full5x5, hoe, d0, dz, ooEmooP_, hasMatchConv, missHits);
    // Look up the ID decision for this electron in 
    // the ValueMap object and store it. We need a Ptr object as the key.
    const Ptr<pat::Electron> elPtr(eleHandle, i);
    //std::cout<<"ValueMap: "<<veto_id_decisions->contains(elPtr.id())<<std::endl;
    bool vidVeto  = (*veto_id_decisions)[ elPtr ];
    bool vidLoose  = (*loose_id_decisions)[ elPtr ];
    bool vidMedium  = (*medium_id_decisions)[ elPtr ];
    bool vidTight = (*tight_id_decisions)[ elPtr ];
    bool vidHEEP  = (*HEEP_id_cutflow_data)[ elPtr ].cutFlowPassed();
    bool vidHEEP_noiso     = (*HEEP_id_cutflow_data)[ elPtr ].getCutFlowResultMasking(maskCuts).cutFlowPassed();
    //passVetoId_.push_back( isPassVeto );
    //passTightId_.push_back( isPassTight );
    if( verboseIdFlag_ ) {
      vid::CutFlowResult fullCutFlowData = (*medium_id_cutflow_data)[elPtr];
      edm::LogInfo("DEBUG:VID") << "CutFlow, full info for cand with pt= " << elPtr->pt();
      printCutFlowResult(fullCutFlowData);
    }

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
    el.addUserFloat("vidVeto",    vidVeto );
    el.addUserFloat("vidLoose",   vidLoose );
    el.addUserFloat("vidMedium",  vidMedium );
    el.addUserFloat("vidTight",   vidTight );
    el.addUserFloat("vidHEEP",    vidHEEP );
    el.addUserFloat("vidHEEPnoiso",    vidHEEP_noiso );



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
  if(abs(eta)>0.0 && abs(eta)<=1.0) effArea = 0.1752;
  if(abs(eta)>1.0 && abs(eta)<=1.479) effArea = 0.1862;
  if(abs(eta)>1.479 && abs(eta)<=2.0) effArea = 0.1411;
  if(abs(eta)>2.0 && abs(eta)<=2.2) effArea = 0.1534;
  if(abs(eta)>2.2 && abs(eta)<=2.3) effArea = 0.1903;
  if(abs(eta)>2.3 && abs(eta)<=2.4) effArea = 0.2243;
  if(abs(eta)>2.4 && abs(eta)<=2.5) effArea = 0.2687;
  return effArea;
}


float ElectronUserData::getEAWithGenWeightOld(float eta)
{
  // The following values refer to EA for cone 0.3 and fixedGridRhoFastjetAll. 
  // They are valid for electrons only, different EA are available for muons.
  float effArea = 0.;
  if(abs(eta)>0.0 && abs(eta)<=0.8) effArea = 0.0958;
  if(abs(eta)>0.8 && abs(eta)<=1.3) effArea = 0.0940;
  if(abs(eta)>1.3 && abs(eta)<=2.0) effArea = 0.0616;
  if(abs(eta)>2.0 && abs(eta)<=2.2) effArea = 0.0708;
  if(abs(eta)>2.2 && abs(eta)<=2.5) effArea = 0.1321;
  return effArea;
}


float ElectronUserData::getEAOld(float eta)
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

bool ElectronUserData::passIDWP(string WP, bool isEB, float dEtaIn, float dPhiIn, float full5x5, float hoe, float d0, float dz, float ooemoop, bool conv, int missHits){
  bool pass = false;

  if(WP == "VETO"){
    if(isEB){
      pass = (fabs(dEtaIn) <  0.0126 ) && (fabs(dPhiIn) <  0.107 ) && (full5x5 < 0.012 ) && (hoe <  0.186 ) && (fabs(d0) < 0.0621 ) && (fabs(dz) <  0.613 ) && (fabs(ooemoop) <  0.239 ) && !conv && (missHits <= 2);
    }
    else{
      pass = (fabs(dEtaIn) <  0.0109 ) && (fabs(dPhiIn) <  0.219 ) && (full5x5 < 0.0339 ) && (hoe <  0.0962 ) && (fabs(d0) < 0.279 ) && (fabs(dz) < 0.947 ) && (fabs(ooemoop) < 0.141 ) && !conv && (missHits <= 3);
    }
  }
  if(WP == "LOOSE"){
    if(isEB){
      pass = (fabs(dEtaIn) < 0.00976 ) && (fabs(dPhiIn) < 0.0929 ) && (full5x5 <  0.0105 ) && (hoe < 0.0765 ) && (fabs(d0) < 0.0227 ) && (fabs(dz) < 0.379 ) && (fabs(ooemoop) <  0.184 ) && !conv && (missHits <= 2);
    }
    else{
      pass = (fabs(dEtaIn) < 0.00952 ) && (fabs(dPhiIn) < 0.181 ) && (full5x5 < 0.0318 ) && (hoe < 0.0824 ) && (fabs(d0) < 0.242 ) && (fabs(dz) < 0.921 ) && (fabs(ooemoop) < 0.125 ) && !conv && (missHits <= 1);
    }
  }

  if(WP == "MEDIUM"){
    if(isEB){
      pass = (fabs(dEtaIn) <  0.0094 ) && (fabs(dPhiIn) <  0.0296 ) && (full5x5 <  0.0101 ) && (hoe <  0.0372 ) && (fabs(d0) <  0.0151 ) && (fabs(dz) <  0.238 ) && (fabs(ooemoop) <  0.118 ) && !conv && (missHits <= 2);
    }
    else{
      pass = (fabs(dEtaIn) <  0.00773 ) && (fabs(dPhiIn) <  0.148 ) && (full5x5 <  0.0287 ) && (hoe <  0.0546 ) && (fabs(d0) <  0.0535 ) && (fabs(dz) <  0.572 ) && (fabs(ooemoop) <  0.104 ) && !conv && (missHits <= 1);
    }
  }

  if(WP == "TIGHT"){
    if(isEB){
      pass = (fabs(dEtaIn) <  0.0095 ) && (fabs(dPhiIn) <  0.0291 ) && (full5x5 <  0.0101 ) && (hoe <  0.0372 ) && (fabs(d0) <  0.0144 ) && (fabs(dz) <  0.323 ) && (fabs(ooemoop) <  0.0174 ) && !conv && (missHits <= 2);
    }
    else{
      pass = (fabs(dEtaIn) <  0.00762 ) && (fabs(dPhiIn) <  0.0439 ) && (full5x5 <  0.0287 ) && (hoe <  0.0544 ) && (fabs(d0) <  0.0377 ) && (fabs(dz) <  0.571 ) && (fabs(ooemoop) <  0.01 ) && !conv && (missHits <= 1);
    }
  }
  return pass;
}

void ElectronUserData::printCutFlowResult(vid::CutFlowResult &cutflow){

  printf("    CutFlow name= %s    decision is %d\n", 
      cutflow.cutFlowName().c_str(),
      (int) cutflow.cutFlowPassed());
  int ncuts = cutflow.cutFlowSize();
  printf(" Index                               cut name              isMasked    value-cut-upon     pass?\n");
  for(int icut = 0; icut<ncuts; icut++){
    printf("  %2d      %50s    %d        %f          %d\n", icut,
        cutflow.getNameAtIndex(icut).c_str(),
        (int)cutflow.isCutMasked(icut),
        cutflow.getValueCutUpon(icut),
        (int)cutflow.getCutResultByIndex(icut));
  }
  printf("    WARNING: the value-cut-upon is bugged in 7.4.7, it is always 0.0 or 1.0\n");

}


#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(ElectronUserData);
