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
  bool isMatchedWithTrigger(const pat::Electron, trigger::TriggerObjectCollection,int&);
  void printCutFlowResult(vid::CutFlowResult &cutflow);
  float getEA(float eta);
  EDGetTokenT< std::vector< pat::Electron > > eleLabel_;
  EDGetTokenT< std::vector< reco::Vertex > > pvLabel_;
  EDGetTokenT< pat::PackedCandidateCollection > packedPFCandsLabel_;
  EDGetTokenT< reco::ConversionCollection > convLabel_;
  EDGetTokenT< reco::BeamSpot > beamLabel_;
  EDGetTokenT< double > rho_;
  EDGetTokenT< edm::TriggerResults > triggerResultsLabel_;
  EDGetTokenT< trigger::TriggerEvent > triggerSummaryLabel_;
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
   eleLabel_(consumes<std::vector<pat::Electron>>(iConfig.getParameter<edm::InputTag>("eleLabel"))), 
   pvLabel_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("pv"))), // "offlinePrimaryVertex"
   packedPFCandsLabel_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("packedPFCands"))), 
   convLabel_(consumes<reco::ConversionCollection>(iConfig.getParameter<edm::InputTag>("conversion"))),   // "offlinePrimaryVertex"
   beamLabel_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"))), 
   rho_(consumes<double>(iConfig.getParameter<edm::InputTag>("rho"))), 
   triggerResultsLabel_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerResults"))),
   triggerSummaryLabel_(consumes<trigger::TriggerEvent>(iConfig.getParameter<edm::InputTag>("triggerSummary"))),
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
  iEvent.getByToken(pvLabel_, vertices);

  //RHO
  edm::Handle<double> rhoHandle;
  iEvent.getByToken(rho_,rhoHandle);
  double rho = *rhoHandle; 

  if(debug_>=1) cout<<"vtx size " << vertices->size()<<endl; 

  reco::TrackBase::Point vtxPoint(0,0, 0);
  if(  vertices->size() >= 1 ) {
    vtxPoint = vertices->at(0).position();
  }

  if(debug_>=1) cout<<"vtxPoint " <<vtxPoint.X()<<" "<< vtxPoint.Y()<<" "<< vtxPoint.Z()<<endl; 
    
  //Electrons
  edm::Handle<std::vector<pat::Electron> > eleHandle;
  iEvent.getByToken(eleLabel_, eleHandle);
  
  auto_ptr<vector<pat::Electron> > eleColl( new vector<pat::Electron> (*eleHandle) );

  //PackedPFCands for Mini-iso
  edm::Handle<pat::PackedCandidateCollection> packedPFCands;
  iEvent.getByToken(packedPFCandsLabel_, packedPFCands);

  Handle<reco::BeamSpot> bsHandle;
  iEvent.getByToken(beamLabel_, bsHandle);
  const reco::BeamSpot &beamspot = *bsHandle.product();
    
  Handle<reco::ConversionCollection> conversions;
  iEvent.getByToken(convLabel_, conversions);
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
  iEvent.getByToken(triggerResultsLabel_, triggerResults);
  if (size_t(triggerBit) < triggerResults->size() )
    if (triggerResults->accept(triggerBit))
      std::cout << "event pass : " << hltPath_ << std::endl;
  // if (triggerResults->accept(triggerBit)){ //crashing here to be understood!
  //   std::cout << "event pass : " << hltPath_ << std::endl;
  // }

  edm::Handle<trigger::TriggerEvent> triggerSummary;
  iEvent.getByToken(triggerSummaryLabel_, triggerSummary);

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


float ElectronUserData::getEA(float eta){
  float effArea = 0.;
  if(abs(eta)>0.0 && abs(eta)<=1.0)   effArea = 0.1752;
  if(abs(eta)>1.0 && abs(eta)<=1.479) effArea = 0.1862;
  if(abs(eta)>1.479 && abs(eta)<=2.0) effArea = 0.1411;
  if(abs(eta)>2.0 && abs(eta)<=2.2)   effArea = 0.1534;
  if(abs(eta)>2.2 && abs(eta)<=2.3)   effArea = 0.1903;
  if(abs(eta)>2.3 && abs(eta)<=2.4)   effArea = 0.2243;
  if(abs(eta)>2.4 && abs(eta)<=2.5)   effArea = 0.2687;
  return effArea;
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
