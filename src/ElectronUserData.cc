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
  EDGetTokenT< double > rho_miniIso_;
  EDGetTokenT< edm::TriggerResults > triggerResultsLabel_;
  EDGetTokenT< trigger::TriggerEvent > triggerSummaryLabel_;
  InputTag hltElectronFilterLabel_;
  TString hltPath_;
  HLTConfigProvider hltConfig;
  int triggerBit;
  int debug_; 
  // ----------member data ---------------------------
  edm::EDGetTokenT<edm::ValueMap<vid::CutFlowResult> > electronHEEPIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<vid::CutFlowResult> > eleVetoIdFullInfoMapToken_;
  edm::EDGetTokenT<edm::ValueMap<vid::CutFlowResult> > eleLooseIdFullInfoMapToken_;
  edm::EDGetTokenT<edm::ValueMap<vid::CutFlowResult> > eleMediumIdFullInfoMapToken_;
  edm::EDGetTokenT<edm::ValueMap<vid::CutFlowResult> > eleTightIdFullInfoMapToken_;
 
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
   rho_miniIso_(consumes<double>(edm::InputTag("fixedGridRhoFastjetCentralNeutral"))),
   triggerResultsLabel_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerResults"))),
   triggerSummaryLabel_(consumes<trigger::TriggerEvent>(iConfig.getParameter<edm::InputTag>("triggerSummary"))),
   hltElectronFilterLabel_ (iConfig.getParameter<edm::InputTag>("hltElectronFilter")),   //trigger objects we want to match
   hltPath_ (iConfig.getParameter<std::string>("hltPath")),
   electronHEEPIdMapToken_(consumes<edm::ValueMap<vid::CutFlowResult> >(iConfig.getParameter<edm::InputTag>("eleHEEPIdIdFullInfoMap"))),
   eleVetoIdFullInfoMapToken_(consumes<edm::ValueMap<vid::CutFlowResult> >(iConfig.getParameter<edm::InputTag>("eleVetoIdFullInfoMap"))),
   eleLooseIdFullInfoMapToken_(consumes<edm::ValueMap<vid::CutFlowResult> >(iConfig.getParameter<edm::InputTag>("eleLooseIdFullInfoMap"))),
   eleMediumIdFullInfoMapToken_(consumes<edm::ValueMap<vid::CutFlowResult> >(iConfig.getParameter<edm::InputTag>("eleMediumIdFullInfoMap"))),
   eleTightIdFullInfoMapToken_(consumes<edm::ValueMap<vid::CutFlowResult> >(iConfig.getParameter<edm::InputTag>("eleTightIdFullInfoMap"))),
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
  edm::Handle<double> rhoHandle, rhoHandle_miniIso;
  iEvent.getByToken(rho_,rhoHandle);
  iEvent.getByToken(rho_miniIso_,rhoHandle_miniIso);
  double rho = *rhoHandle, rho_miniIso = *rhoHandle_miniIso;

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
  edm::Handle<edm::ValueMap<vid::CutFlowResult> > veto_id_cutflow_data;
  edm::Handle<edm::ValueMap<vid::CutFlowResult> > loose_id_cutflow_data;
  edm::Handle<edm::ValueMap<vid::CutFlowResult> > medium_id_cutflow_data;
  edm::Handle<edm::ValueMap<vid::CutFlowResult> > tight_id_cutflow_data;
  edm::Handle<edm::ValueMap<vid::CutFlowResult> > heep_id_cutflow_data;
  iEvent.getByToken(eleVetoIdFullInfoMapToken_,veto_id_cutflow_data);
  iEvent.getByToken(eleLooseIdFullInfoMapToken_,loose_id_cutflow_data);
  iEvent.getByToken(eleMediumIdFullInfoMapToken_,medium_id_cutflow_data);
  iEvent.getByToken(eleTightIdFullInfoMapToken_,tight_id_cutflow_data);
  iEvent.getByToken(electronHEEPIdMapToken_,heep_id_cutflow_data);
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


    //Cut variables     
    float dEtaInSeed = el.superCluster().isNonnull() && el.superCluster()->seed().isNonnull() ? el.deltaEtaSuperClusterTrackAtVtx() - el.superCluster()->eta() + el.superCluster()->seed()->eta() : std::numeric_limits<float>::max();  // the cutID uses the absolute of this variable
    float dPhiIn = el.deltaPhiSuperClusterTrackAtVtx(); // the cutID uses the absolute of this variable
    float full5x5 = el.full5x5_sigmaIetaIeta();
    float hoe = el.hadronicOverEm();
   
    double EA = getEA(el.eta());
    float absiso_EA = pfIso.sumChargedHadronPt + std::max(0.0 , pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - rho * EA );
    float relIsoWithEA_ = absiso_EA/el.pt();
    float missHits = el.gsfTrack()->hitPattern().numberOfHits(HitPattern::MISSING_INNER_HITS);
   

    float ooEmooP_; 
    if( el.ecalEnergy() == 0 ){
      printf("Electron energy is zero!\n");
      ooEmooP_ = 999;
    }else if( !std::isfinite(el.ecalEnergy())){
      printf("Electron energy is not finite!\n");
      ooEmooP_ = 998;
    }else{
      ooEmooP_ = std::abs(1.0/el.ecalEnergy() - el.eSuperClusterOverP()/el.ecalEnergy() );
    }



    //Other variables
    float absisoWithDBeta = pfIso.sumChargedHadronPt + max(0.0 , pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - 0.5 * pfIso.sumPUPt );
    float relIsoWithDBeta_ = absisoWithDBeta/el.pt();
    double miniIso = getPFMiniIsolation(packedPFCands, dynamic_cast<const reco::Candidate *>(&el), 0.05, 0.2, 10., false, true, EA, rho_miniIso);
   
    // Impact parameter
    float dxy   = el.gsfTrack()->dxy(vtxPoint);
    float dz    = el.gsfTrack()->dz(vtxPoint);
    float dB    = el.dB (pat::Electron::PV3D);
    float dBErr = el.edB(pat::Electron::PV3D);


    if(debug_>=1) cout<<" ele " << i <<" pt "<< el.pt()<<" eta "<<el.eta()<<"fabs(1/E-1/P) "<< ooEmooP_ <<" dxy "<< dxy <<" dz " << dz <<" iso " << relIsoWithDBeta_<<endl; 		    
     // conversion rejection match
    bool hasMatchConv = ConversionTools::hasMatchedConversion(el, conversions, beamspot.position());


    // Look up the ID decision for this electron in 
    // the ValueMap object and store it. We need a Ptr object as the key.
    const Ptr<pat::Electron> elPtr(eleHandle, i);
    bool vidVeto  = (*veto_id_cutflow_data)[ elPtr ].cutFlowPassed();
    bool vidLoose  = (*loose_id_cutflow_data)[ elPtr ].cutFlowPassed();
    bool vidMedium  = (*medium_id_cutflow_data)[ elPtr ].cutFlowPassed();
    bool vidTight = (*tight_id_cutflow_data)[ elPtr ].cutFlowPassed();
    bool vidHEEP  = (*heep_id_cutflow_data)[ elPtr ].cutFlowPassed();


    //retrieving bits from fullflowcutData  and masking rel iso EA cut isolation cut
    const int cutIndexToMask = 7;     // this is the relative iso cut index - one can verify with printing out the full info 
    vid::CutFlowResult fullCutFlowData = (*veto_id_cutflow_data)[elPtr]; 
    vid::CutFlowResult maskedCutFlowData = fullCutFlowData.getCutFlowResultMasking(cutIndexToMask);
    bool vidVeto_noiso    =  (int) maskedCutFlowData.cutFlowPassed();

    fullCutFlowData = (*loose_id_cutflow_data)[elPtr];
    maskedCutFlowData = fullCutFlowData.getCutFlowResultMasking(cutIndexToMask);
    bool vidLoose_noiso   =  (int) maskedCutFlowData.cutFlowPassed();

    fullCutFlowData = (*medium_id_cutflow_data)[elPtr];
    maskedCutFlowData = fullCutFlowData.getCutFlowResultMasking(cutIndexToMask);
    bool vidMedium_noiso  =  (int) maskedCutFlowData.cutFlowPassed();

    fullCutFlowData = (*tight_id_cutflow_data)[elPtr];
    maskedCutFlowData = fullCutFlowData.getCutFlowResultMasking(cutIndexToMask);
    bool vidTight_noiso   =  (int) maskedCutFlowData.cutFlowPassed();

    fullCutFlowData = (*heep_id_cutflow_data)[elPtr];
    maskedCutFlowData = fullCutFlowData.getCutFlowResultMasking(7); //7 is for trackIsoPt and 8 for hademd1 check printout
    vid::CutFlowResult maskedCutFlowData2 = fullCutFlowData.getCutFlowResultMasking(cutIndexToMask);
    maskedCutFlowData2 = fullCutFlowData.getCutFlowResultMasking(8);
    bool vidHEEP_noiso     = (int) maskedCutFlowData.cutFlowPassed() && (int) maskedCutFlowData2.cutFlowPassed() ;    

    


    if( verboseIdFlag_ ) {
      vid::CutFlowResult fullCutFlowData = (*medium_id_cutflow_data)[elPtr];
      edm::LogInfo("DEBUG:VID") << "CutFlow, full info for cand with pt= " << elPtr->pt();
      printCutFlowResult(fullCutFlowData);
    }

    el.addUserFloat("dEtaInSeed",  dEtaInSeed);
    el.addUserFloat("dPhiIn",      dPhiIn);
    el.addUserFloat("full5x5",     full5x5);
    el.addUserFloat("hoe",         hoe);
    el.addUserFloat("dxy",         dxy);
    el.addUserFloat("dz",          dz);
    el.addUserFloat("dB",          dB);
    el.addUserFloat("dBErr",       dBErr);
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
    el.addUserFloat("vidVetonoiso",    vidVeto_noiso );
    el.addUserFloat("vidLoosenoiso",    vidLoose_noiso );
    el.addUserFloat("vidMediumnoiso",    vidMedium_noiso );
    el.addUserFloat("vidTightnoiso",    vidTight_noiso );
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
  //These are Effective areas suitable for 80X samples post ICHEP
  float effArea = 0.;
  if(abs(eta)>0.0 && abs(eta)<=1.0)   effArea = 0.1703;
  if(abs(eta)>1.0 && abs(eta)<=1.479) effArea = 0.1715;
  if(abs(eta)>1.479 && abs(eta)<=2.0) effArea = 0.1213;
  if(abs(eta)>2.0 && abs(eta)<=2.2)   effArea = 0.1230;
  if(abs(eta)>2.2 && abs(eta)<=2.3)   effArea = 0.1635;
  if(abs(eta)>2.3 && abs(eta)<=2.4)   effArea = 0.1937;
  if(abs(eta)>2.4 && abs(eta)<=5.0)   effArea = 0.2393;
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
