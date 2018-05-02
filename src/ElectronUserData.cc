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

  //VID now embedded in slimmed electrons, can just be called from ElectronID.
  bool verboseIdFlag_;
  TString vidLooseLabel_,vidMediumLabel_,vidTightLabel_,vidVetoLabel_;
  TString vidLooseNoIsoLabel_,vidMediumNoIsoLabel_,vidTightNoIsoLabel_,vidVetoNoIsoLabel_;
  TString vidLooseMVALabel_,vidMediumMVALabel_,vidTightMVALabel_,vidVetoMVALabel_;
  TString vidLooseMVANoIsoLabel_,vidMediumMVANoIsoLabel_,vidTightMVANoIsoLabel_,vidVetoMVANoIsoLabel_;
  TString vidHEEPLabel_;


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
   verboseIdFlag_(iConfig.getParameter<bool>("eleIdVerbose")),
   vidLooseLabel_(iConfig.getUntrackedParameter<string>("vidLooseLabel","cutBasedElectronID-Fall17-94X-V1-loose")),
   vidMediumLabel_(iConfig.getUntrackedParameter<string>("vidMediumLabel","cutBasedElectronID-Fall17-94X-V1-medium")),
   vidTightLabel_(iConfig.getUntrackedParameter<string>("vidTightLabel","cutBasedElectronID-Fall17-94X-V1-tight")),
   vidVetoLabel_(iConfig.getUntrackedParameter<string>("vidVetoLabel","cutBasedElectronID-Fall17-94X-V1-veto")),

   //No iso available only for MVA at the moment
   vidLooseNoIsoLabel_(iConfig.getUntrackedParameter<string>("vidLooseNoIsoLabel","")),
   vidMediumNoIsoLabel_(iConfig.getUntrackedParameter<string>("vidMediumNoIsoLabel","")),
   vidTightNoIsoLabel_(iConfig.getUntrackedParameter<string>("vidTightNoIsoLabel","")),
   vidVetoNoIsoLabel_(iConfig.getUntrackedParameter<string>("vidVetoNoIsoLabel","")),

   vidLooseMVALabel_(iConfig.getUntrackedParameter<string>("vidLooseMVALabel","mvaEleID-Fall17-iso-V1-wpLoose")),
   vidMediumMVALabel_(iConfig.getUntrackedParameter<string>("vidMediumMVALabel","mvaEleID-Fall17-iso-V1-wp90")),
   vidTightMVALabel_(iConfig.getUntrackedParameter<string>("vidTightMVALabel","mvaEleID-Fall17-iso-V1-wp80")),
   vidVetoMVALabel_(iConfig.getUntrackedParameter<string>("vidVetoMVALabel","")),

   vidLooseMVANoIsoLabel_(iConfig.getUntrackedParameter<string>("vidLooseNoIsoLabel","mvaEleID-Fall17-noIso-V1-wpLoose")),
   vidMediumMVANoIsoLabel_(iConfig.getUntrackedParameter<string>("vidMediumNoIsoLabel","mvaEleID-Fall17-noIso-V1-wp90")),
   vidTightMVANoIsoLabel_(iConfig.getUntrackedParameter<string>("vidTightNoIsoLabel","mvaEleID-Fall17-noIso-V1-wp80")),
   vidVetoMVANoIsoLabel_(iConfig.getUntrackedParameter<string>("vidVetoNoIsoLabel","")),

   vidHEEPLabel_(iConfig.getUntrackedParameter<string>("vidHEEPLabel","heepElectronID-HEEPV70"))  
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
  
  std::unique_ptr<vector<pat::Electron> > eleColl( new vector<pat::Electron> (*eleHandle) );

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
  //passVetoId_.clear();     
  //passTightId_.clear();

  // Cuts to be masked for non-iso version. 
  vector<string> maskCuts;
  maskCuts.push_back("GsfEleTrkPtIsoCut_0"); 
  maskCuts.push_back("GsfEleEmHadD1IsoRhoCut_0");   
  
  //// TRIGGER
  //bool changedConfig = false;
  //if (!hltConfig.init(iEvent.getRun(), iSetup, "HLT", changedConfig)) {
  //  cout << "Initialization of HLTConfigProvider failed!!" << endl;
  //  return;
  //}

  //
  //if (changedConfig){
  //  std::cout << "the curent menu is " << hltConfig.tableName() << " "<< hltConfig.triggerNames().size() <<endl; 
  //  triggerBit = -1;
  //  for (size_t j = 0; j < hltConfig.triggerNames().size(); j++) {
  //    if(debug_ ==99) cout<<"trigName " <<j<<" "<< hltConfig.triggerNames()[j] <<endl; 
  //    if (TString(hltConfig.triggerNames()[j]).Contains(hltPath_)) triggerBit = j;
  //  }
  //  if (triggerBit == -1) cout << "HLT path not found" << endl;
  //}
  //edm::Handle<edm::TriggerResults> triggerResults;
  //iEvent.getByToken(triggerResultsLabel_, triggerResults);
  //if (size_t(triggerBit) < triggerResults->size() )
  //  if (triggerResults->accept(triggerBit))
  //    std::cout << "event pass : " << hltPath_ << std::endl;
  //// if (triggerResults->accept(triggerBit)){ //crashing here to be understood!
  ////   std::cout << "event pass : " << hltPath_ << std::endl;
  //// }

  //edm::Handle<trigger::TriggerEvent> triggerSummary;
  //iEvent.getByToken(triggerSummaryLabel_, triggerSummary);

  ////find the index corresponding to the event
  //if(false){ ///to be done after understanding which trigger(s) to use.
  //size_t ElectronFilterIndex = (*triggerSummary).filterIndex(hltElectronFilterLabel_);
  //trigger::TriggerObjectCollection allTriggerObjects = triggerSummary->getObjects();
  //trigger::TriggerObjectCollection ElectronLegObjects;
  //if (ElectronFilterIndex < (*triggerSummary).sizeFilters()) { //check if the trigger object is present
  //  //save the trigger objects corresponding to ele leg
  //  const trigger::Keys &keysElectrons = (*triggerSummary).filterKeys(ElectronFilterIndex);
  //  for (size_t j = 0; j < keysElectrons.size(); j++) {
  //    trigger::TriggerObject foundObject = (allTriggerObjects)[keysElectrons[j]];
  //    ElectronLegObjects.push_back(foundObject);
  //  }
  //}
  //if( debug_ >=1)  std::cout << "ElectronLegObjects: " << ElectronLegObjects.size() << std::endl;
  ///////////  /////////  /////////  /////////  /////////  /////////  /////////  /////////  /////////  /////////
  //}

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
    float missHits = el.gsfTrack()->hitPattern().numberOfAllHits(HitPattern::MISSING_INNER_HITS);

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
    bool vidVeto   = true; 
    bool vidLoose   = true; 
    bool vidMedium   = true;
    bool vidTight  = true; 

    bool vidVeto_noiso   = true; 
    bool vidLoose_noiso   = true; 
    bool vidMedium_noiso   = true;
    bool vidTight_noiso  = true; 

    bool vidVeto_mva   = true; 
    bool vidLoose_mva   = true; 
    bool vidMedium_mva   = true;
    bool vidTight_mva  = true; 

    bool vidVeto_mva_noiso   = true; 
    bool vidLoose_mva_noiso   = true; 
    bool vidMedium_mva_noiso   = true;
    bool vidTight_mva_noiso  = true; 
    
    bool vidHEEP   = true; 
    bool vidHEEP_noiso   = true; 
    

    // if(gp_mva_val > 0.1 ) cout<<"true ele "<<mvaval<<endl;

    //    cout << " electron "<<
    (vidLooseLabel_=="") ? vidLoose=0 : vidLoose=el.electronID(vidLooseLabel_);
    if(debug_>=1) cout << " vid loose label "<< vidLooseLabel_ << " vidLoose "<< vidLoose <<endl;; 
    (vidMediumLabel_=="") ? vidMedium=0 : vidMedium=el.electronID(vidMediumLabel_);
    if(debug_>=1) cout << " vid medium label "<< vidMediumLabel_ << " vidMedium "<< vidMedium <<endl;; 

    (vidTightLabel_=="") ? vidTight=0 : vidTight=el.electronID(vidTightLabel_);
    (vidVetoLabel_=="") ? vidVeto=0 : vidVeto=el.electronID(vidVetoLabel_);
    
    (vidLooseNoIsoLabel_=="") ? vidLoose_noiso=0 : vidLoose_noiso=el.electronID(vidLooseLabel_);
    (vidMediumNoIsoLabel_=="") ? vidMedium_noiso=0 : vidMedium_noiso=el.electronID(vidMediumLabel_);
    (vidTightNoIsoLabel_=="") ? vidTight_noiso=0 : vidTight_noiso=el.electronID(vidTightLabel_);
    (vidVetoNoIsoLabel_=="") ? vidVeto_noiso=0 : vidVeto_noiso=el.electronID(vidVetoLabel_);
    
    (vidLooseMVALabel_=="") ? vidLoose_mva=0 : vidLoose_mva=el.electronID(vidLooseMVALabel_);
    (vidMediumMVALabel_=="") ? vidMedium_mva=0 : vidMedium_mva=el.electronID(vidMediumMVALabel_);
    (vidTightMVALabel_=="") ? vidTight_mva=0 : vidTight_mva=el.electronID(vidTightMVALabel_);
    (vidVetoMVALabel_=="") ? vidVeto_mva=0 : vidVeto_mva=el.electronID(vidVetoMVALabel_);

    (vidLooseMVANoIsoLabel_=="") ? vidLoose_mva_noiso=0 : vidLoose_mva_noiso=el.electronID(vidLooseMVANoIsoLabel_);
    (vidMediumMVANoIsoLabel_=="") ? vidMedium_mva_noiso=0 : vidMedium_mva_noiso=el.electronID(vidMediumMVANoIsoLabel_);
    (vidTightMVANoIsoLabel_=="") ? vidTight_mva_noiso=0 : vidTight_mva_noiso=el.electronID(vidTightMVANoIsoLabel_);
    (vidVetoMVANoIsoLabel_=="") ? vidVeto_mva_noiso=0 : vidVeto_mva_noiso=el.electronID(vidVetoMVANoIsoLabel_);

    (vidHEEPLabel_=="") ? vidHEEP=0 : vidHEEP=el.electronID(vidHEEPLabel_);


    //    if( verboseIdFlag_ ) {
    //  vid::CutFlowResult fullCutFlowData = (*medium_id_cutflow_data)[elPtr];
    //  edm::LogInfo("DEBUG:VID") << "CutFlow, full info for cand with pt= " << elPtr->pt();
    //  printCutFlowResult(fullCutFlowData);
    // }

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

    el.addUserFloat("vidMVAVeto",    vidVeto_mva );
    el.addUserFloat("vidMVALoose",   vidLoose_mva );
    el.addUserFloat("vidMVAMedium",  vidMedium_mva );
    el.addUserFloat("vidMVATight",   vidTight_mva );

    el.addUserFloat("vidMVAVetonoiso",    vidVeto_mva_noiso );
    el.addUserFloat("vidMVALoosenoiso",   vidLoose_mva_noiso );
    el.addUserFloat("vidMVAMediumnoiso",  vidMedium_mva_noiso );
    el.addUserFloat("vidMVATightnoiso",   vidTight_mva_noiso );
   
  }

  iEvent.put( std::move(eleColl) );

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
