#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDMException.h"
// dR and dPhi                                                                                                                                                              
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"


// Photon dependencies
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "RecoEgamma/EgammaTools/interface/EffectiveAreas.h" 

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/Candidate/interface/Candidate.h"
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

class  PhotonUserData : public edm::EDProducer {
public:
  PhotonUserData( const edm::ParameterSet & );   

private:
  void produce( edm::Event &, const edm::EventSetup & );
  bool isMatchedWithTrigger();
  bool passIDWP();
  float EA25nsPho(float eta);
  float EA25nsNeu(float eta);
  EDGetTokenT< double > rhoLabel_;
  EDGetTokenT< std::vector< pat::Photon > > phoLabel_;
  edm::EDGetToken                      photonsMiniAODToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > phoLooseIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > phoMediumIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > phoTightIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<float> > phoISOCMapToken_;
  edm::EDGetTokenT<edm::ValueMap<float> > phoISOPMapToken_;
  edm::EDGetTokenT<edm::ValueMap<float> > phoISONMapToken_;
  edm::EDGetTokenT<edm::ValueMap<float> > full5x5SigmaIEtaIEtaMapToken_;
  
  InputTag triggerResultsLabel_, triggerSummaryLabel_;
  InputTag hltElectronFilterLabel_;
  TString hltPath_;
  HLTConfigProvider hltConfig;
  int triggerBit;
  int debug_; 
  
};

//Necessary for the  Isolation Rho Correction

PhotonUserData::PhotonUserData(const edm::ParameterSet& iConfig):
  rhoLabel_(consumes<double>(iConfig.getParameter<edm::InputTag>("rho"))), //rhofixedgridRhoFastjet All"
  phoLabel_(consumes<std::vector<pat::Photon>>(iConfig.getParameter<edm::InputTag>("pholabel"))), 
  photonsMiniAODToken_(consumes<edm::View<pat::Photon> >(iConfig.getParameter<edm::InputTag>("pholabel"))),
  phoLooseIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("phoLooseIdMap"))),
  phoMediumIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("phoMediumIdMap"))),
  phoTightIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("phoTightIdMap"))),
  phoISOCMapToken_(consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("phoChgIsoMap"))),
  phoISOPMapToken_(consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("phoPhoIsoMap"))),
  phoISONMapToken_(consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("phoNeuIsoMap"))),
  full5x5SigmaIEtaIEtaMapToken_(consumes <edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("full5x5SigmaIEtaIEtaMap")))
{
  debug_ = iConfig.getUntrackedParameter<int>("debugLevel",int(0));
  produces<std::vector<pat::Photon> >();
}

void PhotonUserData::produce( edm::Event& iEvent, const edm::EventSetup& iSetup) {

  edm::Handle<edm::View<pat::Photon> > photonS;
  edm::Handle< double > rhoH;
  edm::Handle<edm::ValueMap<bool> > loose_id_decisions;
  edm::Handle<edm::ValueMap<bool> > medium_id_decisions;
  edm::Handle<edm::ValueMap<bool> > tight_id_decisions;
  edm::Handle<edm::ValueMap<float> > isoc_idvar;
  edm::Handle<edm::ValueMap<float> > isop_idvar;
  edm::Handle<edm::ValueMap<float> > ison_idvar;
  edm::Handle<edm::ValueMap<float> > full5x5SigmaIEtaIEtaMap;
  iEvent.getByToken(photonsMiniAODToken_,photonS);
  iEvent.getByToken(rhoLabel_,rhoH);
  iEvent.getByToken(phoLooseIdMapToken_ ,loose_id_decisions);
  iEvent.getByToken(phoMediumIdMapToken_,medium_id_decisions);
  iEvent.getByToken(phoTightIdMapToken_ ,tight_id_decisions);
  iEvent.getByToken(phoISOCMapToken_ ,isoc_idvar);
  iEvent.getByToken(phoISOPMapToken_ ,isop_idvar);
  iEvent.getByToken(phoISONMapToken_ ,ison_idvar);
  iEvent.getByToken(full5x5SigmaIEtaIEtaMapToken_, full5x5SigmaIEtaIEtaMap);


  edm::Handle<std::vector<pat::Photon> > photonS2;
  iEvent.getByToken(phoLabel_, photonS2);
  auto_ptr<vector<pat::Photon> > phoColl( new vector<pat::Photon> (*photonS2) );

  //rho
  float rho_;
  rho_ = *rhoH;

  //Photon Loop

  for (size_t i = 0; i < phoColl->size(); ++i){
    const auto pho = photonS->ptrAt(i);
    pat::Photon & phoi = (*phoColl)[i];
   if(phoi.hadTowOverEm() > 0.15 ) continue;
   
    float pho_isoC       = 1*rho_; //(*isoc_idvar)[pho];
  
    
    bool pho_isLoose  = (*loose_id_decisions)[pho];
    bool pho_isMedium = (*medium_id_decisions)[pho];
    bool pho_isTight  = (*tight_id_decisions)[pho];
  
    //Isolations raw and EA corrected
    float pho_isoP       = (*isop_idvar)[pho];
    float pho_isoN       = (*ison_idvar)[pho];
    
    
    float abseta = fabs(pho->eta());    

    float pho_isoCea     = std::max( float(0.0) ,(*isoc_idvar)[pho]                         );
    float pho_isoPea     = std::max( float(0.0) ,(*isop_idvar)[pho] - rho_*EA25nsPho(abseta)) ;
    float pho_isoNea     = std::max( float(0.0) ,(*ison_idvar)[pho] - rho_*EA25nsNeu(abseta));


    //showershapes 
    float pho_r9 = pho->r9();
    float pho_sieie = (*full5x5SigmaIEtaIEtaMap)[pho];
    float pho_hoe = pho->hadTowOverEm();
    

    bool hasPixelSeed = pho->hasPixelSeed(); 

    //Shower shapes
    phoi.addUserInt("hasPixelSeed", hasPixelSeed);
    phoi.addUserFloat("sigmaIetaIeta", pho_sieie);
    phoi.addUserFloat("hoe", pho_hoe);
    phoi.addUserFloat("r9",  pho_r9);

    
    //isolation
    phoi.addUserFloat("isoC",pho_isoC);
  
    
    phoi.addUserFloat("isoP",pho_isoP);
    phoi.addUserFloat("isoN",pho_isoN);
    phoi.addUserFloat("isoC_EAcor",pho_isoCea);
    phoi.addUserFloat("isoP_EAcor",pho_isoPea);
    phoi.addUserFloat("isoN_EAcor",pho_isoNea);

    phoi.addUserInt("isLoose",    pho_isLoose);
    phoi.addUserInt("isMedium",   pho_isMedium);
    phoi.addUserInt("isTight",    pho_isTight);
    




  }
//EOF photons loop
  iEvent.put( phoColl );

}

float
PhotonUserData::EA25nsPho(float eta){
  //Returns the Effective areas for the PF:gamma Isolation

 float effArea = 0; 
  if(abs(eta)>0.0 && abs(eta)<=1.0) effArea = 0.1271;
  if(abs(eta)>1.0 && abs(eta)<=1.479) effArea = 0.1101;
  if(abs(eta)>1.479 && abs(eta)<=2.0) effArea = 0.0756;
  if(abs(eta)>2.0 && abs(eta)<=2.2) effArea = 0.1175;
  if(abs(eta)>2.2 && abs(eta)<=2.3) effArea = 0.1498;
  if(abs(eta)>2.3 && abs(eta)<=2.4) effArea = 0.1857;
  if(abs(eta)>2.4 && abs(eta)<=2.5) effArea = 0.2183;
  return effArea;

}


float
PhotonUserData::EA25nsNeu(float eta){
  //Returns the Effective areas for the PF:h0 Isolation
  float effArea = 0; 
  if(abs(eta)>0.0 && abs(eta)<=1.0) effArea = 0.0599;
  if(abs(eta)>1.0 && abs(eta)<=1.479) effArea = 0.0819;
  if(abs(eta)>1.479 && abs(eta)<=2.0) effArea = 0.0696;
  if(abs(eta)>2.0 && abs(eta)<=2.2) effArea = 0.0360;
  if(abs(eta)>2.2 && abs(eta)<=2.3) effArea = 0.0360;
  if(abs(eta)>2.3 && abs(eta)<=2.4) effArea = 0.0462;
  if(abs(eta)>2.4 && abs(eta)<=2.5) effArea = 0.0656;
  return effArea;

}





bool
PhotonUserData::isMatchedWithTrigger(){
  return true;
}


bool PhotonUserData::passIDWP(){
 
  return true;
}



#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PhotonUserData);
