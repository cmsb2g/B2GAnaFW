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
  float EASP16(float eta,int type);

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

  edm::Handle<edm::View<pat::Photon> > photonHandle;
  edm::Handle< double > rhoH;
  edm::Handle<edm::ValueMap<bool> > loose_id_decisions;
  edm::Handle<edm::ValueMap<bool> > medium_id_decisions;
  edm::Handle<edm::ValueMap<bool> > tight_id_decisions;
  edm::Handle<edm::ValueMap<float> > isoc_idvar;
  edm::Handle<edm::ValueMap<float> > isop_idvar;
  edm::Handle<edm::ValueMap<float> > ison_idvar;
  edm::Handle<edm::ValueMap<float> > full5x5SigmaIEtaIEtaMap;
  iEvent.getByToken(photonsMiniAODToken_,photonHandle);
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
  std::unique_ptr<vector<pat::Photon> > phoColl( new vector<pat::Photon> (*photonS2) );

  //rho
  float rho_;
  rho_ = *rhoH;

  //Photon Loop

  for (size_t i = 0; i < phoColl->size(); ++i){
    //const Ptr<pat::Photon> pho(photonHandle,i); 
    const auto pho = photonHandle->ptrAt(i);
    pat::Photon & phoi = (*phoColl)[i];
  
    bool pho_isLoose  = (*loose_id_decisions)[pho];
    bool pho_isMedium = (*medium_id_decisions)[pho];
    bool pho_isTight  = (*tight_id_decisions)[pho];
  
    //Isolations raw and EA corrected
    float pho_isoP       = (*isop_idvar)[pho];
    float pho_isoN       = (*ison_idvar)[pho];
    float pho_isoC       = (*isoc_idvar)[pho];
    float abseta         = fabs(pho->eta());    
    float pho_isoPea     = std::max( float(0.0) ,(*isop_idvar)[pho] - rho_*EASP16(abseta,2)) ;
    float pho_isoNea     = std::max( float(0.0) ,(*ison_idvar)[pho] - rho_*EASP16(abseta,1));
    float pho_isoCea     = std::max( float(0.0) ,(*ison_idvar)[pho] - rho_*EASP16(abseta,0));


    //showershapes 
    float pho_r9 = pho->r9();
    float pho_sieie = (*full5x5SigmaIEtaIEtaMap)[pho];
    float pho_hoe = pho->hadTowOverEm();
    float pho_sieip = pho->sep(); 
    float pho_sipip = pho->spp(); 
  
    float pho_e1x5     = pho->e1x5(); 
    float pho_e5x5     = pho->e5x5(); 
  

    //Photon Ele discrimination 

    bool hasPixelSeed = pho->hasPixelSeed(); 
    bool pho_eleveto  = pho->passElectronVeto();

    phoi.addUserInt("hasPixelSeed", hasPixelSeed);
    phoi.addUserInt("eleveto", pho_eleveto); //
    phoi.addUserFloat("sigmaIetaIeta", pho_sieie);
    phoi.addUserFloat("sigmaIetaIphi",pho_sieip); //
    phoi.addUserFloat("sigmaIphiIphi",pho_sipip); //
    phoi.addUserFloat("e1x5",pho_e1x5); //
    phoi.addUserFloat("e5x5",pho_e5x5); //
    phoi.addUserFloat("hoe", pho_hoe);
    phoi.addUserFloat("r9",  pho_r9);
    
    
    //isolation
    phoi.addUserFloat("isoC",pho_isoC);
    phoi.addUserFloat("isoP",pho_isoP);
    phoi.addUserFloat("isoN",pho_isoN);
    phoi.addUserFloat("isoP_EAcor",pho_isoPea);
    phoi.addUserFloat("isoN_EAcor",pho_isoNea);
    phoi.addUserFloat("isoC_EAcor",pho_isoCea);

    phoi.addUserInt("isLoose",    pho_isLoose);
    phoi.addUserInt("isMedium",   pho_isMedium);
    phoi.addUserInt("isTight",    pho_isTight);

  }
//EOF photons loop
  iEvent.put( std::move(phoColl) );

}

float
PhotonUserData::EASP16(float eta,int type){
  //Returns the Effective areas for the PF:gamma Isolation
  //Taken from here: https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2#Selection_implementation_details
  
 float effArea = 0; 

 //----type 0 charged hadron  EAs
  if(abs(eta)>0.0 && abs(eta)<=1.0 && type == 0) effArea = 0.0360;
  if(abs(eta)>1.0 && abs(eta)<=1.479 && type == 0) effArea = 0.0377;
  if(abs(eta)>1.479 && abs(eta)<=2.0 && type == 0) effArea = 0.0306;
  if(abs(eta)>2.0 && abs(eta)<=2.2 && type == 0) effArea = 0.0283;
  if(abs(eta)>2.2 && abs(eta)<=2.3 && type == 0) effArea = 0.0254;
  if(abs(eta)>2.3 && abs(eta)<=2.4 && type == 0) effArea = 0.0217;
  if(abs(eta)>2.4 && abs(eta)<=2.5 && type == 0) effArea = 0.0167;

   //----type 1 neutral hadron  EAs
  if(abs(eta)>0.0 && abs(eta)<=1.0 && type == 1) effArea = 0.0597;
  if(abs(eta)>1.0 && abs(eta)<=1.479 && type == 1) effArea = 0.0807;
  if(abs(eta)>1.479 && abs(eta)<=2.0 && type == 1) effArea = 0.0629;
  if(abs(eta)>2.0 && abs(eta)<=2.2 && type == 1) effArea = 0.0197;
  if(abs(eta)>2.2 && abs(eta)<=2.3 && type == 1) effArea = 0.0184;
  if(abs(eta)>2.3 && abs(eta)<=2.4 && type == 1) effArea = 0.0284;
  if(abs(eta)>2.4 && abs(eta)<=2.5 && type == 1) effArea = 0.0591;

   //----type 2 photons  EAs
  if(abs(eta)>0.0 && abs(eta)<=1.0 && type == 2) effArea = 0.1210;
  if(abs(eta)>1.0 && abs(eta)<=1.479 && type == 2) effArea = 0.1107;
  if(abs(eta)>1.479 && abs(eta)<=2.0 && type == 2) effArea = 0.0699;
  if(abs(eta)>2.0 && abs(eta)<=2.2 && type == 2) effArea = 0.1056;
  if(abs(eta)>2.2 && abs(eta)<=2.3 && type == 2) effArea = 0.1457;
  if(abs(eta)>2.3 && abs(eta)<=2.4 && type == 2) effArea = 0.1719;
  if(abs(eta)>2.4 && abs(eta)<=2.5 && type == 2) effArea = 0.1998;

  return effArea;
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PhotonUserData);
