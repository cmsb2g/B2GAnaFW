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

// MiniIsolation
#include "Isolations.h"

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

class  MuonUserData : public edm::EDProducer {
  public:
    MuonUserData( const edm::ParameterSet & );   

  private:
    void produce( edm::Event &, const edm::EventSetup & );
    bool isMatchedWithTrigger(const pat::Muon, trigger::TriggerObjectCollection,int&,double&,double);
    void put( edm::Event& evt, double value, const char* instanceName);
    float getEA(float eta);
    bool isMediumMuon_23Sep2016(const reco::Muon&);

    TH1F* convertTGraph2TH1F(TGraphAsymmErrors* g);
    double getSF_muonID(double, double);
    double getSFerror_muonID(double, double);
    double getSF_muonISO(double, double);
    double getSFerror_muonISO(double, double);
    double getSF_singleMuonHLT(double, double);
    double getSF_doubleMuonHLT(double, double, double, double);

    EDGetTokenT< std::vector< pat::Muon > > muLabel_;
    EDGetTokenT< std::vector< reco::Vertex > > pvLabel_;
    EDGetTokenT< pat::PackedCandidateCollection > packedPFCandsLabel_;
    EDGetTokenT< edm::TriggerResults > triggerResultsLabel_;
    EDGetTokenT< trigger::TriggerEvent > triggerSummaryLabel_;
    EDGetTokenT< double > rho_miniIso_;
    std::string hltPath_;
    double hlt2reco_deltaRmax_;
    TString mainROOTFILEdir_;
    HLTConfigProvider hltConfig;
    int triggerBit;

    double muon_pt_max_;
    double muon_pt_min_;
    TH1F* muonID_0_0p9_;
    TH1F* muonID_0p9_1p2_;
    TH1F* muonID_1p2_2p1_;
    TH1F* muonID_2p1_2p4_;

    TH1F* muonISO_0_0p9_;
    TH1F* muonISO_0p9_1p2_;
    TH1F* muonISO_1p2_2p1_;
    TH1F* muonISO_2p1_2p4_;

};


MuonUserData::MuonUserData(const edm::ParameterSet& iConfig):
  muLabel_(consumes<std::vector<pat::Muon>>(iConfig.getParameter<edm::InputTag>("muonLabel"))), 
  pvLabel_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("pv"))), // "offlinePrimaryVertex"
  packedPFCandsLabel_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("packedPFCands"))), 
  triggerResultsLabel_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerResults"))),
  triggerSummaryLabel_(consumes<trigger::TriggerEvent>(iConfig.getParameter<edm::InputTag>("triggerSummary"))),
  rho_miniIso_(consumes<double>(edm::InputTag("fixedGridRhoFastjetCentralNeutral"))),
  hltPath_            (iConfig.getParameter<std::string>("hltPath")),
  hlt2reco_deltaRmax_ (iConfig.getParameter<double>("hlt2reco_deltaRmax")),
  mainROOTFILEdir_    (iConfig.getUntrackedParameter<std::string>("mainROOTFILEdir",""))
{
  produces<std::vector<pat::Muon> >();


  if (mainROOTFILEdir_!=""){
    TFile* file_muonSF_ID  = new TFile(mainROOTFILEdir_+"MuonEfficiencies_Run2012ReReco_53X.root",     "READ");

    //DATA_over_MC_Tight_eta_pt20-500
    muonID_0_0p9_   = (file_muonSF_ID->IsZombie() ? NULL : convertTGraph2TH1F( (TGraphAsymmErrors*)file_muonSF_ID->Get("DATA_over_MC_Tight_pt_abseta<0.9") ) );
    muonID_0p9_1p2_ = (file_muonSF_ID->IsZombie() ? NULL : convertTGraph2TH1F( (TGraphAsymmErrors*)file_muonSF_ID->Get("DATA_over_MC_Tight_pt_abseta0.9-1.2") ) );
    muonID_1p2_2p1_ = (file_muonSF_ID->IsZombie() ? NULL : convertTGraph2TH1F( (TGraphAsymmErrors*)file_muonSF_ID->Get("DATA_over_MC_Tight_pt_abseta1.2-2.1") ) );
    muonID_2p1_2p4_ = (file_muonSF_ID->IsZombie() ? NULL : convertTGraph2TH1F( (TGraphAsymmErrors*)file_muonSF_ID->Get("DATA_over_MC_Tight_pt_abseta2.1-2.4") ) );

    TFile* file_muonSF_ISO = new TFile(mainROOTFILEdir_+"MuonEfficiencies_ISO_Run_2012ReReco_53X.root","READ");
    //DATA_over_MC_combRelIsoPF04dBeta<012_Tight_eta_pt20-500
    muonISO_0_0p9_   = (file_muonSF_ISO->IsZombie() ? NULL : convertTGraph2TH1F( (TGraphAsymmErrors*)file_muonSF_ISO->Get("DATA_over_MC_combRelIsoPF04dBeta<012_Tight_pt_abseta<0.9") ) );
    muonISO_0p9_1p2_ = (file_muonSF_ISO->IsZombie() ? NULL : convertTGraph2TH1F( (TGraphAsymmErrors*)file_muonSF_ISO->Get("DATA_over_MC_combRelIsoPF04dBeta<012_Tight_pt_abseta0.9-1.2") ) );
    muonISO_1p2_2p1_ = (file_muonSF_ISO->IsZombie() ? NULL : convertTGraph2TH1F( (TGraphAsymmErrors*)file_muonSF_ISO->Get("DATA_over_MC_combRelIsoPF04dBeta<012_Tight_pt_abseta1.2-2.1") ) );
    muonISO_2p1_2p4_ = (file_muonSF_ISO->IsZombie() ? NULL : convertTGraph2TH1F( (TGraphAsymmErrors*)file_muonSF_ISO->Get("DATA_over_MC_combRelIsoPF04dBeta<012_Tight_pt_abseta2.1-2.4") ) );

  }

  //  TFile* file_muonSF_singleMuHLT = new TFile(mainROOTFILEdir_+"SingleMuonTriggerEfficiencies_eta2p1_Run2012ABCD_v5trees.root","READ");
  //  TFile* file_muonSF_doubleMuHLT = new TFile(mainROOTFILEdir_+"MuHLTEfficiencies_Run_2012ABCD_53X_DR03-2","READ");

  if (mainROOTFILEdir_!=""){
    if (muonID_0_0p9_!=NULL) {
      muon_pt_min_ = muonID_0_0p9_->GetXaxis()->GetXmin();
      muon_pt_max_ = muonID_0_0p9_->GetXaxis()->GetXmax();
    }
  }
}


#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
void MuonUserData::produce( edm::Event& iEvent, const edm::EventSetup& iSetup) {


  //  bool isMC = (!iEvent.isRealData());

  //PV
  edm::Handle<std::vector<reco::Vertex> > pvHandle;
  iEvent.getByToken(pvLabel_, pvHandle);
  const reco::Vertex& PV= pvHandle->front();

  //RHO
  edm::Handle<double> rhoHandle_miniIso;
  iEvent.getByToken(rho_miniIso_,rhoHandle_miniIso);
  double rho_miniIso = *rhoHandle_miniIso;

  //Muons
  edm::Handle<std::vector<pat::Muon> > muonHandle;
  iEvent.getByToken(muLabel_, muonHandle);
  std::unique_ptr<vector<pat::Muon> > muonColl( new vector<pat::Muon> (*muonHandle) );

  //PackedPFCands for Mini-isolation
  edm::Handle<pat::PackedCandidateCollection> packedPFCands;
  iEvent.getByToken(packedPFCandsLabel_, packedPFCands);

  ///////////  /////////  /////////  /////////  /////////  /////////  /////////  /////////  /////////
  //// TRIGGER (this is not really needed ...)
  //bool changedConfig = false;
  //bool pathFound = false;
  //if (!hltConfig.init(iEvent.getRun(), iSetup, "HLT", changedConfig)) {
  //  std::cout << ">>>> MuonUserData: Initialization of HLTConfigProvider failed!!" << std::endl;
  //  return;
  //}

  //if (changedConfig){
  //  std::cout << ">>>> MuonUserData: The current menu is " << hltConfig.tableName() << std::endl;
  //  triggerBit = -1;
  //  for (size_t j = 0; j < hltConfig.triggerNames().size(); j++) {
  //    //      std::cout << "hltConfig.triggerNames()[" << j << "]: " << hltConfig.triggerNames()[j] << std::endl;
  //    if (TString(hltConfig.triggerNames()[j]).Contains(hltPath_)) {triggerBit = j;pathFound=true;}
  //  }
  //  if (triggerBit == -1) std::cout << ">>>> MuonUserData: HLT path not found" << std::endl;
  //}

  //edm::Handle<edm::TriggerResults> triggerResults;
  //iEvent.getByToken(triggerResultsLabel_, triggerResults);
  //if (size_t(triggerBit) < triggerResults->size() && pathFound  )
  //  if (triggerResults->accept(triggerBit))
  //    std::cout << ">>>> MuonUserData: event pass : " << hltPath_ << std::endl;

  /////////  /////////  /////////  /////////  /////////  /////////  /////////  /////////  /////////
  // TRIGGER MATCHING

  trigger::TriggerObjectCollection MuonLegObjects;

  //edm::Handle<trigger::TriggerEvent> triggerSummary;

  //if ( triggerSummary.isValid() ) {
  //  iEvent.getByToken(triggerSummaryLabel_, triggerSummary);

  //  // Results from TriggerEvent product - Attention: must look only for
  //  // modules actually run in this path for this event!
  //  if(pathFound){
  //    const unsigned int triggerIndex(hltConfig.triggerIndex(hltPath_));
  //    const vector<string>& moduleLabels(hltConfig.moduleLabels(triggerIndex));
  //    const unsigned int moduleIndex(triggerResults->index(triggerIndex));

  //    for (unsigned int j=0; j<=moduleIndex; ++j) {

  //      const string& moduleLabel(moduleLabels[j]);
  //      const string  moduleType(hltConfig.moduleType(moduleLabel));
  //      // check whether the module is packed up in TriggerEvent product
  //      const unsigned int filterIndex(triggerSummary->filterIndex(InputTag(moduleLabel,"","HLT")));

  //      if (filterIndex<triggerSummary->sizeFilters()) {
  //        //	  cout << " 'L3' filter in slot " << j << " - label/type " << moduleLabel << "/" << moduleType << endl;
  //        TString lable = moduleLabel.c_str();
  //        if (lable.Contains(hltMuonFilterLabel_.label())) {

  //          const trigger::Vids& VIDS (triggerSummary->filterIds(filterIndex));
  //          const trigger::Keys& KEYS(triggerSummary->filterKeys(filterIndex));
  //          const size_type nI(VIDS.size());
  //          const size_type nK(KEYS.size());
  //          assert(nI==nK);
  //          const size_type n(max(nI,nK));
  //          //	    cout << "   " << n  << " accepted TRIGGER objects found: " << endl;
  //          const trigger::TriggerObjectCollection& TOC(triggerSummary->getObjects());
  //          for (size_type i=0; i!=n; ++i) {
  //            const trigger::TriggerObject& TO(TOC[KEYS[i]]);
  //            MuonLegObjects.push_back(TO);	  
  //            //	  cout << "   " << i << " " << VIDS[i] << "/" << KEYS[i] << ": "
  //            //	       << TO.id() << " " << TO.pt() << " " << TO.eta() << " " << TO.phi() << " " << TO.mass()
  //            //	       << endl;
  //          }
  //        }
  //      }
  //    }
  //  }
  //}

  //  std::cout << "----> MuonLegObjects: " << MuonLegObjects.size() << " <--> RECO : " << muonColl->size() << std::endl;
  /////////  /////////  /////////  /////////  /////////  /////////  /////////  /////////  /////////  /////////
  for (size_t i = 0; i< muonColl->size(); i++){
    pat::Muon & m = (*muonColl)[i];

    // muon ID
    bool isSoftMuon  = m.isSoftMuon(PV);
    bool isLooseMuon = m.isLooseMuon();
    bool isMediumMuon = m.isMediumMuon();
    bool isTightMuon = m.isTightMuon(PV);
    bool isHighPtMuon = m.isHighPtMuon(PV);

    // 2016 Medium ID
    // https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2?rev=27#MediumID2016_to_be_used_with_Run
    bool isMediumMuon2016 = isMediumMuon_23Sep2016(m);
    double segmentCompatibility = muon::segmentCompatibility(m);

    // impact parameters
    double dxy = m.muonBestTrack()->dxy(PV.position());     
    double dxyErr = m.muonBestTrack()->dxyError();     
    double dz  = m.muonBestTrack()->dz(PV.position());
    double dzErr  = m.muonBestTrack()->dzError();
    double dB    = m.dB (pat::Muon::PV3D);
    double dBErr = m.edB(pat::Muon::PV3D);

    // isolation (delta beta corrections)
    // https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId#Accessing_PF_Isolation_from_reco
    double sumChargedHadronPt = m.pfIsolationR04().sumChargedHadronPt;
    double sumNeutralHadronPt = m.pfIsolationR04().sumNeutralHadronEt;
    double sumPhotonPt        = m.pfIsolationR04().sumPhotonEt;
    double sumPUPt            = m.pfIsolationR04().sumPUPt;
    double pt                 = m.pt();
    double iso04 = (sumChargedHadronPt+TMath::Max(0.,sumNeutralHadronPt+sumPhotonPt-0.5*sumPUPt))/pt;
    double EA = getEA(m.eta());
    double miniIso = getPFMiniIsolation(packedPFCands, dynamic_cast<const reco::Candidate *>(&m), 0.05, 0.2, 10., false, true, EA, rho_miniIso);

    // trigger matched 

    //int idx       = -1;
    //double deltaR = -1.;
    //bool isMatched2trigger = isMatchedWithTrigger(m, MuonLegObjects, idx, deltaR, hlt2reco_deltaRmax_) ;
    //double hltEta = ( isMatched2trigger ? MuonLegObjects[0].eta()    : -999.);
    //double hltPhi = ( isMatched2trigger ? MuonLegObjects[0].phi()    : -999.);
    //double hltPt  = ( isMatched2trigger ? MuonLegObjects[0].pt()     : -999.);
    //double hltE   = ( isMatched2trigger ? MuonLegObjects[0].energy() : -999.);

    // SF from Muon POG
    // https://twiki.cern.ch/twiki/bin/viewauth/CMS/MuonReferenceEffs
    double eta = m.eta();
    double muonIDweight  = getSF_muonID (pt,eta);
    double muonISOweight = getSF_muonISO(pt,eta);
    double muonHLTweight = getSF_singleMuonHLT(pt,eta);


    m.addUserFloat("isSoftMuon",  isSoftMuon);
    m.addUserFloat("isLooseMuon", isLooseMuon);
    m.addUserFloat("isMediumMuon", isMediumMuon);
    m.addUserFloat("isMediumMuon2016", isMediumMuon2016);
    m.addUserFloat("isTightMuon", isTightMuon);
    m.addUserFloat("isHighPtMuon", isHighPtMuon);
    m.addUserFloat("segmentCompatibility", segmentCompatibility);
    m.addUserFloat("dxy",         dxy);
    m.addUserFloat("dxyErr",      dxyErr);
    m.addUserFloat("dz",          dz);
    m.addUserFloat("dzErr",       dzErr);
    m.addUserFloat("dB",          dB);
    m.addUserFloat("dBErr",       dBErr);
    m.addUserFloat("iso04",       iso04);
    m.addUserFloat("miniIso",     miniIso);

    //m.addUserFloat("HLTmuonEta",   hltEta);
    //m.addUserFloat("HLTmuonPhi",   hltPhi);
    //m.addUserFloat("HLTmuonPt",    hltPt);
    //m.addUserFloat("HLTmuonE",     hltE);
    //m.addUserFloat("HLTmuonDeltaR",deltaR);

    m.addUserFloat("muonIDweight",  muonIDweight );
    m.addUserFloat("muonISOweight", muonISOweight);
    m.addUserFloat("muonHLTweight", muonHLTweight);

    /*
    //     **** DEBUG dB ****
    double dB_     = m.dB ();
    double dB_PV2D = m.dB (pat::Muon::PV2D);
    double dB_PV3D = m.dB (pat::Muon::PV3D);
    double dB_BS2D = m.dB (pat::Muon::BS2D);
    double dB_BS3D = m.dB (pat::Muon::BS3D);

    m.addUserFloat("dB",    dB_);
    m.addUserFloat("dBPV2D",dB_PV2D);
    m.addUserFloat("dBPV3D",dB_PV3D);
    m.addUserFloat("dBBS2D",dB_BS2D);
    m.addUserFloat("dBBS3D",dB_BS3D);
    */

  }

  iEvent.put( std::move(muonColl) );

}

// ------------ method called once each job just after ending the event loop  ------------

  bool
MuonUserData::isMatchedWithTrigger(const pat::Muon p, trigger::TriggerObjectCollection triggerObjects, int& index, double& deltaR, double deltaRmax = 0.1)
{
  for (size_t i = 0 ; i < triggerObjects.size() ; i++){
    //float dR = sqrt(pow(triggerObjects[i].eta()-p.eta(),2)+ pow(acos(cos(triggerObjects[i].phi()-p.phi())),2)) ;

    float dR = 0;
    float DeltaPhi = 0;

    if((triggerObjects[i].phi()-p.phi()) > acos(-1)) {
      DeltaPhi = pow((triggerObjects[i].phi() - p.phi())- 2*acos(-1),2);
    }
    else if((triggerObjects[i].phi()-p.phi()) <= -acos(-1)) {
      DeltaPhi = pow((triggerObjects[i].phi() - p.phi()) + 2*acos(-1),2);
    }
    else{
      DeltaPhi = pow((triggerObjects[i].phi()-p.phi()),2);
    }

    dR=sqrt(pow(triggerObjects[i].eta()-p.eta(),2)+ pow(DeltaPhi,2)) ;

    if (dR<deltaRmax) {
      deltaR = dR;
      index  = i;
      return true;
    }
  }
  return false;
}

  void 
MuonUserData::put(edm::Event& evt, double value, const char* instanceName)
{
  std::unique_ptr<double> varPtr(new double(value));
  evt.put(std::move(varPtr), instanceName);
}

  double
MuonUserData::getSF_muonID(double pt, double eta)
{  
  double SF = 1.;
  if (mainROOTFILEdir_!=""){
    if(pt < muon_pt_min_) pt = muon_pt_min_;
    if(pt > muon_pt_max_) pt = muon_pt_max_;
    eta=fabs(eta);

    if(eta<=0.9) SF = muonID_0_0p9_->GetBinContent(muonID_0_0p9_->FindBin(pt));
    if(eta>0.9 && eta<=1.2) SF = muonID_0p9_1p2_->GetBinContent(muonID_0p9_1p2_->FindBin(pt));
    if(eta>1.2 && eta<=2.1) SF = muonID_1p2_2p1_->GetBinContent(muonID_1p2_2p1_->FindBin(pt));
    if(eta>2.1 && eta<=2.4) SF = muonID_2p1_2p4_->GetBinContent(muonID_2p1_2p4_->FindBin(pt));
  }
  return SF;
}
  double
MuonUserData::getSFerror_muonID(double pt, double eta)
{  
  double SF = 1.;
  if(pt < muon_pt_min_) pt = muon_pt_min_;
  if(pt > muon_pt_max_) pt = muon_pt_max_;
  eta=fabs(eta);

  if (mainROOTFILEdir_!=""){
    if(eta<=0.9) SF = muonID_0_0p9_->GetBinError(muonID_0_0p9_->FindBin(pt));
    if(eta>0.9 && eta<=1.2) SF = muonID_0p9_1p2_->GetBinError(muonID_0p9_1p2_->FindBin(pt));
    if(eta>1.2 && eta<=2.1) SF = muonID_1p2_2p1_->GetBinError(muonID_1p2_2p1_->FindBin(pt));
    if(eta>2.1 && eta<=2.4) SF = muonID_2p1_2p4_->GetBinError(muonID_2p1_2p4_->FindBin(pt));
  }
  return SF;
}

  double
MuonUserData::getSF_muonISO(double pt, double eta)
{
  double SF = 1.;
  if(pt < muon_pt_min_) pt = muon_pt_min_;
  if(pt > muon_pt_max_) pt = muon_pt_max_;
  eta=fabs(eta);

  if (mainROOTFILEdir_!=""){
    if(eta<=0.9) SF = muonISO_0_0p9_->GetBinContent(muonISO_0_0p9_->FindBin(pt));
    if(eta>0.9 && eta<=1.2) SF = muonISO_0p9_1p2_->GetBinContent(muonISO_0p9_1p2_->FindBin(pt));
    if(eta>1.2 && eta<=2.1) SF = muonISO_1p2_2p1_->GetBinContent(muonISO_1p2_2p1_->FindBin(pt));
    if(eta>2.1 && eta<=2.4) SF = muonISO_2p1_2p4_->GetBinContent(muonISO_2p1_2p4_->FindBin(pt));
  }
  return SF;  
}
  double
MuonUserData::getSFerror_muonISO(double pt, double eta)
{
  double SF = 1.;
  if(pt < muon_pt_min_) pt = muon_pt_min_;
  if(pt > muon_pt_max_) pt = muon_pt_max_;
  eta=fabs(eta);

  if (mainROOTFILEdir_!=""){
    if(eta<=0.9) SF = muonISO_0_0p9_->GetBinError(muonISO_0_0p9_->FindBin(pt));
    if(eta>0.9 && eta<=1.2) SF = muonISO_0p9_1p2_->GetBinError(muonISO_0p9_1p2_->FindBin(pt));
    if(eta>1.2 && eta<=2.1) SF = muonISO_1p2_2p1_->GetBinError(muonISO_1p2_2p1_->FindBin(pt));
    if(eta>2.1 && eta<=2.4) SF = muonISO_2p1_2p4_->GetBinError(muonISO_2p1_2p4_->FindBin(pt));
  }
  return SF;  
}

  double
MuonUserData::getSF_singleMuonHLT(double pt, double eta)
{
  double SF = 1.;
  return SF;  
}

  double
MuonUserData::getSF_doubleMuonHLT(double pt1, double eta1, double pt2, double eta2)
{
  double SF = 1.;
  return SF;  
}

TH1F*
MuonUserData::convertTGraph2TH1F(TGraphAsymmErrors* g) {
  size_t n=g->GetN();
  float x[n+1];
  for(size_t i=0; i<n; i++) 
    x[i]=g->GetX()[i]-g->GetEXlow()[i];
  x[n]=g->GetX()[n-1]+g->GetEXhigh()[n-1];

  TH1F* h=new TH1F(g->GetName(), g->GetTitle(), n, x);
  h->Sumw2();
  for(size_t i=0; i<n; i++) {
    h->SetBinContent(i+1, g->GetY()[i]);
    h->SetBinError(i+1, g->GetEYhigh()[i]);
  }
  return h;
}

float MuonUserData::getEA(float eta){
  float effArea = 0.;
  if(abs(eta)>0.0 && abs(eta)<=0.8) effArea = 0.0735;
  if(abs(eta)>0.8 && abs(eta)<=1.3) effArea = 0.0619;
  if(abs(eta)>1.3 && abs(eta)<=2.0) effArea = 0.0465;
  if(abs(eta)>2.0 && abs(eta)<=2.2) effArea = 0.0433;
  if(abs(eta)>2.2 && abs(eta)<=2.5) effArea = 0.0577;
  return effArea;
}

bool MuonUserData::isMediumMuon_23Sep2016(const reco::Muon & recoMu) {
  bool goodGlob = recoMu.isGlobalMuon() && 
    recoMu.globalTrack()->normalizedChi2() < 3 && 
    recoMu.combinedQuality().chi2LocalPosition < 12 && 
    recoMu.combinedQuality().trkKink < 20; 
  bool isMedium = muon::isLooseMuon(recoMu) && 
    recoMu.innerTrack()->validFraction() > 0.49 && 
    muon::segmentCompatibility(recoMu) > (goodGlob ? 0.303 : 0.451); 
  return isMedium; 
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(MuonUserData);
