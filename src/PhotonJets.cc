#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDMException.h"

// dR and dPhi                                                                                                        

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"

// some basic root definitions
#include <TLorentzVector.h>
#include <vector>
#include <TROOT.h>


// Photon dependencies
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

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

// Fastjet (for creating subjets)                                          
#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>
#include "fastjet/tools/Filter.hh"
#include <fastjet/ClusterSequence.hh>
#include <fastjet/ClusterSequenceArea.hh>
#include <fastjet/tools/Pruner.hh>



// Jets
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"





#include<vector>

using namespace reco;
using namespace edm;
using namespace std;
using namespace fastjet;
using namespace trigger;

typedef std::vector<pat::Jet> PatJetCollection;

class  PhotonJets : public edm::EDProducer {
public:
  PhotonJets( const edm::ParameterSet & );   

private:
  void produce( edm::Event &, const edm::EventSetup & );
  bool isMatchedWithTrigger();

  EDGetTokenT< std::vector< pat::Photon > > phoLabel_;
  EDGetTokenT< std::vector< pat::Jet > > jLabel_;
  EDGetTokenT< std::vector< reco::Vertex > > pvLabel_;
  EDGetTokenT< std::vector< pat::PackedCandidate > > pckPFCdsLabel_;
  EDGetTokenT< double > rhoLabel_;
  edm::EDGetTokenT<edm::ValueMap<bool> > phoLooseIdMapToken_;
  edm::EDGetTokenT<EcalRecHitCollection> ebReducedRecHitCollection_;
  edm::EDGetTokenT<EcalRecHitCollection> eeReducedRecHitCollection_;
  InputTag triggerResultsLabel_, triggerSummaryLabel_;
  InputTag hltElectronFilterLabel_;
  TString hltPath_;
  HLTConfigProvider hltConfig;
  int triggerBit;
  int debug_; 
  
};



//Necessary for the Photon Footprint removal
template <class T, class U>
bool isInFootprint(const T& thefootprint, const U& theCandidate) {
  for ( auto itr = thefootprint.begin(); itr != thefootprint.end(); ++itr ) {

    if( itr.key() == theCandidate.key() ) return true;
    
  }
  return false;
}


PhotonJets::PhotonJets(const edm::ParameterSet& iConfig):
   phoLabel_(consumes<std::vector<pat::Photon>>(iConfig.getParameter<edm::InputTag>("phoLabel"))),
   jLabel_(consumes<std::vector<pat::Jet>>(iConfig.getParameter<edm::InputTag>("jetLabel"))),
   pvLabel_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("pv"))), // "offlinePrimaryVertex"
   pckPFCdsLabel_(consumes<std::vector<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("packedPFCands"))), 
   rhoLabel_(consumes<double>(iConfig.getParameter<edm::InputTag>("rho"))), 
   ebReducedRecHitCollection_(consumes <EcalRecHitCollection> (iConfig.getParameter<edm::InputTag>("ebReducedRecHitCollection"))), //Lazy tool additions
   eeReducedRecHitCollection_(consumes <EcalRecHitCollection> (iConfig.getParameter<edm::InputTag>("eeReducedRecHitCollection")))  // Lazy tool additions  


{
  debug_ = iConfig.getUntrackedParameter<int>("debugLevel",int(0));
  produces<std::vector<pat::Jet> >();
}

void PhotonJets::produce( edm::Event& iEvent, const edm::EventSetup& iSetup) {

  
  //PV
  edm::Handle<std::vector<reco::Vertex> > vertices;
  iEvent.getByToken(pvLabel_, vertices);


  //Photons
  edm::Handle<std::vector<pat::Photon> > phoHandle;
  iEvent.getByToken(phoLabel_, phoHandle);
  std::unique_ptr<vector<pat::Photon> > phoColl( new vector<pat::Photon> (*phoHandle) );
  
  
  //Packed PF Cands
  edm::Handle<std::vector<pat::PackedCandidate>> pfCnd1Handle;
  iEvent.getByToken(pckPFCdsLabel_,pfCnd1Handle); 
  std::unique_ptr<vector<pat::PackedCandidate> > CandColl( new vector<pat::PackedCandidate> (*pfCnd1Handle) );


  //edm::Handle< edm::View<reco::Candidate>> pfCndHandle;
  edm::Handle< edm::View<vector<pat::PackedCandidate>>> pfCndHandle;
  iEvent.getByToken(pckPFCdsLabel_,pfCndHandle);

  //Jet collection
  edm::Handle<pat::JetCollection> fatjets;
  iEvent.getByToken(jLabel_,fatjets);

  //Jet collection
  edm::Handle<std::vector<pat::Jet> > jetHandle, packedjetHandle;
  iEvent.getByToken(jLabel_, jetHandle);
  std::unique_ptr<vector<pat::Jet> > jetColl( new vector<pat::Jet> (*jetHandle) );

  edm::Handle< double > rhoH;
  iEvent.getByToken(rhoLabel_,rhoH);



  if(debug_>=1) cout<<"vtx size " << vertices->size()<<endl; 

  reco::TrackBase::Point vtxPoint(0,0, 0);
  if(  vertices->size() >= 1 ) {
    vtxPoint = vertices->at(0).position();
  }

  if(debug_>=1) cout<<"vtxPoint " <<vtxPoint.X()<<" "<< vtxPoint.Y()<<" "<< vtxPoint.Z()<<endl; 

  // Lets check the photon 

  
  int ijet = -1; 
  for(const pat::Jet &jet : *fatjets) {
    ijet++;
    pat::Jet & jett = (*jetColl)[ijet];
    
    //cout<<" ------------- NEW JET ---------------"<<endl;

    // //cout<<jett.tau1CHS()<<endl;
    
    //Defining necessary variables: 
    float spt0  = -99;     float spt1  = -99;    float spt2  = -99;
    float sphi0 = -99;     float sphi1 = -99;    float sphi2 = -99;
    float seta0 = -99;     float seta1 = -99;    float seta2 = -99;
    float sen0  = -99;     float sen1  = -99;    float sen2  = -99;

    float photon_subjet_frac = -1;
    double jet_eta = jet.eta(); 
    double jet_phi = jet.phi();
    int subindx  = -1;
    int phoindx  = -1;
    float max_g_pt = -99; 
    
    //    if(jet_pt < 200 ) continue; 
    
    //Loop over photons: 
    for(size_t ipho = 0; ipho < (*phoColl).size();ipho++){
	
	pat::Photon & pho = (*phoColl)[ipho];
	double pho_pt  = pho.pt(); 
	double pho_eta = pho.eta(); 
	double pho_phi = pho.phi(); 
	
      if(pho_pt < 20 || fabs(pho_eta) > 2.5 ) continue;
      double dr_pho_jet = deltaR(pho_eta,pho_phi,jet_eta,jet_phi);
      if(dr_pho_jet < 0.8 ){ 
	if(max_g_pt < pho_pt){
	  max_g_pt = pho_pt; 	
	  phoindx  = ipho; 
	}
	
      }//in dr cone      
    }//EOF loop over photons
           
    //in case indeed photon is inside lets play
    if(phoindx != -1 ){
      phoindx = -1; 
      std::vector<reco::Candidate const *> Jetconsts; 
      std::vector<fastjet::PseudoJet> FJparticles;

      
      for (unsigned k = 0; k < jet.numberOfDaughters(); ++k){
	reco::Candidate const * cand = jet.daughter(k);       
	if(cand->numberOfDaughters() == 0){
	  Jetconsts.push_back(cand);
	}else{
	  for(unsigned jda = 0; jda < cand->numberOfDaughters(); ++jda){
	    reco::Candidate const *cand2 = cand->daughter(jda);
	    Jetconsts.push_back(cand2);
	    if( cand2->numberOfDaughters() > 0 ) cout<<"!!!!!!!!!!"<<endl;
	    
	  }
	}//check if they have subdaughters 
      }//Loop over Jet Costituents
      
          
     for(uint ijc = 0; ijc < Jetconsts.size() ; ijc ++){ // Jet Constituens Loop 
	reco::Candidate const * cand = Jetconsts.at(ijc);
	FJparticles.push_back( fastjet::PseudoJet(cand->px(),cand->py(), cand->pz(), cand->energy()) );
     }//Loop of jet constituents
     //Now add the photon in
      //Take the resulting cands and built again the jet , recluster in exlusive 3 subjets

      fastjet::JetDefinition jet_def_kt8(fastjet::kt_algorithm,0.8); 
      fastjet::ClusterSequence clust_seq_08(FJparticles, jet_def_kt8);
      std::vector<fastjet::PseudoJet>  jet_new = sorted_by_pt(clust_seq_08.exclusive_jets_up_to(1));
     	
      Pruner pruner(fastjet::kt_algorithm,0.1,0.5);
      PseudoJet jetnew = pruner(jet_new[0]);  	


      TVector3 subj1,subj2,subj3;
      float subjf1 = -1; 
      float subjf2 = -1; 
      float subjf3 = -1; 
      
      int ipho1 = -1; 
      int ipho2 = -1; 
      int ipho3 = -1; 

      subj1.SetPtEtaPhi(0,0,0);
      subj2.SetPtEtaPhi(0,0,0);
      subj3.SetPtEtaPhi(0,0,0);
    
      if(jetnew.has_pieces()){

        fastjet::PseudoJet sub0; 
        fastjet::PseudoJet sub1; 
        fastjet::PseudoJet suba; 
        fastjet::PseudoJet subb;
        fastjet::PseudoJet subc; 
	      
	jetnew.has_parents(sub0,sub1);
	if(sub0.m() > sub1.m() && sub0.has_pieces() && sub0.m() > 0 ) {
	  sub0.has_parents(suba,subb);
	  if(suba.pt() > 10 && subb.pt() >  10 && sub1.pt() > 10 ){
	    subj1.SetPtEtaPhi(suba.pt(),suba.eta(),suba.phi());
	    subj2.SetPtEtaPhi(subb.pt(),subb.eta(),subb.phi());
	    subj3.SetPtEtaPhi(sub1.pt(),sub1.eta(),sub1.phi());
	  }else if( sub1.has_pieces() && sub1.m() > 0 ){ 
	    sub1.has_parents(suba,subb);
	    if(suba.pt() > 10 && subb.pt() > 10 && sub0.pt() > 10 ){
	      subj1.SetPtEtaPhi(suba.pt(),suba.eta(),suba.phi());
	      subj2.SetPtEtaPhi(subb.pt(),subb.eta(),subb.phi());
	      subj3.SetPtEtaPhi(sub0.pt(),sub0.eta(),sub0.phi());
	    }
	  }
	}else if( sub1.has_pieces() && sub1.m() > 0 && sub1.m() > sub0.m() ){
	   sub1.has_parents(suba,subb);
	   if(suba.pt() >  10 && subb.pt() >  10 && sub0.pt() > 10 ){
	    subj1.SetPtEtaPhi(suba.pt(),suba.eta(),suba.phi());
	    subj2.SetPtEtaPhi(subb.pt(),subb.eta(),subb.phi());
	    subj3.SetPtEtaPhi(sub0.pt(),sub0.eta(),sub0.phi());
	   }else if( sub0.has_pieces() && sub0.m() > 0){ 
	    sub0.has_parents(suba,subb);
	    if(suba.pt() > 10 && subb.pt() > 10 && sub1.pt() > 10 ){
	      subj1.SetPtEtaPhi(suba.pt(),suba.eta(),suba.phi());
	      subj2.SetPtEtaPhi(subb.pt(),subb.eta(),subb.phi());
	      subj3.SetPtEtaPhi(sub1.pt(),sub1.eta(),sub1.phi());
	    }
	   }
	}
	
	if(subj1.Pt() > 10 && subj2.Pt() > 10 && subj3.Pt() > 10  ){
	  for(size_t ipho = 0; ipho < (*phoColl).size();ipho++){
	    pat::Photon & pho = (*phoColl)[ipho];
	    double pho_pt  = pho.pt(); 
	    double pho_eta = pho.eta();
            double pho_hoe = pho.hadTowOverEm();
            double pho_sie = pho.full5x5_sigmaIetaIeta();

	    if(pho_pt < 20  || fabs(pho_eta) > 2.4 || pho_hoe > 0.09 ) continue;
            if(fabs(pho_eta) <1.479 && pho_sie  > 0.0105) continue;
            if(fabs(pho_eta) >1.479 && pho_sie  > 0.0305) continue;

	    double pho_phi = pho.phi(); 
	    TVector3 phovect; 
	    phovect.SetPtEtaPhi(pho_pt,pho_eta,pho_phi);
	    float dr1 = phovect.DeltaR(subj1); 
	    float dr2 = phovect.DeltaR(subj2); 
	    float dr3 = phovect.DeltaR(subj3); 
	    double dr_pho_jet = deltaR(pho_eta,pho_phi,jet_eta,jet_phi);
	    if(dr_pho_jet >  0.8 ) continue; 
	    if(dr1 < dr2 && dr1 <dr3 ){
	      if( (pho_pt/subj1.Pt()) > subjf1 ) subjf1 =  pho_pt / subj1.Pt();
	      ipho1 = ipho; 
	      	    }
	    if(dr2 < dr1 && dr2 <dr3 ){
	      if( (pho_pt/subj2.Pt()) > subjf2 ) subjf2 =  pho_pt / subj2.Pt();
	      ipho2 = ipho; 
	    }
	    if(dr3 < dr2 && dr3 <dr1 ){
	      if( (pho_pt/subj3.Pt()) > subjf3 ) subjf3 =  pho_pt / subj3.Pt();
	      ipho3 = ipho; 
	    }
	  }
	}
      }

      subjf1 =-1; //<<" "<<ppt1 <<endl;
      subjf2 =-1; //<<" "<<ppt2 <<endl;
      subjf3 =-1; //<<" "<<ppt3 <<endl;
      subj1.SetPtEtaPhi(0,0,0);
      subj2.SetPtEtaPhi(0,0,0);
      subj3.SetPtEtaPhi(0,0,0);
      
    


      // Have found probable subjets and probable photons. 
      // New remove the 3-2-1 photons footprints and reclster to get the final subjf restricted [0,1]
      if(ipho1 != -1 || ipho2 != -1 || ipho3 != -1){
	
	std::vector<reco::Candidate const *> Jetconsts; 
	std::vector<fastjet::PseudoJet> FJparticlesUPDATED;
	for (unsigned k = 0; k < jet.numberOfDaughters(); ++k){
	  reco::Candidate const * cand = jet.daughter(k);       
	  if(cand->numberOfDaughters() == 0){
	    Jetconsts.push_back(cand);
	  }else{
	    for(unsigned jda = 0; jda < cand->numberOfDaughters(); ++jda){
	      reco::Candidate const *cand2 = cand->daughter(jda);
	      Jetconsts.push_back(cand2);
	      // DEBUGING FOR PF CANDS IN MINIAOD if( cand2->numberOfDaughters() > 0 ) cout<<"!!!!!!!!!!"<<endl;
	    }
	  }//check if they have subdaughters 
	}//Loop over Jet Costituents
	for(uint ijc = 0; ijc < Jetconsts.size() ; ijc ++){ // Jet Constituens Loop 
	  reco::Candidate const * cand = Jetconsts.at(ijc);
	  for(uint ika = 0; ika < CandColl->size() ;ika++){ //PFCand Loop
	    pat::PackedCandidate & this_constituent = (*CandColl)[ika];
	    const auto& iCand = pfCndHandle->ptrAt(ika);
	    if( this_constituent.px() == cand->px() && this_constituent.py() == cand->py() && this_constituent.eta() == cand->eta() && this_constituent.phi() == cand->phi() &&  this_constituent.energy() == cand->energy()){ 
	      
	      bool isF1 = 0; 
	      bool isF2 = 0; 
	      bool isF3 = 0; 
	      if(ipho1 != -1){
		pat::Photon & pho1 = (*phoColl)[ipho1];
		isF1 = isInFootprint(pho1.associatedPackedPFCandidates(),iCand); 
	      }
	      if(ipho2 != -1){
		pat::Photon & pho2 = (*phoColl)[ipho2];
		isF2 = isInFootprint(pho2.associatedPackedPFCandidates(),iCand); 
	      }
	      if(ipho3 != -1){
		pat::Photon & pho3 = (*phoColl)[ipho3];
		isF3 = isInFootprint(pho3.associatedPackedPFCandidates(),iCand); 
	      }
	      if(isF1 == 0 && isF2 == 0 && isF3 == 0){
		FJparticlesUPDATED.push_back( fastjet::PseudoJet(this_constituent.px(), this_constituent.py(), this_constituent.pz(), this_constituent.energy()) );
	      }
	      break;
	    }//eof check of same pfcand
	  }
	}//Loop of jet constituents
	
	//Now add the photons back in as single objects   
	if(ipho1 != -1 ){
	  pat::Photon & pho = (*phoColl)[ipho1];
	  fastjet::PseudoJet PhoC(pho.px(),pho.py(),pho.pz(),pho.energy());
	  PhoC.set_user_index(99991);
	  FJparticlesUPDATED.push_back( PhoC );
	}
	if(ipho2 != -1 ){
	  pat::Photon & pho = (*phoColl)[ipho2];
	  fastjet::PseudoJet PhoC(pho.px(),pho.py(),pho.pz(),pho.energy());
	  PhoC.set_user_index(99992);
	  FJparticlesUPDATED.push_back( PhoC );
	}
	if(ipho3 != -1 ){
	  pat::Photon & pho = (*phoColl)[ipho3];
	  fastjet::PseudoJet PhoC(pho.px(),pho.py(),pho.pz(),pho.energy());
	  PhoC.set_user_index(99993);
	  FJparticlesUPDATED.push_back( PhoC );
	}
	//Find the new jet
	fastjet::ClusterSequence clust_seq_08_Updated(FJparticlesUPDATED, jet_def_kt8);
	std::vector<fastjet::PseudoJet>  jetUpdated = sorted_by_pt(clust_seq_08_Updated.exclusive_jets_up_to(1));
	
	if(jetUpdated[0].has_pieces()){
          fastjet::PseudoJet sub0;
          fastjet::PseudoJet sub1;
          fastjet::PseudoJet suba;
          fastjet::PseudoJet subb;
          fastjet::PseudoJet subc;

	  jetUpdated[0].has_parents(sub0,sub1);
	  if(sub0.m() > sub1.m() && sub0.has_pieces() && sub0.m() > 0 ) {
	    sub0.has_parents(suba,subb);
	    if(suba.pt() > 10 && subb.pt() >  10 && sub1.pt() > 10 ){
	      subj1.SetPtEtaPhi(suba.pt(),suba.eta(),suba.phi());
	      subj2.SetPtEtaPhi(subb.pt(),subb.eta(),subb.phi());
	      subj3.SetPtEtaPhi(sub1.pt(),sub1.eta(),sub1.phi());
	      subc = sub1;
	    }else if( sub1.has_pieces() && sub1.m() > 0 ){ 
	      sub1.has_parents(suba,subb);
	      if(suba.pt() > 10 && subb.pt() > 10  && sub0.pt() > 10 ){
		subj1.SetPtEtaPhi(suba.pt(),suba.eta(),suba.phi());
		subj2.SetPtEtaPhi(subb.pt(),subb.eta(),subb.phi());
		subj3.SetPtEtaPhi(sub0.pt(),sub0.eta(),sub0.phi());
		subc = sub0;
	      }
	    }
	  }else if( sub1.has_pieces() && sub1.m() > 0 && sub1.m() > sub0.m() ){
	    sub1.has_parents(suba,subb);
	    if(suba.pt() >  10 && subb.pt() >  10 && sub0.pt() > 10 ){
	      subj1.SetPtEtaPhi(suba.pt(),suba.eta(),suba.phi());
	      subj2.SetPtEtaPhi(subb.pt(),subb.eta(),subb.phi());
	      subj3.SetPtEtaPhi(sub0.pt(),sub0.eta(),sub0.phi());
	      subc = sub0;
	    }else if( sub0.has_pieces() && sub0.m() > 0){ 
	      sub0.has_parents(suba,subb);
	      if(suba.pt() > 10 && subb.pt() > 10 && sub1.pt() > 10 ){
		subj1.SetPtEtaPhi(suba.pt(),suba.eta(),suba.phi());
		subj2.SetPtEtaPhi(subb.pt(),subb.eta(),subb.phi());
		subj3.SetPtEtaPhi(sub1.pt(),sub1.eta(),sub1.phi());
		subc = sub1;
	      }
	    }
	  }
	  //Arrange by pt  , 1 is always bigger than 2 here :)	  
	  float spt_1  = subj1.Pt();
	  float spt_2  = subj2.Pt();
	  float spt_3  = subj3.Pt();
	
	  if(spt_2 > spt_3){
	    spt0  = subj1.Pt();   spt1  = subj2.Pt();   spt2  = subj3.Pt();
	    sphi0 = subj1.Phi();  sphi1 = subj2.Phi();  sphi2 = subj3.Phi();
	    seta0 = subj1.Eta();  seta1 = subj2.Eta();  seta2 = subj3.Eta();
	  }
	  
	    
	  if(spt_2 < spt_3 && spt_3 < spt_1){
	    spt0  = subj1.Pt();   spt1  = subj3.Pt();   spt2  = subj2.Pt();
	    sphi0 = subj1.Phi();  sphi1 = subj3.Phi();  sphi2 = subj2.Phi();
	    seta0 = subj1.Eta();  seta1 = subj3.Eta();  seta2 = subj2.Eta();
	  }

	  
	  if(spt_3 >  spt_1 ){
	    spt0  = subj3.Pt();   spt1  = subj1.Pt();   spt2  = subj2.Pt();
	    sphi0 = subj3.Phi();  sphi1 = subj1.Phi();  sphi2 = subj2.Phi();
	    seta0 = subj3.Eta();  seta1 = subj1.Eta();  seta2 = subj2.Eta();
	  }

	
	
	int pho_s1 = -1; 
	int pho_s2 = -1; 
	int pho_s3 = -1; 
	float sbja = -1; 
	float sbjb = -1; 
	float sbjc = -1; 
	
	//cout<<endl;
	//cout<<"finding photons again"<<endl;
	//If indeed photon has convered you Find it as the subjet componen
	if(suba.has_pieces()){
	  vector <fastjet::PseudoJet> Sconst = suba.constituents();
	  for(unsigned int jconst = 0; jconst < Sconst.size(); jconst++){
	    
	    if(Sconst[jconst].user_index() == 99991 && ipho1 != -1){
	      pat::Photon & pho = (*phoColl)[ipho1];
	      if(sbja < pho.pt() / suba.pt() ){
		sbja = (pho.pt() / suba.pt() );
		pho_s1 = ipho1; 
	      }
	    }
	    if(Sconst[jconst].user_index() == 99992 && ipho2 != -1){	    
	      pat::Photon & pho = (*phoColl)[ipho2];
	      if(sbja < pho.pt() / suba.pt() ){
		sbja = (pho.pt() / suba.pt() );
		pho_s1 = ipho2; 
	      }
	    }
	    if(Sconst[jconst].user_index() == 99993 && ipho3 != -1){	    
	      pat::Photon & pho = (*phoColl)[ipho3];
	      if(sbja < pho.pt() / suba.pt() ){
		sbja = (pho.pt() / suba.pt() );
		pho_s1 = ipho3; 
	      }
	    }
	    
	  }
	}
	if(subb.has_pieces() ){
	  vector <fastjet::PseudoJet> Sconst = subb.constituents();
	  for(unsigned int jconst = 0; jconst < Sconst.size(); jconst++){
	    if(Sconst[jconst].user_index() == 99991 && ipho1 != -1){
	      pat::Photon & pho = (*phoColl)[ipho1];
	      if(sbjb < pho.pt() / subb.pt() ){
		sbjb = (pho.pt() / subb.pt() );
		pho_s2 = ipho1; 
	      }
	    }
	    if(Sconst[jconst].user_index() == 99992 && ipho2 != -1){	    
	      pat::Photon & pho = (*phoColl)[ipho2];
	      if(sbjb < pho.pt() / subb.pt() ){
		sbjb = (pho.pt() / subb.pt() );
		pho_s2 = ipho2; 
	      }
	    }
	    if(Sconst[jconst].user_index() == 99993 && ipho3 != -1){	    
	      pat::Photon & pho = (*phoColl)[ipho3];
	      if(sbjb < pho.pt() / subb.pt() ){
		sbjb = (pho.pt() / subb.pt() );
		pho_s2 = ipho3; 
	      }
	    }
	  }
	}
	if(subc.has_pieces()){
	  vector <fastjet::PseudoJet>  Sconst = subc.constituents();	
	  for(unsigned int jconst = 0; jconst < Sconst.size(); jconst++){
	    if(Sconst[jconst].user_index() == 99991 && ipho1 != -1){
	      pat::Photon & pho = (*phoColl)[ipho1];
	      if(sbjc < pho.pt() / subc.pt() ){
		sbjc = (pho.pt() / subc.pt() );
		pho_s3 = ipho1; 
	      }
	    }
	    if(Sconst[jconst].user_index() == 99992 && ipho2 != -1){
	      pat::Photon & pho = (*phoColl)[ipho2];
	      if(sbjc < pho.pt() / subc.pt() ){
		sbjc  = (pho.pt() / subc.pt() );
		pho_s3 = ipho2; 
	      }
	    }
	    if(Sconst[jconst].user_index() == 99993 && ipho3 != -1){
	      pat::Photon & pho = (*phoColl)[ipho3];
	      if(sbjc < pho.pt() / subc.pt() ){
		sbjc  = (pho.pt() / subc.pt() );
		pho_s3 = ipho3; 
	      }
	    }
	  }	  	
	}


	//cout<<"fracs:"<<sbja<<" "<<sbjb<<" "<<sbjc<<endl;

	if(sbja > sbjb && sbja > sbjc) {
	  photon_subjet_frac = sbja;
	  if(suba.pt() > subc.pt()  )   subindx = 1; 
	  if(suba.pt() < subc.pt()  )   subindx = 2; 
	  phoindx = pho_s1; 
	}

	if(sbjb > sbja && sbjb > sbjc) {
	  photon_subjet_frac = sbjb;
	  if(subb.pt() > subc.pt()  )   subindx = 2; 
	  if(suba.pt() < subc.pt()  )   subindx = 3; 
	  if(subb.pt() < subc.pt()  )     subindx = 3;  
	  phoindx = pho_s2; 
 	}

	if(sbjc > sbja && sbjc > sbjb) {
	  photon_subjet_frac = sbjc;
	  if(subb.pt() > subc.pt()  )   subindx = 3; 
	  if(suba.pt() < subc.pt()  )   subindx = 1; 
	  if(subb.pt() < subc.pt() && suba.pt() > subc.pt() ) subindx = 2; 
	  phoindx = pho_s3; 
	}
	  
      }//eof good subjets
    }    

  }//EOF FOUND photon in the jet
    
    jett.addUserInt("phoIndex",phoindx);
    jett.addUserInt("subIndex",subindx);
    jett.addUserFloat("phoSubjetPtFrac",photon_subjet_frac);
        
    jett.addUserFloat("SubPt0", spt0);
    jett.addUserFloat("SubPt1", spt1);
    jett.addUserFloat("SubPt2", spt2);
    jett.addUserFloat("SubEta0",seta0);
    jett.addUserFloat("SubEta1",seta1);
    jett.addUserFloat("SubEta2",seta2); 
    jett.addUserFloat("SubPhi0",sphi0);
    jett.addUserFloat("SubPhi1",sphi1);
    jett.addUserFloat("SubPhi2",sphi2);
    jett.addUserFloat("SubEne0",sen0);
    jett.addUserFloat("SubEne1",sen1);
    jett.addUserFloat("SubEne2",sen2);
  
  }//EOF loop over jets
  iEvent.put( std::move(jetColl) );

}


bool
PhotonJets::isMatchedWithTrigger(){
  return true;
}




#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PhotonJets);
