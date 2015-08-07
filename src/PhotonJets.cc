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
  bool passIDWP();
  float IsoCalc();

  InputTag phoLabel_;
  InputTag pvLabel_, convLabel_;
  InputTag rhoLabel_;
  InputTag pckPFCdsLabel_;
  InputTag jLabel_;
  edm::EDGetTokenT<EcalRecHitCollection> ebReducedRecHitCollection_;
  edm::EDGetTokenT<EcalRecHitCollection> eeReducedRecHitCollection_;
  // edm::EDGetTokenT<std::vector<pat::Jet> > jetToken_;
  
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





// The Cut based Photon ID's - the current is PHYS 14 first Iteration



PhotonJets::PhotonJets(const edm::ParameterSet& iConfig):
   phoLabel_(iConfig.getParameter<edm::InputTag>("phoLabel")),
   pvLabel_(iConfig.getParameter<edm::InputTag>("pv")),   // "offlinePrimaryVertex"
   rhoLabel_(iConfig.getParameter<edm::InputTag>("rho")), //rhofixedgridRhoFastjet All"
   pckPFCdsLabel_(iConfig.getParameter<edm::InputTag>("packedPFCands")),
   jLabel_       (iConfig.getParameter<edm::InputTag>("jetLabel")),
   ebReducedRecHitCollection_(consumes <EcalRecHitCollection> (iConfig.getParameter<edm::InputTag>("ebReducedRecHitCollection"))), //Lazy tool additions
   eeReducedRecHitCollection_(consumes <EcalRecHitCollection> (iConfig.getParameter<edm::InputTag>("eeReducedRecHitCollection")))  // Lazy tool additions  


{
  debug_ = iConfig.getUntrackedParameter<int>("debugLevel",int(0));
  produces<std::vector<pat::Jet> >();
}

void PhotonJets::produce( edm::Event& iEvent, const edm::EventSetup& iSetup) {

  
  //PV
  edm::Handle<std::vector<reco::Vertex> > vertices;
  iEvent.getByLabel(pvLabel_, vertices);


  //Photons
  edm::Handle<std::vector<pat::Photon> > phoHandle;
  iEvent.getByLabel(phoLabel_, phoHandle);
  auto_ptr<vector<pat::Photon> > phoColl( new vector<pat::Photon> (*phoHandle) );
  
  
  //Packed PF Cands
  edm::Handle<std::vector<pat::PackedCandidate>> pfCnd1Handle;
  iEvent.getByLabel(pckPFCdsLabel_,pfCnd1Handle); 
  auto_ptr<vector<pat::PackedCandidate> > CandColl( new vector<pat::PackedCandidate> (*pfCnd1Handle) );


  edm::Handle< edm::View<reco::Candidate>> pfCndHandle;
  iEvent.getByLabel(pckPFCdsLabel_,pfCndHandle);

  //Jet collection
  edm::Handle<pat::JetCollection> fatjets;
  iEvent.getByLabel(jLabel_,fatjets);

  //Jet collection
  edm::Handle<std::vector<pat::Jet> > jetHandle, packedjetHandle;
  iEvent.getByLabel(jLabel_, jetHandle);
  auto_ptr<vector<pat::Jet> > jetColl( new vector<pat::Jet> (*jetHandle) );

  edm::Handle< double > rhoH;
  iEvent.getByLabel(rhoLabel_,rhoH);



  if(debug_>=1) cout<<"vtx size " << vertices->size()<<endl; 

  reco::TrackBase::Point vtxPoint(0,0, 0);
  if(  vertices->size() >= 1 ) {
    vtxPoint = vertices->at(0).position();
  }

  if(debug_>=1) cout<<"vtxPoint " <<vtxPoint.X()<<" "<< vtxPoint.Y()<<" "<< vtxPoint.Z()<<endl; 
  
  int ijet = -1; 
  for(const pat::Jet &jet : *fatjets) {
    ijet++;
    pat::Jet & jett = (*jetColl)[ijet];
    
    //Defining necessary variables: 
    float spt0  = -99;     float spt1  = -99;    float spt2  = -99;
    float sphi0 = -99;     float sphi1 = -99;    float sphi2 = -99;
    float seta0 = -99;     float seta1 = -99;    float seta2 = -99;
    float sen0  = -99;     float sen1  = -99;    float sen2  = -99;
    
    
    float photon_subjet_frac = 0;
    double jet_eta = jet.eta(); 
    double jet_phi = jet.phi();
    int subindx  = -99;
    int phoindx  = -99;
    float max_g_pt = -99; 
    //Loop over photons: 
    for(size_t ipho = 0; ipho < (*phoColl).size();ipho++){
	
	pat::Photon & pho = (*phoColl)[ipho];
	double pho_pt  = pho.pt(); 
	double pho_eta = pho.eta(); 
	double pho_phi = pho.phi(); 
	
      if(pho_pt < 30 || fabs(pho_eta) > 2.5 ) continue;
      double dr_pho_jet = deltaR(pho_eta,pho_phi,jet_eta,jet_phi);
      if(dr_pho_jet < 0.8 ){ 
	if(max_g_pt < pho_pt){
	  max_g_pt = pho_pt; 	
	  phoindx  = ipho; 
	}
	
      }//in dr cone      
    }//EOF loop over photons


        
    //in case indeed photon is inside lets play
    if(phoindx != -99 ){
      pat::Photon & pho = (*phoColl)[phoindx];
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
	if((pho.associatedPackedPFCandidates()).size() == 0){
	  FJparticles.push_back( fastjet::PseudoJet(cand->px(), cand->py(), cand->pz(), cand->energy()) );
	}else{
	  for(uint ika = 0; ika < CandColl->size() ;ika++){ //PFCand Loop
	    pat::PackedCandidate & this_constituent = (*CandColl)[ika];
	    const auto& iCand = pfCndHandle->ptrAt(ika);
	    if( this_constituent.px() == cand->px() && this_constituent.py() == cand->py() && this_constituent.eta() == cand->eta() && this_constituent.phi() == cand->phi() &&  this_constituent.energy() == cand->energy()){ 
	      bool isF = isInFootprint(pho.associatedPackedPFCandidates(),iCand); 
	      if(isF == 0){
		FJparticles.push_back( fastjet::PseudoJet(this_constituent.px(), this_constituent.py(), this_constituent.pz(), this_constituent.energy()) );
	      }
	      break;
	    }//eof check of same pfcand
	  }
	}
      }//Loop of jet constituents
      
    
            
      if((pho.associatedPackedPFCandidates()).size() > 0){
	fastjet::PseudoJet PhoC(pho.px(),pho.py(),pho.pz(),pho.energy());
	PhoC.set_user_index(0456);
	FJparticles.push_back( PhoC );
      }
      //Take the resulting cands and built again the jet , recluster in exlusive 3 subjets
      // fastjet::JetDefinition jet_def_3(fastjet::antikt_algorithm,1.549); 
      fastjet::JetDefinition jet_def_ca8(fastjet::cambridge_algorithm,0.8); 
      //  fastjet::ClusterSequence clust_seq_3(FJparticles, jet_def_3);
      fastjet::ClusterSequence clust_seq_08(FJparticles, jet_def_ca8);
      
      int nSubJets = 3;
      //std::vector<fastjet::PseudoJet> subjets = sorted_by_pt(clust_seq_3.exclusive_jets_up_to(nSubJets));
      std::vector<fastjet::PseudoJet> subjets = sorted_by_pt(clust_seq_08.exclusive_jets_up_to(nSubJets));
    
  
      if(subjets.size() == 3){ 
	// int sub_wgamma = -99; 
	spt0 = subjets[0].perp();  spt1 = subjets[1].perp();  spt2 = subjets[2].perp();
	sphi0 = subjets[0].phi();  sphi1 = subjets[1].phi();  sphi2 = subjets[2].phi();
	seta0 = subjets[0].eta();  seta1 = subjets[1].eta();  seta2 = subjets[2].eta();
	sen0 = subjets[0].e();     sen1 = subjets[1].e();     sen2 = subjets[2].e();
	
	
	//Looping over reconstructed subjets
	int sindx = -99;
	float psmindr = 99;
	for(int isc = 0 ; isc < nSubJets; isc++){
	  vector<fastjet::PseudoJet> Sconst = subjets[isc].constituents();
	  float subeta = subjets[isc].eta(); 
	  float subphi = subjets[isc].phi();
	  
	  if((pho.associatedPackedPFCandidates()).size() > 0 ){
	    //If indeed photon has convered you Find it as the subjet componen
	    for(unsigned int jconst = 0; jconst < Sconst.size(); jconst++){
	      if(Sconst[jconst].user_index() == 0456){
		sindx = isc;
		break;
	      }
	    }//eof subjet const loop
	  }else{
	    // If it did not disintegrated Make dr matching to subjets
	    float drsp = deltaR(pho.eta(),pho.phi(),subeta,subphi);
	    if(psmindr > drsp ){
	      psmindr = drsp;
	      sindx = isc;
	    }	    
	  }
	}//eof subjets loop  	
	
	if(sindx != -99) {
	  photon_subjet_frac = max_g_pt/subjets[sindx].pt();
	}
      }// If subjet size is 3      
    }//EOF FOUND photon in the jet
    
  
    //store the reconstructed quantities
    jett.addUserInt("jetIndex",ijet);
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
  iEvent.put( jetColl );

}

float
PhotonJets::IsoCalc(){

  return 2.0;

}


bool
PhotonJets::isMatchedWithTrigger(){
  return true;
}


bool PhotonJets::passIDWP(){
 
  return true;
}



#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PhotonJets);
