#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDMException.h"

#include "DataFormats/PatCandidates/interface/Jet.h"

// dR and dPhi
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"

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

#include <vector>
#include <sstream>

using namespace reco;
using namespace edm;
using namespace std;
using namespace trigger;

typedef std::vector<pat::Jet> PatJetCollection;

struct orderByPt {
  const std::string mCorrLevel;
  orderByPt(const std::string& fCorrLevel) : mCorrLevel(fCorrLevel) {}
  bool operator ()(PatJetCollection::const_iterator const& a, PatJetCollection::const_iterator const& b) {
    if( mCorrLevel=="Uncorrected" ) return a->correctedJet("Uncorrected").pt() > b->correctedJet("Uncorrected").pt();
    else return a->pt() > b->pt();
  }
};


class PatJetUserData : public edm::EDProducer {
  public:
    PatJetUserData( const edm::ParameterSet & );

  private:
    void produce( edm::Event &, const edm::EventSetup & );
    bool isMatchedWithTrigger(const pat::Jet&, trigger::TriggerObjectCollection,int&,double&,double);
    double getResolutionRatio(double eta);
    double getJERup(double eta);
    double getJERdown(double eta);

    void matchPackedJets(const edm::Handle<PatJetCollection>& jets,
        const edm::Handle<PatJetCollection>& packedjets,
        std::vector<int>& matchedIndices) ;

    edm::EDGetTokenT<std::vector<pat::Jet> >     jetToken_;
    edm::EDGetTokenT<std::vector<reco::Vertex> > pvToken_;

    //InputTag jetLabel_;
    InputTag jLabel_, packedjLabel_, sjLabel_;
    InputTag pvLabel_;
    InputTag triggerResultsLabel_, triggerSummaryLabel_;
    InputTag hltJetFilterLabel_;
    std::string hltPath_;
    double hlt2reco_deltaRmax_;
    HLTConfigProvider hltConfig;
    int triggerBit;
    bool doSubjets_ ;
};


PatJetUserData::PatJetUserData(const edm::ParameterSet& iConfig) :

  jLabel_(iConfig.getParameter<edm::InputTag>("jetLabel")),
  packedjLabel_(iConfig.getParameter<edm::InputTag>("packedjetLabel")),
  sjLabel_(iConfig.getParameter<edm::InputTag>("subjetLabel")),
  pvLabel_(iConfig.getParameter<edm::InputTag>("pv")),   // "offlinePrimaryVertex"

  triggerResultsLabel_(iConfig.getParameter<edm::InputTag>("triggerResults")),
  triggerSummaryLabel_(iConfig.getParameter<edm::InputTag>("triggerSummary")),
  hltJetFilterLabel_  (iConfig.getParameter<edm::InputTag>("hltJetFilter")),   //trigger objects we want to match
  hltPath_            (iConfig.getParameter<std::string>("hltPath")),
  hlt2reco_deltaRmax_ (iConfig.getParameter<double>("hlt2reco_deltaRmax")),
  doSubjets_          (iConfig.getParameter<bool>("doSubjets"))
{
  produces<vector<pat::Jet> >();
}

void PatJetUserData::produce( edm::Event& iEvent, const edm::EventSetup& iSetup) {

  //Jets
  edm::Handle<std::vector<pat::Jet> > jetHandle, packedjetHandle, subjetHandle;
  iEvent.getByLabel(jLabel_, jetHandle);

  auto_ptr<vector<pat::Jet> > jetColl( new vector<pat::Jet> (*jetHandle) );

  int subjet0Id = -1, subjet1Id = -1 ;
  std::vector<int> matchedjetIndices;
  if (doSubjets_) {
    iEvent.getByLabel(packedjLabel_, packedjetHandle);
    iEvent.getByLabel(sjLabel_, subjetHandle);
    if( packedjetHandle->size() > jetHandle->size() )
      edm::LogError("TooManyGroomedJets") << "Groomed packed jet collection size " << packedjetHandle->size()
        << " bigger than original fat jet collection size " << jetHandle->size() ;
    matchPackedJets(jetHandle,packedjetHandle,matchedjetIndices);

    for (size_t ii = 0; ii < jetColl->size(); ii++){
      pat::Jet*  jet = &(jetColl->at(ii)) ;
      int pfjIdx = matchedjetIndices.at(ii) ; //// Get index of matched packed jet
      int nSubjets = 0 ;
      std::vector<PatJetCollection::const_iterator> itsubjets ;
      if (pfjIdx >= 0) {
        nSubjets = packedjetHandle->at(pfjIdx).numberOfDaughters() ;
        const edm::Ptr<reco::Jet> originalObjRef = edm::Ptr<reco::Jet>( packedjetHandle->at(pfjIdx).originalObjectRef() );
        for (int isj = 0; isj < nSubjets; ++isj) {
          for (std::vector<pat::Jet>::const_iterator itsj = subjetHandle->begin(); itsj != subjetHandle->end(); ++itsj) {
            if ( originalObjRef->daughterPtr(isj) == itsj->originalObjectRef() ) {
              itsubjets.push_back(itsj) ;
            }
          }
        } //// Loop over subjets
      }
      std::sort(itsubjets.begin(), itsubjets.end(), orderByPt("Uncorrected"));
      if ( itsubjets.size() > 1 ) {
        subjet0Id = itsubjets.at(0) - subjetHandle->begin() ;
        subjet1Id = itsubjets.at(1) - subjetHandle->begin() ;
        jet->addUserInt("subjetIndex0", subjet0Id) ;
        jet->addUserInt("subjetIndex1", subjet1Id) ;
        std::cout << " 2 subjet found with indices " << subjet0Id << " and " << subjet1Id << "\n" ;
      }
    } //// Looping over all jets

  }

  iEvent.put( jetColl );

}

// ------------ method called once each job just after ending the event loop  ------------
  bool
PatJetUserData::isMatchedWithTrigger(const pat::Jet& p, trigger::TriggerObjectCollection triggerObjects, int& index, double& deltaR, double deltaRmax = 0.2)
{
  for (size_t i = 0 ; i < triggerObjects.size() ; i++){
    float dR = sqrt(pow(triggerObjects[i].eta()-p.eta(),2)+ pow(acos(cos(triggerObjects[i].phi()-p.phi())),2)) ;
    //    std::cout << "dR: " << dR << std::endl;
    if (dR<deltaRmax) {
      deltaR = dR;
      index  = i;
      return true;
    }
  }
  return false;
}

  double
PatJetUserData::getResolutionRatio(double eta)
{
  eta=fabs(eta);
  if(eta>=0.0 && eta<0.5) return 1.079; // +-0.005 +-0.026
  if(eta>=0.5 && eta<1.1) return 1.099; // +-0.005 +-0.028
  if(eta>=1.1 && eta<1.7) return 1.121; // +-0.005 +-0.029
  if(eta>=1.7 && eta<2.3) return 1.208; // +-0.013 +-0.045
  if(eta>=2.3 && eta<2.8) return 1.254; // +-0.026 +-0.056
  if(eta>=2.8 && eta<3.2) return 1.395; // +-0.036 +-0.051
  if(eta>=3.2 && eta<5.0) return 1.056; // +-0.048 +-0.185
  return -1.;
}

  double
PatJetUserData::getJERup(double eta)
{
  eta=fabs(eta);
  if(eta>=0.0 && eta<0.5) return 1.053 ;
  if(eta>=0.5 && eta<1.1) return 1.071 ;
  if(eta>=1.1 && eta<1.7) return 1.092 ;
  if(eta>=1.7 && eta<2.3) return 1.162 ;
  if(eta>=2.3 && eta<2.8) return 1.192 ;
  if(eta>=2.8 && eta<3.2) return 1.332 ;
  if(eta>=3.2 && eta<5.0) return 0.865 ;
  return -1.;
}

  double
PatJetUserData::getJERdown(double eta)
{
  eta=fabs(eta);
  if(eta>=0.0 && eta<0.5) return 1.105 ;
  if(eta>=0.5 && eta<1.1) return 1.127 ;
  if(eta>=1.1 && eta<1.7) return 1.150 ;
  if(eta>=1.7 && eta<2.3) return 1.254 ;
  if(eta>=2.3 && eta<2.8) return 1.316 ;
  if(eta>=2.8 && eta<3.2) return 1.458 ;
  if(eta>=3.2 && eta<5.0) return 1.247 ;
  return -1.;
}

// ------------ method that matches packed and original jets based on minimum dR ------------
void PatJetUserData::matchPackedJets(const edm::Handle<PatJetCollection>& jets,
    const edm::Handle<PatJetCollection>& packedJets, std::vector<int>& matchedIndices) {

  std::vector<bool> jetLocks(jets->size(),false);
  std::vector<int>  jetIndices;

  for(size_t pj=0; pj<packedJets->size(); ++pj) {
    double matchedDR = 1e9;
    int matchedIdx = -1;

    for(size_t ij = 0; ij < jets->size(); ++ij) {
      if( jetLocks.at(ij) ) continue; // skip jets that have already been matched

      double tempDR = reco::deltaR( jets->at(ij).rapidity(), jets->at(ij).phi(), packedJets->at(pj).rapidity(), packedJets->at(pj).phi() );
      if( tempDR < matchedDR ) {
        matchedDR = tempDR;
        matchedIdx = ij;
      }
    }

    if( matchedIdx >= 0 ) jetLocks.at(matchedIdx) = true;
    jetIndices.push_back(matchedIdx);
  }

  if( std::find( jetIndices.begin(), jetIndices.end(), -1 ) != jetIndices.end() )
    edm::LogError("JetMatchingFailed") << "Matching groomed to original jets failed. Please check that the two jet collections belong to each other.";

  for(size_t ij = 0; ij < jets->size(); ++ij) {
    std::vector<int>::iterator matchedIndex = std::find( jetIndices.begin(), jetIndices.end(), ij );
    matchedIndices.push_back( matchedIndex != jetIndices.end() ? std::distance(jetIndices.begin(),matchedIndex) : -1 );
  }
}

#include "FWCore/Framework/interface/MakerMacros.h"


DEFINE_FWK_MODULE(PatJetUserData);
