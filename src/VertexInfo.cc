/*
 *\Author: A. Orso M. Iorio 
 *
 *
 *\version  $Id: SingleTopVertexInfoDumper.cc,v 1.1.2.1 2011/09/21 13:19:39 oiorio Exp $ 
 */

// Single Top producer: produces a top candidate made out of a Lepton, a B jet and a MET
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include <FWCore/Framework/interface/Run.h>

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h" 

#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
//#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Common/interface/View.h"

#include "DataFormats/Candidate/interface/NamedCompositeCandidate.h"

//#include "TLorentzVector.h"
//#include "TopQuarkAnalysis/SingleTop/interface/EquationSolver.h"

#include "DataFormats/Candidate/interface/CandAssociation.h"

//#include "TopQuarkAnalysis/SingleTop/interface/VertexInfo.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include <vector>
#include <memory>

#include "DataFormats/Math/interface/LorentzVector.h"


//using namespace pat;

class VertexInfo : public edm::EDProducer {
public:
  explicit VertexInfo(const edm::ParameterSet & iConfig);
  ~VertexInfo();
  virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup);
  //       static void fillDescriptions(edm::ConfigurationDescriptions & descriptions);
private:
  edm::EDGetTokenT< std::vector< reco::Vertex > > src_;
  
};
//}

VertexInfo::VertexInfo(const edm::ParameterSet& iConfig):
   src_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("src")))
{
  // initialize the configurables
  produces<int>("npv");
  produces<std::vector<float>>("vx");
  produces<std::vector<float>>("vy");
  produces<std::vector<float>>("vz"); 
  produces< std::vector<float> >("chi");
  produces< std::vector<float> >("rho");
  produces< std::vector<float> >("z");
  produces< std::vector< int> >("ndof");
 
}

void VertexInfo::produce(edm::Event & iEvent, const edm::EventSetup & iEventSetup){
  
  
  edm::Handle<std::vector<reco::Vertex> > vertices;
  iEvent.getByToken(src_, vertices);

  
  // primary vertices

  std::unique_ptr<int> npv (new int() );
  std::unique_ptr<std::vector<float>> vx (new std::vector<float> );
  std::unique_ptr<std::vector<float>> vy (new std::vector<float> );
  std::unique_ptr<std::vector<float>> vz (new std::vector<float> );

  std::unique_ptr< std::vector< float > >chi_(new std::vector< float >) ;
  std::unique_ptr< std::vector< float > >rho_ (new std::vector< float >) ;
  std::unique_ptr< std::vector< float > >z_ (new std::vector< float >) ;
  std::unique_ptr< std::vector< int > >ndof_ (new std::vector< int >) ;

  *npv = vertices->size();

  for(size_t v = 0; v<vertices->size();++v){;
    vx -> push_back(vertices->at(v).x()) ;
    vy -> push_back(vertices->at(v).y()) ;
    vz -> push_back(vertices->at(v).z()) ;
    chi_->push_back(vertices->at(v).chi2());
    ndof_->push_back(vertices->at(v).ndof());
    z_->push_back(vertices->at(v).position().z());
    rho_->push_back(vertices->at(v).position().rho());

  }


  iEvent.put( std::move(npv), "npv"); 
  iEvent.put( std::move(vx), "vx"); 
  iEvent.put( std::move(vy), "vy"); 
  iEvent.put( std::move(vz), "vz"); 

  iEvent.put( std::move(chi_),"chi");
  iEvent.put( std::move(rho_),"rho");
  iEvent.put( std::move(z_),"z");
  iEvent.put( std::move(ndof_),"ndof");
}

VertexInfo::~VertexInfo(){;}


DEFINE_FWK_MODULE( VertexInfo );
