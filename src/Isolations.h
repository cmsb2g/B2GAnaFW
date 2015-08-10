#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

double getPFMiniIsolation(edm::Handle<pat::PackedCandidateCollection> pfcands,
			  const reco::Candidate* ptcl,  
			  double r_iso_min, double r_iso_max, double kt_scale,
			  bool charged_only);
