//  Njettiness Package
//  Version 0.5.1 (September 19, 2012)
//  Questions/Comments?  jthaler@jthaler.net

// Copyright (c) 2011-12, Jesse Thaler, Ken Van Tilburg, and Christopher K.
// Vermilion
//
//----------------------------------------------------------------------
// This file is part of the N-jettiness package ("N-jettiness").
//
//  N-jettiness is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 3 of the License, or
//  (at your option) any later version.
//
//  SpartyJet is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with SpartyJet; if not, write to the Free Software
//  Foundation, Inc.:
//      59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//----------------------------------------------------------------------


#ifndef __NSUBJETTINESS_HH__
#define __NSUBJETTINESS_HH__

#include "Njettiness.hh"

#include "fastjet/FunctionOfPseudoJet.hh"

#include <string>
#include <climits>

#ifndef G__DICTIONARY
typedef double Double32_t; // ROOT will store as 32-bit, but in code is double
#endif


FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

/// Nsubjettiness extends the concept of Njettiness to a jet shape, but other
/// than the set of particles considered, they are identical.  This class
/// wraps the core Njettiness code to provide the fastjet::FunctionOfPseudoJet
/// interface for convenience in larger analyses.  See NjettinessPlugin.hh for
/// definitions of tau_N and the constructor options.

class Nsubjettiness : public FunctionOfPseudoJet<Double32_t> {
public:

   Nsubjettiness(int N, Njettiness::AxesMode mode, double beta, double R0, double Rcutoff=std::numeric_limits<double>::max());
   
   /// returns tau_N, measured on the constituents of this jet
   Double32_t result(const PseudoJet& jet) const;

private:

   int _N;
   mutable Njettiness _njettinessFinder; // should muck with this so result can be const without this mutable

};

inline Nsubjettiness::Nsubjettiness(int N, Njettiness::AxesMode mode, double beta, double R0, double Rcutoff)
  : _N(N), _njettinessFinder(mode, NsubParameters(beta, R0, Rcutoff))
{}

inline Double32_t Nsubjettiness::result(const PseudoJet& jet) const
{
   std::vector<fastjet::PseudoJet> particles = jet.constituents();
   return _njettinessFinder.getTau(_N, particles);
}


FASTJET_END_NAMESPACE      // defined in fastjet/internal/base.hh

#endif  // __NSUBJETTINESS_HH__
