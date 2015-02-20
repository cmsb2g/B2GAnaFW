import math
import ROOT

def solve_nu_tmass(bjet, lepton, vnu, tmass=175.0) : 

    #
    # Purpose: Solve for the neutrino longitudinal z-momentum that makes
    #          the leptonic top have mass TMASS.
    # Copied from TQAF TopHitFit routine from
    # TopQuarkAnalysis/TopHitFit/src/Top_Decaykin.cc
    #
    # Inputs:
    #   tmass -       The desired top mass.
    #
    # Outputs:
    #   nuz1 -        First solution (smaller absolute value).
    #   nuz2 -        Second solution.
    #
    # Returns:
    #   True if there was a real solution.  False if there were only
    #   imaginary solutions.  (In that case, we just set the imaginary
    #   part to zero.)
    #

    discrim_flag = True

    nuz1 = 0.0
    nuz2 = 0.0

    cprime = lepton + bjet
    alpha1 = tmass*tmass - cprime.M2()
    a = 2. * 4. * (cprime.Pz()*cprime.Pz() - cprime.E()*cprime.E())
    alpha = alpha1 + 2.*(cprime.Px()*vnu.Px() + cprime.Py()*vnu.Py())
    b = 4. * alpha * cprime.Pz()
    c = alpha*alpha - 4. * cprime.E()*cprime.E() * vnu.Vect().Perp2()
    d = b*b - 2.*a*c
    if d < 0 :
        discrim_flag = False
        d = 0


    dd = math.sqrt (d)
    nuz1 = (-b + dd)/a
    nuz2 = (-b - dd)/a
    if abs (nuz1) > abs (nuz2) :
        nuz1,nuz2 = nuz2,nuz1

    return discrim_flag, nuz1, nuz2



def solve_nu(vlep, vnu, wmass=80.4) : 
    #
    # Purpose: Solve for the neutrino longitudinal z-momentum that makes
    #          the leptonic W have mass WMASS.
    # Copied from TQAF TopHitFit routine from
    # TopQuarkAnalysis/TopHitFit/src/Top_Decaykin.cc
    #
    # Inputs:
    #   vlep -        The lepton p4
    #   vnu -         The neutrino p4
    #   wmass -       The desired W mass.
    #
    # Outputs:
    #   nuz1 -        First solution (smaller absolute value).
    #   nuz2 -        Second solution.
    #
    # Returns:
    #   True if there was a real solution.  False if there were only
    #   imaginary solutions.  (In that case, we just set the imaginary
    #   part to zero.)
    #
    
    discrim_flag = True

    nuz1 = 0.0
    nuz2 = 0.0

    x = vlep.Px()*vnu.Px() + vlep.Py()*vnu.Py() + wmass*wmass/2.
    a = vlep.Pz()*vlep.Pz() - vlep.E()*vlep.E()
    b = 2.*x*vlep.Pz()
    c = x*x - vnu.Vect().Perp2() * vlep.E()*vlep.E()

    d = b*b - 4.*a*c
    if d < 0 :
        d = 0
        discrim_flag = False


    nuz1 = (-b + math.sqrt (d))/2./a
    nuz2 = (-b - math.sqrt (d))/2./a
    if abs (nuz1) > abs (nuz2) :
        nuz1, nuz2 = nuz2, nuz1

    return discrim_flag, nuz1, nuz2
