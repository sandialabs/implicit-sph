#pragma once
#ifndef __FUNCTOR_ELECTROSTATIC_FORCE_H__
#define __FUNCTOR_ELECTROSTATIC_FORCE_H__

#include <iostream>
#include <string>
#include "utils.h"
#include "functor.h"

namespace LAMMPS_NS {

  template<class PairIsph>
  class FunctorOuterElectrostaticForce : public FunctorOuter<PairIsph> {
  public:
    FunctorOuterElectrostaticForce(PairIsph *isph, 
                                   double **phigrad,
                                   double *psi, double **psigrad,
                                   double **f)
      : FunctorOuter<PairIsph>(isph),
        _phigrad(phigrad),
        _psi(psi),
        _psigrad(psigrad),
        _f(f),
        _ezcb(isph->pb.ezcb),
        _psiref(isph->pb.psiref),
        _e(isph->pb.e),
        _gamma(isph->pb.gamma) 
    { }

    void operator()(const int ii);

  protected:
    FUNCTOR_REMOVE_THIS_POINTER;
    double **_phigrad,*_psi,**_psigrad,**_f;
    double _ezcb,_psiref,*_e,_gamma;
  };

  template<class PairIsph> inline void
  FunctorOuterElectrostaticForce<PairIsph>::operator()(const int ii) {
    const int i = _ilist[ii];

    // copy appropriate applied efield
    double e[3];
    if (_phigrad == NULL) {
      e[0] = _e[0]; e[1] = _e[1]; e[2] = _e[2];
    } else {
      const double *phigrad_at_i = _phigrad[i];
      e[0] = -phigrad_at_i[0]; e[1] = -phigrad_at_i[1]; e[2] = -phigrad_at_i[2];
    }

    // compute force
    for (int k=0;k<_dim;++k)
      _f[i][k] -= _ezcb*2.0*
        sinh(_psi[i])/(1.0 + 2.0*_gamma*pow(sinh(_psi[i]/2.0),2))
        *(-_psiref*_psigrad[i][k] + e[k]);
  }

}

#endif

