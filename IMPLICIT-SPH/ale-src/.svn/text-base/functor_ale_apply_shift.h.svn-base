#pragma once
#ifndef __FUNCTOR_ALE_APPLY_SHIFT_H__
#define __FUNCTOR_ALE_APPLY_SHIFT_H__

#include <iostream>
#include <string>
#include "utils.h"
#include "functor.h"

namespace LAMMPS_NS {

  template<class PairIsph>
  class FunctorOuterAleApplyShift : public FunctorOuter<PairIsph> {
  public:
    FunctorOuterAleApplyShift(PairIsph *isph, 
                              double dt, double gamma, 
                              double **dr,
                              double **xdot)
      : FunctorOuter<PairIsph>(isph),
        _dt(dt), _gamma(gamma),
        _dr(dr),
        _xdot(xdot) { }

    void enterFor();
    void operator()(const int ii);

  protected:
    FUNCTOR_REMOVE_THIS_POINTER;

    double _dt,_gamma,**_dr,**_xdot;

  };

  template<class PairIsph>
  inline void
  FunctorOuterAleApplyShift<PairIsph>::enterFor() {
    if (_filter != NULL)
      _pair->error->all(FLERR, "FunctorAleApplyShift:: Filter is not allowed");
  }

  template<class PairIsph>
  inline void
  FunctorOuterAleApplyShift<PairIsph>::operator()(const int ii) {
    int i = _ilist[ii]; 
    int itype = _type[i]; 

    // shifting is allowed for fluid particles only
    if (_pair->isParticleFixed(itype))
      return;
    
    // after time is advanced, xdot include new extrapolated velocity field
    
    // update xdot
    for (int k=0;k<_dim;++k) 
      _xdot[i][k] += _gamma/_dt*_dr[i][k];

    // update position
    for (int k=0;k<_dim;++k) 
      _x[i][k] += _dr[i][k];
  }

}

#endif

