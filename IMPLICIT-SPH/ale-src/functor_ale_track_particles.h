#pragma once
#ifndef __FUNCTOR_ALE_TRACK_PARTICLES_H__
#define __FUNCTOR_ALE_TRACK_PARTICLES_H__

#include <iostream>
#include <string>
#include "utils.h"
#include "functor.h"

namespace LAMMPS_NS {

  template<class PairIsph>
  class FunctorOuterAleTrackParticles : public FunctorOuter<PairIsph> {
  public:
    FunctorOuterAleTrackParticles(PairIsph *isph, 
                                  double dt, double gamma,
                                  double **xdot) 
      : FunctorOuter<PairIsph>(isph),
        _dt(dt),_gamma(gamma),
        _xdot(xdot) { }
    void enterFor();
    void operator()(const int ii);

  protected:
    FUNCTOR_REMOVE_THIS_POINTER;
    double _dt,_gamma,**_xdot;

  };
  
  template<class PairIsph> 
  inline void
  FunctorOuterAleTrackParticles<PairIsph>::enterFor() {
    if (_filter != NULL)
      _pair->error->all(FLERR, "FunctorAleTrackParticles:: Filter is not allowed");
  }

  template<class PairIsph> 
  inline void
  FunctorOuterAleTrackParticles<PairIsph>::operator()(const int ii) {
    int i = _ilist[ii]; 
    int itype = _type[i];


    if (!_pair->isParticleFixed(itype)) 
      for (int k=0;k<_dim;++k) 
        _x[i][k] += _dt*_xdot[i][k];

    // _x should include BDF right hand side
    for (int k=0;k<_dim;++k)     
      _x[i][k] /= _gamma;
  }

}

#endif

