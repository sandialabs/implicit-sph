#pragma once
#ifndef __FUNCTOR_CONTINUUM_SURFACE_FORCE_H__
#define __FUNCTOR_CONTINUUM_SURFACE_FORCE_H__

#include <iostream>
#include <string>
#include "utils.h"
#include "functor.h"

namespace LAMMPS_NS {

  template<class PairIsph,
           template<typename> class Divergence>
  class FunctorOuterContinuumSurfaceForce : public FunctorOuter<PairIsph> {
  public:
    FunctorOuterContinuumSurfaceForce(PairIsph *isph, double **normal, double *mag,
                                      double **f)
      : FunctorOuter<PairIsph>(isph),
        _op_normal(isph, normal, mag),
        _normal(normal),
        _mag(mag),
        _f(f) { }

    void enterFor();
    void operator()(const int ii);

  protected:
    Divergence<PairIsph> _op_normal;
    FilterBinary _filter_op_normal;

    double **_normal, *_mag,  **_f;

    FUNCTOR_REMOVE_THIS_POINTER;
  };

  
  template<class PairIsph,
           template<typename> class Divergence>
  inline void
  FunctorOuterContinuumSurfaceForce<PairIsph,Divergence>::enterFor() {
    if (_filter == NULL)
      _pair->error->all(FLERR, "FunctorContinuumSurfaceForce:: Filter is required");

    // use the same filter
    _op_normal.setFilter(_filter);
    _op_normal.enterFor();
  }

  template<class PairIsph,
           template<typename> class Divergence>
  inline void
  FunctorOuterContinuumSurfaceForce<PairIsph,Divergence>::operator()(const int ii) {
    const int i = _ilist[ii];

    if (_mag[i] > ISPH_EPSILON) {
      _op_normal(ii);
      const double kappa = _op_normal.getDivergence();
      const double sign  = (kappa > 0.0 ? 1.0 : -1.0);
      const double alpha = _pair->st.csf.alpha*(1.0 - exp(-_pair->st.csf.kappa/(sign*kappa)));

      for (int k=0;k<_dim;++k)
        _f[i][k] -= alpha*kappa*_normal[i][k]*_mag[i];
    }
  }

}

#endif

