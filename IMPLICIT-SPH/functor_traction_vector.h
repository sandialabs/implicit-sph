#pragma once
#ifndef __FUNCTOR_TRACTION_VECTOR_H__
#define __FUNCTOR_TRACTION_VECTOR_H__

#include <iostream>
#include <string>
#include "utils.h"
#include "functor.h"

namespace LAMMPS_NS {

  template<class PairIsph,
           template<typename> class Gradient>
  class FunctorOuterTractionVector : public FunctorOuter<PairIsph> {
  public:
    FunctorOuterTractionVector(PairIsph *isph, 
                             double *rho, double *nu,
                             double *p, double **v,
                             double **normal = NULL, 
                             double **traction = NULL)
      : FunctorOuter<PairIsph>(isph),
        _op_v(isph, v),
        _rho(rho),
        _nu(nu),
        _p(p),
        _normal(normal),
        _traction(traction)
    { }

    void enterFor();
    void operator()(const int ii);

    double *getStressTensor(const int k = 0) const;
    double *getTractionVector() const;

  protected:
    FUNCTOR_REMOVE_THIS_POINTER;

    Gradient<PairIsph> _op_v;

    double *_rho, *_nu, *_p, **_normal, **_traction;
    double _stress_at_i[3][3], _traction_at_i[3];
  };

  template<class PairIsph,
           template<typename> class Gradient> 
  inline double*
  FunctorOuterTractionVector<PairIsph,Gradient>::getStressTensor(const int k) const {
    return _stress_at_i[k];
  }

  template<class PairIsph,
           template<typename> class Gradient> 
  inline double*
  FunctorOuterTractionVector<PairIsph,Gradient>::getTractionVector() const {
    return _traction_at_i;
  }

  template<class PairIsph,
           template<typename> class Gradient> 
  inline void
  FunctorOuterTractionVector<PairIsph,Gradient>::enterFor() {
    // use the same filter 
    _op_v.setFilter(_filter);
  }

  template<class PairIsph,
           template<typename> class Gradient> 
  inline void
  FunctorOuterTractionVector<PairIsph,Gradient>::operator()(const int ii) {
    const int i = _ilist[ii];
    const int itype = _type[i];
    const int ikind = _pair->getParticleKind(itype);

    memset(&_stress_at_i[0][0], 0, sizeof(double)*9);
    memset(&_traction_at_i[0],  0, sizeof(double)*3);

    if (_filter != NULL && !_filter->yes(ikind)) {
      // do nothing
    } else {

      // stress = -p I + T;
      // T = 2 mu tau
      // tau = 1/2 * (\grad(v) + \grad(v)')

      const double mu = _nu[i]*_rho[i];
      
      _op_v(ii);
      double *stress = &_stress_at_i[0][0];
      const double *gradv = _op_v.getVectorGradient();
      
      for (int k2=0;k2<_dim;++k2)
        for (int k1=0;k1<_dim;++k1)
          VIEW2(stress, 3, k1, k2) = (VIEW2(gradv, 3, k1, k2) +
                                      VIEW2(gradv, 3, k2, k1))*mu;

      for (int k=0;k<_dim;++k)
        VIEW2(stress, 3, k, k) -= _p[i];

      if (_normal != NULL) 
        for (int k=0;k<_dim;++k)
          _traction_at_i[k] = util.dotVectors(3, _stress_at_i[k], _normal[i]); 
      
      if (_traction != NULL)
        memcpy(_traction[i], _traction_at_i, sizeof(double)*_dim);
    }

  }

}

#endif

