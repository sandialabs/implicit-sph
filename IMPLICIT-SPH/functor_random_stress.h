#pragma once
#ifndef __FUNCTOR_RANDOM_STRESS_H__
#define __FUNCTOR_RANDOM_STRESS_H__

#include "math.h"
#include <iostream>
#include <string>
#include "utils.h"
#include "functor.h"

namespace LAMMPS_NS {

  template<class PairIsph,
           template<typename> class Divergence>
  class FunctorOuterRandomStress : public FunctorOuter<PairIsph> {
  public:
    FunctorOuterRandomStress(PairIsph *isph, 
                             double dt, double *nu, double *rho,
                             double **f,
                             double **rstress_x,
                             double **rstress_y,
                             double **rstress_z) 
      : FunctorOuter<PairIsph>(isph),
        _dt(dt), _nu(nu),_rho(rho), _f(f),
        _op_force_x(isph, rstress_x, -1.0, NULL),
        _op_force_y(isph, rstress_y, -1.0, NULL),
        _op_force_z(isph, rstress_z, -1.0, NULL) { }

    void enterFor();
    void operator()(const int ii);

  protected:
    Divergence<PairIsph> _op_force_x, _op_force_y, _op_force_z;
    double _dt, *_nu, *_rho, **_f, _kBT;

    FUNCTOR_REMOVE_THIS_POINTER;
  };

  template<class PairIsph,
           template<typename> class Divergence>
  inline void 
  FunctorOuterRandomStress<PairIsph,Divergence>::enterFor() {
    // use the same filter inherited from this functor
    _op_force_x.setFilter(_filter);
    _op_force_y.setFilter(_filter);
    _op_force_z.setFilter(_filter);

    _kBT = _pair->ns.kBT;
  }

  template<class PairIsph,
           template<typename> class Divergence>
  inline void 
  FunctorOuterRandomStress<PairIsph,Divergence>::operator()(const int ii) {
    const int i = _ilist[ii];
    const int itype = _type[i];
    const int ikind = _pair->getParticleKind(itype);

    // here we execute divergence for each particle and store its contribution to force
    // this will filter other particles than fluid
    if (_filter != NULL && _filter->yes(ikind)) {
      _op_force_x(ii);
      _op_force_y(ii);
      
      const double sq_variance = sqrt(2.0*_kBT*_nu[i]*_rho[i]/_dt/_vfrac[i]);
      _f[i][0] += _op_force_x.getDivergence()*sq_variance;
      _f[i][1] += _op_force_y.getDivergence()*sq_variance;
      
      if (_dim > 2) {
        _op_force_z(ii);
        _f[i][2] += _op_force_z.getDivergence()*sq_variance;      
      }
    }
  }
  
}

#endif

