#pragma once
#ifndef __FUNCTOR_CORRECT_VELOCITY_H__
#define __FUNCTOR_CORRECT_VELOCITY_H__


#include <iostream>
#include <string>
#include "utils.h"
#include "functor.h"
#include "functor_gradient.h"
//#include "functor_mls_gradient.h"

namespace LAMMPS_NS {

  template<class PairIsph,
           template<typename> class Gradient>
  class FunctorOuterCorrectVelocity : public FunctorOuter<PairIsph> {
  public:
    FunctorOuterCorrectVelocity(PairIsph *isph, 
                                double dt, double *rho, 
                                double *dp, double **vstar)
      : FunctorOuter<PairIsph>(isph),
        _op_dp(isph, dp, 1.0, NULL),
        _dt(dt),_rho(rho),
        _vstar(vstar) { }
    void enterFor();
    void operator()(const int ii);
    void exitFor();

  protected:

    Gradient<PairIsph> _op_dp;
    FilterBinary _filter_op_dp;

    double _dt, *_rho, **_vstar;

    FUNCTOR_REMOVE_THIS_POINTER;
  };

  template<class PairIsph,
           template<typename> class Gradient>
  inline void 
  FunctorOuterCorrectVelocity<PairIsph,Gradient>::enterFor() {
    if (_filter != NULL)
      _pair->error->all(FLERR, "FunctorCorrectVelocity:: Filter is not allowed - internal filter is consistent to its member operators");    

    // gradient operator only consider fluid particles
    _filter_op_dp.setPairYes(PairIsph::Fluid, PairIsph::Fluid);
    _op_dp.setFilter(&_filter_op_dp);
  }

  template<class PairIsph,
           template<typename> class Gradient>
  inline void 
  FunctorOuterCorrectVelocity<PairIsph,Gradient>::operator()(const int ii) {
    int i = _ilist[ii];
    int itype = _type[i];

    if (!_filter_op_dp.yes(_pair->getParticleKind(itype)))
      return;

    // minus sign is accounted in evaluating laplacian
    _op_dp(ii);
    double *grad_dp = _op_dp.getScalarGradient();
    
    // velocity correction
    for (int k=0;k<_dim;++k)
      _vstar[i][k] -= _dt/_rho[i]*grad_dp[k];
  }

  template<class PairIsph,
           template<typename> class Gradient>
  inline void 
  FunctorOuterCorrectVelocity<PairIsph,Gradient>::exitFor() {
    _pair->comm_variable = PairIsph::Vstar;
    _pair->comm_forward = 3;
    _pair->comm->forward_comm_pair(_pair);
  }

}

#endif

