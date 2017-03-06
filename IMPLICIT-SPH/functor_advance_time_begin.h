#pragma once
#ifndef __FUNCTOR_ADVANCE_TIME_BEGIN_H__
#define __FUNCTOR_ADVANCE_TIME_BEGIN_H__

#include <iostream>
#include <string>
#include "utils.h"
#include "functor.h"

namespace LAMMPS_NS {

  // compute dx and update pressure
  template<class PairIsph,
           template<typename> class Gradient>
  class FunctorOuterAdvanceTimeBegin : public FunctorOuter<PairIsph> {
  public:
    FunctorOuterAdvanceTimeBegin(PairIsph *isph, double dt, 
                                 double **v, double *p) 
      : FunctorOuter<PairIsph>(isph),
        _op_p(isph, p), 
        _dt(dt),
        _v(v),_vnp1(isph->vstar),
        _p(p),_dp(isph->dp) { }
    void enterFor();
    void operator()(const int ii);
    void exitFor();

  protected:
    Gradient<PairIsph> _op_p; 
    double _dt,**_v,**_vnp1,*_p,*_dp;

    FilterBinary _filter_op_p;

    FUNCTOR_REMOVE_THIS_POINTER;
  };

  template<class PairIsph,
           template<typename> class Gradient>
  inline void
  FunctorOuterAdvanceTimeBegin<PairIsph,Gradient>::enterFor() {
    if (_filter != NULL)
      _pair->error->all(FLERR, "FunctorAdvanceTimeBegin:: Filter is not allowed");

    _filter_op_p.setPairYes(PairIsph::Fluid, 
                            PairIsph::Fluid);
    _op_p.setFilter(&_filter_op_p);
  }   

  template<class PairIsph,
           template<typename> class Gradient>
  inline void
  FunctorOuterAdvanceTimeBegin<PairIsph,Gradient>::operator()(const int ii) {
    int i = _ilist[ii]; 
    int itype = _type[i];

    // store dx
    double dx[3] = {};
    for (int k=0;k<_dim;++k)
      dx[k] = 0.5*_dt*(_vnp1[i][k] + _v[i][k]);
    
    _op_p(ii);
    double *grad_p = _op_p.getScalarGradient();
    
    if (_filter_op_p.yes(_pair->getParticleKind(itype))) 
      _dp[i] = util.dotVectors(_dim, grad_p, dx);
    else 
      _dp[i] = 0.0;
  }

  template<class PairIsph,
           template<typename> class Gradient> 
  inline void
  FunctorOuterAdvanceTimeBegin<PairIsph,Gradient>::exitFor() {
    _pair->comm_variable = PairIsph::DeltaP;
    _pair->comm_forward = 1;
    _pair->comm->forward_comm_pair(_pair);
  }

}

#endif

