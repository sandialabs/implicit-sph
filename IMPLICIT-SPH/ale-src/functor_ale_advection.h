#pragma once
#ifndef __FUNCTOR_ALE_ADVECTION_H__
#define __FUNCTOR_ALE_ADVECTION_H__

#include <iostream>
#include <string>
#include "utils.h"
#include "functor.h"

namespace LAMMPS_NS {
  
  template<class PairIsph,
           template<typename> class GradientOperator>
  class FunctorOuterAleAdvection : public FunctorOuter<PairIsph> {
  public:
    FunctorOuterAleAdvection(PairIsph *isph, 
                             double **v, double **xdot,  
                             double alpha = 1.0,
                             double **advection = NULL) 
      : FunctorOuter<PairIsph>(isph), 
        _op_grad(isph, alpha),
        _v(v), 
        _xdot(xdot),
        _advection(advection) { }
    
    void enterFor();
    void operator()(const int ii);
    
    double* getAleAdvection();

  protected:
    FUNCTOR_REMOVE_THIS_POINTER;
    
    GradientOperator<PairIsph> _op_grad;
    double **_v, **_xdot;
    double _advection_at_i[3];

    double **_advection;
  };

  template<class PairIsph,
           template<typename> class GradientOperator>
  inline double*
  FunctorOuterAleAdvection<PairIsph,GradientOperator>::getAleAdvection() {
    return _advection_at_i;
  }

  template<class PairIsph,
           template<typename> class GradientOperator>
  inline void
  FunctorOuterAleAdvection<PairIsph,GradientOperator>::enterFor() {
    // use the same filter 
    _op_grad.setFilter(_filter);
  }

  template<class PairIsph,
           template<typename> class GradientOperator>
  inline void 
  FunctorOuterAleAdvection<PairIsph,GradientOperator>::operator()(const int ii) {
    int i = _ilist[ii];
    int itype = _type[i];
    
    memset(_advection_at_i, 0, sizeof(double)*_dim);
    
    if (_filter != NULL && !_filter->yes(_pair->getParticleKind(itype)))
      return;
    
    // compute gradient operator
    _op_grad(ii);
    
    for (int jj=0;jj<_op_grad.getSize();++jj) {
      for (int k=0;k<_dim;++k)
        _advection_at_i[k] += (_v[i][k] - _xdot[i][k])*_op_grad.getValue(jj, k)*_v[i][k];
    }
    
    if (_advection != NULL)
      memcpy(_advection[i], _advection_at_i, sizeof(double)*_dim);
  }
  
}
#endif

