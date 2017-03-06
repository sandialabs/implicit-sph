#pragma once
#ifndef __FUNCTOR_ALE_CORRECT_VELOCITY_H__
#define __FUNCTOR_ALE_CORRECT_VELOCITY_H__

#include <iostream>
#include <string>
#include "utils.h"
#include "functor.h"

namespace LAMMPS_NS {
  
  template<class PairIsph,
           template<typename> class Gradient>
  class FunctorOuterAleCorrectVelocity : public FunctorOuter<PairIsph> {
  public:
    FunctorOuterAleCorrectVelocity(PairIsph *isph, 
                                   double dt, double gamma, 
                                   double *rho,
                                   double *p, 
                                   double **vstar) 
      : FunctorOuter<PairIsph>(isph),
        _op_grad_p(isph, p, dt/gamma),
        _rho(rho),
        _vstar(vstar) { }

    void enterFor();
    void operator()(const int ii);
    void exitFor();

  protected:
    FUNCTOR_REMOVE_THIS_POINTER;

    Gradient<PairIsph> _op_grad_p;
    FilterBinary _filter_op_grad_p;

    double *_rho,**_vstar;
    int _extend;
  };

  template<class PairIsph,
           template<typename> class Gradient> 
  inline void
  FunctorOuterAleCorrectVelocity<PairIsph,Gradient>::enterFor() {
    if (_filter != NULL)
      _pair->error->all(FLERR, "FunctorAleCorrectVelocity:: Filter is not allowed - internal filter is consistent to its member operators");

    switch (_pair->boundary_particle) {
    case PairIsph::Solid: {
      _filter_op_grad_p.setPairYes(PairIsph::Fluid, PairIsph::Fluid);
      break;
    }
    case PairIsph::Boundary: {
      _filter_op_grad_p.setPairYes(PairIsph::Fluid | PairIsph::Boundary, 
                                   PairIsph::Fluid | PairIsph::Boundary);
      break;
    }
    }
    _op_grad_p.setFilter(&_filter_op_grad_p);  
  }

  template<class PairIsph,
           template<typename> class Gradient> 
  inline void
  FunctorOuterAleCorrectVelocity<PairIsph,Gradient>::operator()(const int ii) {
    int i = _ilist[ii]; 

    // compute graient p
    _op_grad_p(ii);
    double *grad_p = _op_grad_p.getScalarGradient();

    for (int k=0;k<_dim;++k)
      _vstar[i][k] -= grad_p[k]/_rho[i];
  }

  template<class PairIsph,
           template<typename> class Gradient>
  inline void
  FunctorOuterAleCorrectVelocity<PairIsph,Gradient>::exitFor() {
    _pair->comm_variable = PairIsph::Vstar;
    _pair->comm_forward = 3;
    _pair->comm->forward_comm_pair(_pair);
  }

}

#endif

