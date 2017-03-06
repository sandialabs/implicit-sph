#pragma once
#ifndef __FUNCTOR_APPLY_SHIFT_H__
#define __FUNCTOR_APPLY_SHIFT_H__

#include <iostream>
#include <string>
#include "utils.h"
#include "functor.h"

namespace LAMMPS_NS {

  template<class PairIsph,
           template<typename> class Gradient>
  class FunctorOuterApplyShift : public FunctorOuter<PairIsph> {
  public:
    FunctorOuterApplyShift(PairIsph *isph, double **dr, double **v, double *p)
      : FunctorOuter<PairIsph>(isph),
        _op_pv(isph, p, v), 
        _dr(dr),_v(v),_p(p),
        _c(NULL) { 
      for (int cid=0;cid<ISPH_MAX_CONCENTRATION;++cid)
        _op_c[cid] = NULL;
    }
    void enterFor();
    void operator()(const int ii);
    void exitFor();

  protected:
    double **_dr;
    FilterBinary _filter_op;

    double *_p,**_v;
    Gradient<PairIsph> _op_pv; 

    // optional 
    double **_c;
    Gradient<PairIsph> *_op_c[ISPH_MAX_CONCENTRATION];

    FUNCTOR_REMOVE_THIS_POINTER;
  };

  template<class PairIsph,
           template<typename> class Gradient>
  inline void
  FunctorOuterApplyShift<PairIsph,Gradient>::enterFor() {
    if (_filter != NULL)
      _pair->error->all(FLERR, "FunctorApplyShift:: Filter is not allowed - internal filter only considers fluid/fluid interaction");

    // shifting is allowed for fluid particles only
    // solid particles have zero pressure and may have non-zero velocity field 
    // to apply dirichlet-like conditions
    // here we compute both gradient for pressure and velocity
    // if they are computed separate, correct filter should be
    // pressure : F + F, velocity : F + A
    // as pressure on solid particles are zero, it does not matter to
    // apply F + A for both
    _filter_op.setPairYes(PairIsph::Fluid, PairIsph::All);

    // pressure and velocity shift
    _op_pv.setFilter(&_filter_op);

    // optional concentration shift
    if (_pair->tr.is_enabled) {
      _c = _pair->atom->concentration;
      for (int cid=0;cid<ISPH_MAX_CONCENTRATION;++cid) 
        if (_pair->tr.mask[cid]) {
          _op_c[cid] = new Gradient<PairIsph>(_pair, _c[cid]);
          _op_c[cid]->setFilter(&_filter_op);
        } 
    }
  }
  
  template<class PairIsph,
           template<typename> class Gradient>
  inline void
  FunctorOuterApplyShift<PairIsph,Gradient>::operator()(const int ii) {
    int i = _ilist[ii]; 
    int itype = _type[i];

    // shifting is allowed for fluid particles only
    if (_pair->isParticleFixed(itype))
      return;

    // compute gradients
    _op_pv(ii);

    // update pressure
    double *grad_p =  _op_pv.getScalarGradient();
    _p[i] += util.dotVectors(_dim, grad_p, _dr[i]);

    // update velocity
    for (int k=0;k<_dim;++k) {
      double *grad_v = _op_pv.getVectorGradient(k);
      _v[i][k] += util.dotVectors(_dim, grad_v, _dr[i]);
    }

    // update concentration
    for (int cid=0;cid<ISPH_MAX_CONCENTRATION;++cid) {
      if (_op_c[cid] != NULL) {
        (*_op_c[cid])(ii);
        double *grad_c = _op_c[cid]->getScalarGradient();
        _c[cid][i] += util.dotVectors(_dim, grad_c, _dr[i]);
      }
    }
    
    // update position
    for (int k=0;k<_dim;++k) 
      _x[i][k] += _dr[i][k];
  }

  template<class PairIsph,
           template<typename> class Gradient>
  inline void
  FunctorOuterApplyShift<PairIsph,Gradient>::exitFor() {
    if (_pair->tr.is_enabled) 
      for (int cid=0;cid<ISPH_MAX_CONCENTRATION;++cid) 
        if (_op_c[cid] != NULL)
          delete _op_c[cid];
  }
}

#endif

