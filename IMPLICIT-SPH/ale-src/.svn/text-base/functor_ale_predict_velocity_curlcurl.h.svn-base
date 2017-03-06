#pragma once
#ifndef __FUNCTOR_ALE_PREDICT_VELOCITY_CURLCURL_H__
#define __FUNCTOR_ALE_PREDICT_VELOCITY_CURLCURL_H__

#include <iostream>
#include <string>
#include "utils.h"
#include "pair_for.h"
#include "functor.h"

namespace LAMMPS_NS {

  template<class PairIsph,
           template<typename> class GradientOperator,
           template<typename> class CurlCurl>
  class FunctorOuterAlePredictVelocityCurlCurl : public FunctorOuter<PairIsph> {
  public:
    FunctorOuterAlePredictVelocityCurlCurl(PairIsph *isph, 
                                           double dt, double gamma, 
                                           double *nu, double *rho,
                                           double **v, double **xdot, 
                                           double **f,
                                           double **vstar) 
      : FunctorOuter<PairIsph>(isph),
        _op_ale_advection(isph, v, xdot),
        _op_curlcurl_v(isph, v),
        _dt(dt),_gamma(gamma),
        _nu(nu),_rho(rho),
        _v(v),_f(f),
        _vstar(vstar) { }
    
    void enterFor();
    void operator()(const int ii);

  protected:
    FUNCTOR_REMOVE_THIS_POINTER;

    FunctorOuterAleAdvection<PairIsph,GradientOperator> _op_ale_advection;
    CurlCurl<PairIsph> _op_curlcurl_v;

    FilterBinary _filter_op;

    double _dt,_gamma,*_nu,*_rho,**_v,**_f;
    double **_vstar;
  };

  template<class PairIsph,
           template<typename> class GradientOperator,
           template<typename> class CurlCurl>
  inline void
  FunctorOuterAlePredictVelocityCurlCurl<PairIsph,GradientOperator,CurlCurl>::enterFor() {
    // --- Error check --------------------------------------------------------------------- 
    if (_filter != NULL)
      _pair->error->all(FLERR, "FunctorAlePredictVelocityCurlCurl:: Filter is not allowed");
    // --- Error check --------------------------------------------------------------------- 

    switch (_pair->boundary_particle) {
    case PairIsph::Solid: {      // sph
      _pair->error->all(FLERR, "FunctorAlePredictVelocityCurlCurl:: SPH CurlCurl does not exist");
      break;
    }
    case PairIsph::Boundary: {   // mls
      _filter_op.setPairYes(PairIsph::Fluid | PairIsph::Boundary, 
                            PairIsph::Fluid | PairIsph::Boundary);
      break;
    }
    }

    _op_ale_advection.setFilter(&_filter_op);
    _op_ale_advection.enterFor();

    _op_curlcurl_v.setFilter(&_filter_op);
    _op_curlcurl_v.enterFor();
  }

  template<class PairIsph,
           template<typename> class GradientOperator,
           template<typename> class CurlCurl>
  inline void 
  FunctorOuterAlePredictVelocityCurlCurl<PairIsph,GradientOperator,CurlCurl>::operator()(const int ii) {
    int i = _ilist[ii];
    int itype = _type[i];
    int ikind = _pair->getParticleKind(itype);

    switch (ikind) {
    case PairIsph::Boundary: 
    case PairIsph::Fluid: {
      // _vstar is already defined as BDF right hand side contribution
      {
        _op_curlcurl_v(ii);
        double *curlcurl_v = _op_curlcurl_v.getCurlCurl();

        _op_ale_advection(ii);
        double *ale_advection = _op_ale_advection.getAleAdvection();

        for (int k=0;k<_dim;++k) 
          _vstar[i][k] += _dt*(-_nu[i]*curlcurl_v[k] - ale_advection[k] + _f[i][k] + _pair->ns.g[k]);

      }

      {
        for (int k=0;k<_dim;++k) 
          _vstar[i][k] /= _gamma;
      }
      break;      
    }
    default:
      _pair->error->all(FLERR, "FunctorAlePredictVelocityCurlCurl:: Particle types are not supported");     
    }
    
  }

}

#endif

