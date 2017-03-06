#pragma once
#ifndef __FUNCTOR_ALE_INCOMP_NAVIER_STOKES_HELMHOLTZ_H__
#define __FUNCTOR_ALE_INCOMP_NAVIER_STOKES_HELMHOLTZ_H__

#include <iostream>
#include <string>
#include "utils.h"
#include "pair_for.h"
#include "functor.h"

namespace LAMMPS_NS {

  template<class PairIsph,
           template<typename> class LaplacianMatrix,
           template<typename> class GradientOperator,
           template<typename> class Gradient,
           template<typename> class GradientInner>
  class FunctorOuterAleIncompNavierStokesHelmholtz : public FunctorOuter<PairIsph> {
  public:
    FunctorOuterAleIncompNavierStokesHelmholtz(PairIsph *isph, 
                                               double dt, double gamma,
                                               double *nu, double *rho, 
                                               double **vstar,             // improved estimation
                                               double **v, double **xdot,  // extrapolation
                                               double *b, int lda)         
      : FunctorOuter<PairIsph>(isph),
        _op_ale_advection(isph, v, xdot),
        _op_curlcurl_v(isph, v), 
        _dt(dt),_gamma(gamma),
        _nu(nu),_rho(rho),
        _vstar(vstar),
        _v(v),_xdot(xdot),
        _b(b),_lda(lda) { }
    
    void enterFor();
    void operator()(const int ii);
    void exitFor();

  protected:
    FUNCTOR_REMOVE_THIS_POINTER;

    FunctorOuterAleAdvection<PairIsph,GradientOperator> _op_ale_advection;
    FunctorOuterCurlCurl<PairIsph,Gradient,GradientInner> _op_curlcurl_v;
    FilterBinary _filter_op;

    double _dt,_gamma,*_nu,*_rho,**_vstar,**_v,**_xdot,*_b;
    int _lda;
  };
  
  template<class PairIsph,
           template<typename> class LaplacianMatrix,
           template<typename> class GradientOperator,
           template<typename> class Gradient,
           template<typename> class GradientInner>
  inline void
  FunctorOuterAleIncompNavierStokesHelmholtz<PairIsph,
                                             LaplacianMatrix,
                                             GradientOperator,
                                             Gradient,
                                             GradientInner>::enterFor() {
    // --- Error check --------------------------------------------------------------------- 
    if (_filter != NULL)
      _pair->error->all(FLERR, "FunctorAleIncompNavierStokesHelmholtz:: Filter is not allowed");

    if (_pair->A.is_filled == 1) 
      _pair->error->all(FLERR, "FunctorIncompNavierStokesHelmholtz:: A is already filled");      
    // --- Error check --------------------------------------------------------------------- 
    
    // compute laplacin/advection matrix
    switch (_pair->boundary_particle) {
    case PairIsph::Solid: {      // sph
      _filter_op.setPairYes(PairIsph::Fluid,
                            PairIsph::Fluid | PairIsph::Solid);
      break;
    }
    case PairIsph::Boundary: {   // mls
      _pair->error->all(FLERR, "FunctorIncompNavierStokesHelmholtz:: Boundary particles are not allowed");
      break;
    }
    }

    _pair->A.crs->PutScalar(0.0);
    
    LaplacianMatrix<PairIsph> op_laplace(_pair, -_dt, _nu);
    op_laplace.setFilter(&_filter_op);
    
    PairFor(op_laplace, op_laplace.getNumberOfWork());
    
    FunctorOuterAleAdvectionMatrix<PairIsph,GradientOperator> op_ale_advection(_pair, _vstar, _xdot, _dt);
    op_ale_advection.setFilter(&_filter_op);
    
    PairFor
      (op_ale_advection, op_ale_advection.getNumberOfWork());
    
    // complete left hand side
    _pair->A.crs->ExtractDiagonalCopy(*_pair->A.scaled_laplace_diagonal);   

    _op_ale_advection.setFilter(&_filter_op);
    _op_ale_advection.enterFor();

    _op_curlcurl_v.setFilter(&_filter_op);
    _op_curlcurl_v.enterFor();
  }

  template<class PairIsph,
           template<typename> class LaplacianMatrix,
           template<typename> class GradientOperator,
           template<typename> class Gradient,
           template<typename> class GradientInner>
  inline void 
  FunctorOuterAleIncompNavierStokesHelmholtz<PairIsph,
                                             LaplacianMatrix,
                                             GradientOperator,
                                             Gradient,
                                             GradientInner>::operator()(const int ii) {
    int i = _ilist[ii];
    int itype = _type[i];
    int ikind = _pair->getParticleKind(itype);

    // ** modify diagonal 
    double &diag = (*_pair->A.diagonal)[i];

    switch (ikind) {
    case PairIsph::Solid: {
      diag = 1.0;
      
      // solid particle take extrapolated velocity (here const)
      for (int k=0;k<_dim;++k)
        VIEW2(_b, _lda, i, k) = _v[i][k];
      break;
    }
    case PairIsph::Fluid: {
      // compute diagonal 
      diag = _gamma + (*_pair->A.scaled_laplace_diagonal)[i];

      // compute right hand side: advection + curlcurl + vstar*gamma
      _op_ale_advection(ii);
      double *ale_advection_at_i = _op_ale_advection.getAleAdvection();

      _op_curlcurl_v(ii);
      double *curlcurl_v_at_i = _op_curlcurl_v.getCurlCurl();
      
      for (int k=0;k<_dim;++k)
        VIEW2(_b, _lda, i, k) = _gamma*_vstar[i][k] + _dt*(ale_advection_at_i[k] + _nu[i]*curlcurl_v_at_i[k]);

      break;
    }
    default:
      _pair->error->all(FLERR, "FunctorAleIncompNavierStokesHelmholtz:: Particle types are not supported");     
    }
  }
  
  template<class PairIsph,
           template<typename> class LaplacianMatrix,
           template<typename> class GradientOperator,
           template<typename> class Gradient,
           template<typename> class GradientInner>
  inline void
  FunctorOuterAleIncompNavierStokesHelmholtz<PairIsph,
                                             LaplacianMatrix,
                                             GradientOperator,
                                             Gradient,
                                             GradientInner>::exitFor() {
    _pair->A.crs->ReplaceDiagonalValues(*_pair->A.diagonal);
  }  

}

#endif

