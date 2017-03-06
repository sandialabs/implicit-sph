#pragma once
#ifndef __FUNCTOR_INCOMP_NAVIER_STOKES_HELMHOLTZ_H__
#define __FUNCTOR_INCOMP_NAVIER_STOKES_HELMHOLTZ_H__

#include <iostream>
#include <string>
#include "utils.h"
#include "pair_for.h"
#include "functor.h"

namespace LAMMPS_NS {

  template<class PairIsph,
           template<typename> class LaplacianMatrix,
           template<typename> class Gradient>
  class FunctorOuterIncompNavierStokesHelmholtz : public FunctorOuter<PairIsph> {
  public:

    FunctorOuterIncompNavierStokesHelmholtz(PairIsph *isph, 
                                            double dt, double theta, 
                                            double *nu, double *rho, double *p,
                                            double **f,
                                            double *b, int lda, 
                                            double *w) 
      : FunctorOuter<PairIsph>(isph),
        _op_grad_p(isph, p),
        _dt(dt),_theta(theta),
        _nu(nu),_rho(rho),_f(f),
        _lda(lda),
        _b(b),  // b is assumed to have v^{n}
        _w(w) { }
    
    void enterFor();
    void operator()(const int ii);
    void exitFor();

  protected:
    Gradient<PairIsph> _op_grad_p;
    FilterBinary _filter_op_grad_p;

    double _dt,_theta,*_nu,*_rho,**_f;
    int _lda;
    double *_b,*_w;

    FUNCTOR_REMOVE_THIS_POINTER;
  };
  
  template<class PairIsph,
           template<typename> class LaplacianMatrix,
           template<typename> class Gradient>
  inline void
  FunctorOuterIncompNavierStokesHelmholtz<PairIsph,LaplacianMatrix,Gradient>::enterFor() {
    // --- Error check --------------------------------------------------------------------- 
    if (_filter != NULL)
      _pair->error->all(FLERR, "FunctorIncompNavierStokesHelmholtz:: Filter is not allowed");

    if (_pair->A.is_filled == 1) 
      _pair->error->all(FLERR, "FunctorIncompNavierStokesHelmholtz:: A is already filled");      
    // --- Error check --------------------------------------------------------------------- 

    _pair->A.crs->PutScalar(0.0);

    // compute laplacin matrix
    FilterBinary filter_op_laplace;

    // inlet boundary does not include dirichlet particles (inlet)
    // but include neumann buffer (outlet)
    filter_op_laplace.setPairYes(PairIsph::Fluid, //- PairIsph::BufferDirichlet, 
                                 PairIsph::All);

    // compute mu; use the same workspace _w
    const int nsize = (_pair->atom->nlocal + _pair->atom->nghost);
    Epetra_SerialDenseVector mu(View, _w, nsize);
    for (int i=0;i<nsize;++i)
      mu[i] = _nu[i]*_rho[i];
    
    LaplacianMatrix<PairIsph> op_laplace(_pair, _dt, mu.Values());
    op_laplace.setFilter(&filter_op_laplace);
    
    PairFor(op_laplace, op_laplace.getNumberOfWork());

    // apply inverse of rho to A
    Epetra_Vector tmp(View, *_pair->nodalmap, _w);
    tmp.Reciprocal(Epetra_Vector(View, *_pair->nodalmap, _rho));

    _pair->A.crs->LeftScale(tmp);

    // b and workspace
    Epetra_MultiVector b(View, *_pair->nodalmap, _b, _lda, _dim);
    Epetra_MultiVector w(View, *_pair->nodalmap, _w, _lda, _dim);
    
    // w = A u, where u is be and A is laplacian dt nu
    _pair->A.crs->Multiply(false, b, w);
    w.Scale(1.0 - _theta);
    
    // complete left hand side
    _pair->A.crs->Scale(-_theta);
    _pair->A.crs->ExtractDiagonalCopy(*_pair->A.scaled_laplace_diagonal);   

    // initiate gradient of pressure
    _filter_op_grad_p.setPairYes(PairIsph::Fluid, // - PairIsph::BufferDirichlet, 
                                 PairIsph::Fluid);
    _op_grad_p.setFilter(&_filter_op_grad_p);
  }

  template<class PairIsph,
           template<typename> class LaplacianMatrix,
           template<typename> class Gradient>
  inline void 
  FunctorOuterIncompNavierStokesHelmholtz<PairIsph,LaplacianMatrix,Gradient>::operator()(const int ii) {
    int i = _ilist[ii];
    int itype = _type[i];
    int ikind = _pair->getParticleKind(itype);

    // ** modify diagonal
    double &diag = (*_pair->A.diagonal)[i];

    switch (ikind) {
    // case PairIsph::BufferDummy:
    // case PairIsph::BufferDirichlet:
    case PairIsph::Solid: {
      diag = 1.0; 
      break;
    }
      //case PairIsph::BufferDummy:
    case PairIsph::BufferDirichlet:
    case PairIsph::BufferNeumann:
    case PairIsph::Fluid: {
      // compute diagonal (I + (-theta nu dt \laplace))
      diag = 1.0 + (*_pair->A.scaled_laplace_diagonal)[i];
      
      // compute corrected gradient of p; all functor use ii index
      {
        for (int k=0;k<_dim;++k) {
          VIEW2(_b, _lda, i, k) += VIEW2(_w, _lda, i, k);
          VIEW2(_b, _lda, i, k) += _dt*(_f[i][k]/_rho[i] + _pair->ns.g[k]);
        }

        if (_pair->ns.is_incremental_pressure_used) {
          _op_grad_p(ii);
          double *grad_p = _op_grad_p.getScalarGradient();
          for (int k=0;k<_dim;++k) 
            VIEW2(_b, _lda, i, k) += _dt*(-1.0/_rho[i]*grad_p[k]);
        }
      }
      break;      
    }
    default:
      _pair->error->all(FLERR, "FunctorIncompNavierStokesHelmholtz:: Particle types are not supported");     
    }
  }
  
  template<class PairIsph,
           template<typename> class LaplacianMatrix,
           template<typename> class Gradient>
  inline void 
  FunctorOuterIncompNavierStokesHelmholtz<PairIsph,LaplacianMatrix,Gradient>::exitFor() {
    _pair->A.crs->ReplaceDiagonalValues(*_pair->A.diagonal);
  }
}

#endif

