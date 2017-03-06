#pragma once
#ifndef __FUNCTOR_INCOMP_NAVIER_STOKES_BLOCK_HELMHOLTZ_H__
#define __FUNCTOR_INCOMP_NAVIER_STOKES_BLOCK_HELMHOLTZ_H__

#include <iostream>
#include <string>
#include "utils.h"
#include "pair_for.h"
#include "functor.h"
#include "modify.h"

namespace LAMMPS_NS {

  template<class PairIsph,
           template<typename> class LaplacianMatrix,
           template<typename> class Gradient,
           template<typename> class NaturalBoundary>
  class FunctorOuterIncompNavierStokesBlockHelmholtz : public FunctorOuter<PairIsph> {
  public:

    FunctorOuterIncompNavierStokesBlockHelmholtz(PairIsph *isph,
                                            double dt, double theta, double beta,
                                            double *nu, double *rho, double *p,
                                            double **f,
                                            double *b, int lda, 
                                            double *w,
                                            double** normal = NULL)
      : FunctorOuter<PairIsph>(isph),
        _op_grad_p(isph, p),
        _dt(dt),_theta(theta), _beta(beta),
        _nu(nu),_rho(rho),_f(f),
        _normal(normal),
        _lda(lda),
        _b(b),  // b is assumed to have v^{n}
        _w(w){}
    
    void enterFor();
    void operator()(const int ii);
    void exitFor();

  protected:
    Gradient<PairIsph> _op_grad_p;
    FilterBinary _filter_op_grad_p;

    double _dt,_theta, _beta, *_nu,*_rho,**_f, **_normal;
    int _lda;
    double *_b,*_w;

    FUNCTOR_REMOVE_THIS_POINTER;
  };
  
  template<class PairIsph,
           template<typename> class LaplacianMatrix,
           template<typename> class Gradient,
           template<typename> class NaturalBoundary>
  inline void
  FunctorOuterIncompNavierStokesBlockHelmholtz<PairIsph,LaplacianMatrix,Gradient,NaturalBoundary>::enterFor() {
    // --- Error check --------------------------------------------------------------------- 
    if (_filter != NULL)
      _pair->error->all(FLERR, "FunctorIncompNavierStokesBlockHelmholtz:: Filter is not allowed");



    // --- Error check --------------------------------------------------------------------- 

    for(int ib=0; ib < _dim; ++ib)
      for(int jb=0; jb < _dim; ++jb)
      _pair->A_blk(ib,jb)->PutScalar(0.0);

    for(int ib=0; ib < _dim; ++ib) {
    
    // compute Laplacian matrix
    FilterBinary filter_op;
    //  put the following code either here or below depending on whether b.c. should be prescribed according to the time scheme or fully implicitly
    //*
    if(ib == 0) {
      filter_op.setPairYes(PairIsph::Fluid, PairIsph::Solid);
      LaplacianMatrix<PairIsph> op_laplace(_pair, _dt, _nu, ib, _normal);
      op_laplace.setFilter(&filter_op);
      PairFor(op_laplace, op_laplace.getNumberOfWork());
    } else {
      NaturalBoundary<PairIsph> op_bd(_pair, -_beta*_dt, _normal, _rho, false, ib);
      PairFor(op_bd, op_bd.getNumberOfWork());
    }//*/

    filter_op.setPairYes(PairIsph::Fluid, PairIsph::Fluid);
    LaplacianMatrix<PairIsph> op_laplace(_pair, _dt, _nu, ib);
    op_laplace.setFilter(&filter_op);
    PairFor(op_laplace, op_laplace.getNumberOfWork());
    }

    for(int ib=0; ib < _dim; ++ib) {
    // b and workspace
    Epetra_Vector b(View, *_pair->nodalmap, _b+ib*_lda);
    Epetra_Vector w(View, *_pair->nodalmap, _w+ib*_lda);
    Epetra_Vector tmp(*_pair->nodalmap);
    
    // w = A u, where u is be and A is laplacian dt nu
    w.PutScalar(0);
    for(int jb=0; jb < _dim; ++jb) {
      _pair->A_blk(ib,jb)->Multiply(false, b, tmp);
      w.Update(1.0, tmp, 1.0); //w+= tmp;
    }

    w.Scale(1.0 - _theta);
    
    // complete left hand side
    for(int jb=0; jb < _dim; ++jb)
      _pair->A_blk(ib,jb)->Scale(-_theta);

    for(int ib=0; ib < _dim; ++ib) {
      // compute Laplacian matrix
      FilterBinary filter_op;
  /* fully implicitly
      if(1){//ib == 0) {
        filter_op.setPairYes(PairIsph::Fluid, PairIsph::Solid);
        LaplacianMatrix<PairIsph> op_laplace(_pair, -_dt, _nu, ib, _normal);
        op_laplace.setFilter(&filter_op);
        PairFor(op_laplace, op_laplace.getNumberOfWork());
      } else {
        NaturalBoundary<PairIsph> op_bd(_pair, _beta*_dt, _normal, _rho, false, ib);
        PairFor(op_bd, op_bd.getNumberOfWork());
      }//*/
    }

    _pair->A_blk(ib,ib)->ExtractDiagonalCopy(*_pair->A_blk.scaled_laplace_diagonal[ib]);
    }

    // initiate gradient of pressure
    _filter_op_grad_p.setPairYes(PairIsph::Fluid, PairIsph::Fluid);
    _op_grad_p.setFilter(&_filter_op_grad_p);
  }

  template<class PairIsph,
           template<typename> class LaplacianMatrix,
           template<typename> class Gradient,
           template<typename> class NaturalBoundary>
  inline void 
  FunctorOuterIncompNavierStokesBlockHelmholtz<PairIsph,LaplacianMatrix,Gradient,NaturalBoundary>::operator()(const int ii) {
    int i = _ilist[ii];
    int itype = _type[i];
    int ikind = _pair->getParticleKind(itype);

    if (_pair->ns.is_incremental_pressure_used)
      _op_grad_p(ii);

    for(int ib=0; ib < _dim; ++ib) {
      // ** modify diagonal
      double &diag = (*_pair->A_blk.diagonal[ib])[i];

      switch (ikind) {
      case PairIsph::Solid: {
        diag = 1.0;
        break;
      }
      case PairIsph::BufferDirichlet:
      case PairIsph::BufferNeumann:
      case PairIsph::Fluid: {
        // compute diagonal (I + (-theta nu dt \laplace))
        diag = 1.0 + (*_pair->A_blk.scaled_laplace_diagonal[ib])[i];

        // compute corrected gradient of p; all functor use ii index
        {
          VIEW2(_b, _lda, i, ib) += VIEW2(_w, _lda, i, ib);
          VIEW2(_b, _lda, i, ib) += _dt*(_f[i][ib]/_rho[i] + _pair->ns.g[ib]);
        }
        if (_pair->ns.is_incremental_pressure_used) {
          double *grad_p = _op_grad_p.getScalarGradient();
          VIEW2(_b, _lda, i, ib) += _dt*(-1.0/_rho[i]*grad_p[ib]);
        }
        break;
      }
      default:
        _pair->error->all(FLERR, "FunctorIncompNavierStokesBlockHelmholtz:: Particle types are not supported");
      }
    }
  }
  
  template<class PairIsph,
           template<typename> class LaplacianMatrix,
           template<typename> class Gradient,
           template<typename> class NaturalBoundary>
  inline void 
  FunctorOuterIncompNavierStokesBlockHelmholtz<PairIsph,LaplacianMatrix,Gradient,NaturalBoundary>::exitFor() {
    for(int ib=0; ib < _dim; ++ib)
      this->_pair->A_blk(ib,ib)->ReplaceDiagonalValues(*_pair->A_blk.diagonal[ib]);
  }
}

#endif

