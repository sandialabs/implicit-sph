#pragma once
#ifndef __FUNCTOR_INCOMP_NAVIER_STOKES_POISSON_H__
#define __FUNCTOR_INCOMP_NAVIER_STOKES_POISSON_H__

#include "Epetra_CrsMatrix.h"

#include <iostream>
#include <string>
#include "utils.h"
#include "pair_for.h"
#include "functor.h"

namespace LAMMPS_NS {

  template<class PairIsph, 
           template<typename> class LaplacianMatrix, 
           template<typename> class Divergence,
           template<typename> class GradientOperator>
  class FunctorOuterIncompNavierStokesPoisson : public FunctorOuter<PairIsph> {
  public:
    
    FunctorOuterIncompNavierStokesPoisson(PairIsph *isph, 
                                          double dt, double **normal, 
                                          double *rho, double **vstar, double *p,
                                          double *b,
                                          double *w) 
      : FunctorOuter<PairIsph>(isph),
        _op_vstar(isph, vstar, 1.0, NULL),
        _dt(dt),_normal(normal),_rho(rho),_p(p),
        _b(b),_w(w) { }
    
    void enterFor();
    void operator()(const int ii);
    void exitFor();
    
  protected:
    Divergence<PairIsph> _op_vstar;            
    FilterBinary _filter_op_vstar;
    
    bool _is_once,_apply_homogeneous_neumann_boundary;
    double _dt,**_normal,*_rho,*_p,*_b,*_w;
    
    FUNCTOR_REMOVE_THIS_POINTER;
  };
  

  template<class PairIsph,
           template<typename> class LaplacianMatrix,
           template<typename> class Divergence,
           template<typename> class GradientOperator>
  inline void
  FunctorOuterIncompNavierStokesPoisson<PairIsph,LaplacianMatrix,Divergence,GradientOperator>::enterFor() {
    // --- Error check --------------------------------------------------------------------- 
    if (_filter != NULL)
      _pair->error->all(FLERR, "FunctorNavierIncompNavierStokesPoisson:: Filter is not allowed");

    if (_pair->A.is_filled == 1) 
      _pair->error->all(FLERR, "FunctorIncompNavierStokesPoisson:: A is already filled"); 
    // --- Error check --------------------------------------------------------------------- 

    _pair->A.crs->PutScalar(0.0);
    
    // compute laplace matrix
    // in principle, this poisson operator should be restricted to F + F
    // as pressure on the solid particles is zero, including solid particles
    // does not matter and it resolve the singular poisson problem.
    // also an appropriate boundary condition is homogeneous neumann boundary
    // hence, this tweaking is okay.
    {
      FilterBinary filter_op_laplace;

      switch (_pair->ns.singular_poisson) {
      case PairIsph::NotSingularPoisson:
        // set pressure on solid zero; solid particle has unit diagonals
        filter_op_laplace.setPairYes(PairIsph::Fluid, //- PairIsph::BufferNeumann, 
                                     PairIsph::All);
        _apply_homogeneous_neumann_boundary = false;
        break;
      default:
        // apply homogeneous neumann boundary; normal grad p = 0
        filter_op_laplace.setPairYes(PairIsph::Fluid, //- PairIsph::BufferNeumann, 
                                     PairIsph::Fluid);
        _apply_homogeneous_neumann_boundary = true;
        break;
      }

      const int nsize = (_pair->atom->nlocal + _pair->atom->nghost);
      Epetra_SerialDenseVector tmp(View, _w, nsize);
      for (int i=0;i<nsize;++i)
        tmp[i] = 1.0/_rho[i]; 
      
      LaplacianMatrix<PairIsph> op_laplace(_pair, -_dt, tmp.Values());
      op_laplace.setFilter(&filter_op_laplace);
      
      PairFor(op_laplace, op_laplace.getNumberOfWork());
    }

    if (_apply_homogeneous_neumann_boundary) {
      FilterBinary filter_op_grad_dot;
      filter_op_grad_dot.setPairYes(PairIsph::Solid,
                                    PairIsph::All);

      FunctorOuterGradientDotOperatorMatrix<PairIsph,GradientOperator> op_grad_dot(_pair, _normal, -_dt);
      op_grad_dot.setFilter(&filter_op_grad_dot);

      PairFor(op_grad_dot, op_grad_dot.getNumberOfWork());
    }

    _pair->A.crs->ExtractDiagonalCopy(*_pair->A.scaled_laplace_diagonal); 
    
    // op_vstar (div) should include dummy solid particles which may have non-zero
    // velocities
    _filter_op_vstar.setPairYes(PairIsph::Fluid, PairIsph::All);
    _op_vstar.setFilter(&_filter_op_vstar);

    _op_vstar.enterFor();

    _is_once = false;
  }

  template<class PairIsph,
           template<typename> class LaplacianMatrix,
           template<typename> class Divergence,
           template<typename> class GradientOperator> 
  inline void 
  FunctorOuterIncompNavierStokesPoisson<PairIsph,LaplacianMatrix,Divergence,GradientOperator>::operator()(const int ii) {
    const int i = _ilist[ii];
    const int itype = _type[i];
    const int ikind = _pair->getParticleKind(itype);

    // ** modify diagonal and right hand side
    double &diag = (*_pair->A.diagonal)[i];
    
    switch (ikind) {
    // case PairIsph::BufferDummy:
    // case PairIsph::BufferNeumann: {
    //   diag = 1.0;
    //   _b[i] = 0.0;
    //   break;
    // }
    case PairIsph::Solid: {
      if (_apply_homogeneous_neumann_boundary) {
        const double norm = util.dotVectors(_dim, _normal[i], _normal[i]);
        if (norm < 0.5)
          diag = 1.0;
      } else {
        diag = 1.0;
      }
      _b[i] = 0.0;
      break;
    }
      //case PairIsph::BufferDummy:
    case PairIsph::BufferNeumann: 
    case PairIsph::BufferDirichlet: 
    case PairIsph::Fluid: {
      diag = (*_pair->A.scaled_laplace_diagonal)[i]; 

      _op_vstar(ii);
      const double div_at_i = _op_vstar.getDivergence();

      _b[i] = -div_at_i;

      if (_pair->comm->me == 0 && !_is_once) {
        _pair->modifySingularMatrix(_tag[i], diag, _b[i]);
        _is_once = true;
      }
      break;
    }
    default:
      _pair->error->all(FLERR, "FunctorIncompNavierStokesPoisson:: Particle types are not supported");
    }

  }
  template<class PairIsph,
           template<typename> class LaplacianMatrix,
           template<typename> class Divergence,
           template<typename> class GradientOperator> 
  inline void 
  FunctorOuterIncompNavierStokesPoisson<PairIsph,LaplacianMatrix,Divergence,GradientOperator>::exitFor() {
    _pair->A.crs->ReplaceDiagonalValues(*_pair->A.diagonal);
  }

}

#endif

