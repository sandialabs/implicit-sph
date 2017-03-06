#pragma once
#ifndef __FUNCTOR_ALE_INCOMP_NAVIER_STOKES_POISSON_BOUNDARY_H__
#define __FUNCTOR_ALE_INCOMP_NAVIER_STOKES_POISSON_BOUNDARY_H__

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
           template<typename> class GradientOperator,
           template<typename> class NeumannBoundary>
  class FunctorOuterAleIncompNavierStokesPoissonBoundary : public FunctorOuter<PairIsph> {
  public:    
    FunctorOuterAleIncompNavierStokesPoissonBoundary(PairIsph *isph, 
                                                     double dt, double gamma,
                                                     double **normal,
                                                     double *rho, 
                                                     double **w, double **vstar, 
                                                     double *b) 
      : FunctorOuter<PairIsph>(isph),
        _op_vstar(isph, vstar, gamma),            
        _op_neumann(isph, normal, w, vstar, gamma),
        _dt(dt),
        _normal(normal),
        _rho(rho),
        _b(b) { }
    
    void enterFor();
    void operator()(const int ii);
    void exitFor();    

  protected:
    FUNCTOR_REMOVE_THIS_POINTER;

    // for ith particle, the operator computes: 
    // \gamma \div v
    Divergence<PairIsph> _op_vstar;            
    FilterBinary _filter_op_vstar;

    // for ith particle, the operator computes:
    // \sum_j \gamma ((w[i] - v[j]) \dot n[i]) 
    // - w is given velocity in the boundary region
    // - v is predicted velocity
    //     note: particles shuld be computed in the boundary (ghost) region
    // - n normal vector at i
    NeumannBoundary<PairIsph> _op_neumann;
    FilterBinary _filter_op_neumann;
    
    bool _is_once;
    double _dt,**_normal,*_rho,*_b;
  };
  
  template<class PairIsph,
           template<typename> class LaplacianMatrix,
           template<typename> class Divergence,
           template<typename> class GradientOperator,
           template<typename> class NeumannBoundary>
  inline void
  FunctorOuterAleIncompNavierStokesPoissonBoundary<PairIsph,
                                                   LaplacianMatrix,
                                                   Divergence,
                                                   GradientOperator,
                                                   NeumannBoundary>::enterFor() {
    // --- Error check ---------------------------------------------------------------------
    if (_filter != NULL)
      _pair->error->all(FLERR, "FunctorIncompNavierStokesPoissonBoundary:: Filter is not allowed");
    
    if (_pair->boundary_particle == PairIsph::Solid) 
      _pair->error->all(FLERR, "FunctorIncompNavierStokesPoissonBoundary:: SPH is not allowed");
    
    if (_pair->A.is_filled == 1) 
      _pair->error->all(FLERR, "FunctorIncompNavierStokesPoissonBoundary:: A is already filled"); 
    // --- Error check ---------------------------------------------------------------------
    
    _pair->A.crs->PutScalar(0.0);
    
    {
      FilterBinary filter_op_laplace;
      filter_op_laplace.setPairYes (PairIsph::Fluid, 
                                    PairIsph::Fluid | PairIsph::Boundary);      
      
      LaplacianMatrix<PairIsph> op_laplace(_pair, -_dt);
      op_laplace.setFilter(&filter_op_laplace);
      
      PairFor(op_laplace, op_laplace.getNumberOfWork());

      _filter_op_vstar.setPairYes  (PairIsph::Fluid, 
                                    PairIsph::Fluid | PairIsph::Boundary);
      _op_vstar.setFilter(&_filter_op_vstar);      
      
    }
    if (_pair->boundary_particle == PairIsph::Boundary) {
      FilterBinary filter_op_grad_dot;
      filter_op_grad_dot.setPairYes (PairIsph::Boundary, 
                                     PairIsph::Fluid | PairIsph::Boundary);
      
      FunctorOuterGradientDotOperatorMatrix<PairIsph,GradientOperator> op_grad_dot(_pair, _normal, -_dt);
      op_grad_dot.setFilter(&filter_op_grad_dot);
      
      PairFor
        (op_grad_dot, op_grad_dot.getNumberOfWork());

      _filter_op_neumann.setPairYes(PairIsph::Boundary);
      _op_neumann.setFilter(&_filter_op_neumann);

      _op_neumann.enterFor();
    }

    _pair->A.crs->ExtractDiagonalCopy(*_pair->A.diagonal); 

    _is_once = false;
  }
  
  template<class PairIsph,
           template<typename> class LaplacianMatrix,
           template<typename> class Divergence,
           template<typename> class GradientOperator,
           template<typename> class NeumannBoundary>
  inline void 
  FunctorOuterAleIncompNavierStokesPoissonBoundary<PairIsph,
                                                   LaplacianMatrix,
                                                   Divergence,
                                                   GradientOperator,
                                                   NeumannBoundary>::operator()(const int ii) {
    int i = _ilist[ii];
    int itype = _type[i];
    int ikind = _pair->getParticleKind(itype);
    
    switch (ikind) {
    case PairIsph::Boundary: {
      _op_neumann(ii);
      double neumann_at_i =  _op_neumann.getNeumannBoundary();
      
      _b[i] = (_rho[i]*neumann_at_i);
      break;
    }
    case PairIsph::Fluid: {
      _op_vstar(ii);
      double div_at_i = _op_vstar.getDivergence();

      _b[i] = (-_rho[i]*div_at_i);

      if (_pair->comm->me == 0 && !_is_once) {
        double &diag = (*_pair->A.diagonal)[i];
        _pair->modifySingularMatrix(_tag[i], diag, _b[i]); 
        _is_once = true;
      }
      break;
    }
    default:
      _pair->error->all(FLERR, "FunctorIncompNavierStokesPoissonBoundary:: Particle types are not supported");
    }
  }

  template<class PairIsph,
           template<typename> class LaplacianMatrix,
           template<typename> class Divergence,
           template<typename> class GradientOperator,
           template<typename> class NeumannBoundary>
  inline void 
  FunctorOuterAleIncompNavierStokesPoissonBoundary<PairIsph,
                                                   LaplacianMatrix,
                                                   Divergence,
                                                   GradientOperator,
                                                   NeumannBoundary>::exitFor() {
    _pair->A.crs->ReplaceDiagonalValues(*_pair->A.diagonal);
  }      

}

#endif

