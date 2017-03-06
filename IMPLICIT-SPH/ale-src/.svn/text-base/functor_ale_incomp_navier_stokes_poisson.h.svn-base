#pragma once
#ifndef __FUNCTOR_ALE_INCOMP_NAVIER_STOKES_POISSON_H__
#define __FUNCTOR_ALE_INCOMP_NAVIER_STOKES_POISSON_H__

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
           template<typename> class NeumannBoundary>
  class FunctorOuterAleIncompNavierStokesPoisson : public FunctorOuter<PairIsph> {
  public:
    
    FunctorOuterAleIncompNavierStokesPoisson(PairIsph *isph, 
                                             double dt, double gamma,
                                             double **normal,
                                             double *rho, 
                                             double **w, double **vstar, 
                                             double *b) 
      : FunctorOuter<PairIsph>(isph),
        _op_vstar(isph, vstar, gamma),            
        _op_neumann(isph, normal, w, vstar, gamma),
        _dt(dt),_rho(rho),
        _normal(normal),
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
    // \sum_j \gamma ((w[j] - v[j]) \dot n[i]) (r_{ij} \dot n[i])
    // - w is given velocity in the boundary region
    // - v is predicted velocity
    //     note: particles shuld be computed in the boundary (ghost) region
    // - n normal vector at i
    NeumannBoundary<PairIsph> _op_neumann;
    FilterBinary _filter_op_neumann;
    
    bool _is_once;
    double _dt,*_rho,**_normal,*_b;
  };
  

  template<class PairIsph,
           template<typename> class LaplacianMatrix,
           template<typename> class Divergence,
           template<typename> class NeumannBoundary>
  inline void
  FunctorOuterAleIncompNavierStokesPoisson<PairIsph,
                                           LaplacianMatrix,
                                           Divergence,
                                           NeumannBoundary>::enterFor() {
    // --- Error check ---------------------------------------------------------------------
    if (_filter != NULL)
      _pair->error->all(FLERR, "FunctorIncompNavierStokesPoisson:: Filter is not allowed");

    if (_pair->boundary_particle == PairIsph::Boundary) 
      _pair->error->all(FLERR, "FunctorIncompNavierStokesPoisson:: MLS is not allowed");

    if (_pair->A.is_filled == 1) 
      _pair->error->all(FLERR, "FunctorIncompNavierStokesPoisson:: A is already filled"); 
    // --- Error check ---------------------------------------------------------------------    

    _pair->A.crs->PutScalar(0.0);
    
    FilterBinary filter_op_laplace;
    filter_op_laplace.setPairYes (PairIsph::Fluid, PairIsph::Fluid);      
    
    LaplacianMatrix<PairIsph> op_laplace(_pair, -_dt);
    op_laplace.setFilter (&filter_op_laplace);
    
    PairFor(op_laplace, op_laplace.getNumberOfWork());
    _pair->A.crs->ExtractDiagonalCopy(*_pair->A.scaled_laplace_diagonal); 
    
    _filter_op_neumann.setPairYes(PairIsph::Fluid, PairIsph::Solid);
    _filter_op_vstar.setPairYes  (PairIsph::Fluid, PairIsph::Fluid);
    
    _op_neumann.setFilter(&_filter_op_neumann);
    _op_vstar.setFilter  (&_filter_op_vstar);      

    _op_neumann.enterFor();
    
    _is_once = false;
  }

  template<class PairIsph,
           template<typename> class LaplacianMatrix,
           template<typename> class Divergence,
           template<typename> class NeumannBoundary>
  inline void 
  FunctorOuterAleIncompNavierStokesPoisson<PairIsph,
                                           LaplacianMatrix,
                                           Divergence,
                                           NeumannBoundary>::operator()(const int ii) {
    int i = _ilist[ii];
    int itype = _type[i];
    int ikind = _pair->getParticleKind(itype);

    // modify diagonal for solid particles
    double &diag = (*_pair->A.diagonal)[i];

    switch (ikind) {
    case PairIsph::Solid: {
      diag = -1.0; 
      _b[i] = 0.0;
      break;
    }
    case PairIsph::Fluid: {
      diag = (*_pair->A.scaled_laplace_diagonal)[i]; 

      _op_vstar(ii);
      double div_at_i = _op_vstar.getDivergence();

      double neumann_at_i = 0.0;

      bool is_neumann = util.dotVectors(_dim, _normal[i], _normal[i]) > 0.5;
      if (is_neumann) {
        _op_neumann(ii);
        neumann_at_i = _op_neumann.getNeumannBoundary();
      }
      
      _b[i] = -_rho[i]*(div_at_i + neumann_at_i);
      
      if (_pair->comm->me == 0 && !is_neumann && !_is_once) {
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
           template<typename> class NeumannBoundary>
  inline void 
  FunctorOuterAleIncompNavierStokesPoisson<PairIsph,
                                           LaplacianMatrix,
                                           Divergence,
                                           NeumannBoundary>::exitFor() {
    _pair->A.crs->ReplaceDiagonalValues(*_pair->A.diagonal);
  }
  
}

#endif

