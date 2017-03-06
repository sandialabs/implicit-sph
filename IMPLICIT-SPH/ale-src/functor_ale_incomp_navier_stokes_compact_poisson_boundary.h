#pragma once
#ifndef __FUNCTOR_ALE_INCOMP_NAVIER_STOKES_COMPACT_POISSON_BOUNDARY_H__
#define __FUNCTOR_ALE_INCOMP_NAVIER_STOKES_COMPACT_POISSON_BOUNDARY_H__

#include "Epetra_CrsMatrix.h"

#include <iostream>
#include <string>
#include "utils.h"
#include "pair_for.h"
#include "functor.h"

namespace LAMMPS_NS {
  
  template<class PairIsph, 
           template<typename> class LaplacianMatrixCompactPoisson>
  class FunctorOuterAleIncompNavierStokesCompactPoissonBoundary : public FunctorOuter<PairIsph> {
  public:    
    FunctorOuterAleIncompNavierStokesCompactPoissonBoundary(PairIsph *isph, 
                                                            double *rho, 
                                                            double *f, double **normal, double *g,
                                                            double *b) 
      : FunctorOuter<PairIsph>(isph),
        _rho(rho),
        _f(f),
        _normal(normal),
        _g(g),
        _b(b) { }
    
    void enterFor();
    void operator()(const int ii);
    void exitFor();    

  protected:
    FUNCTOR_REMOVE_THIS_POINTER;

    bool _is_once;
    double *_rho, *_f, **_normal, *_g, *_b;
  };
  
  template<class PairIsph,
           template<typename> class LaplacianMatrixCompactPoisson>
  inline void
  FunctorOuterAleIncompNavierStokesCompactPoissonBoundary<PairIsph,
                                                          LaplacianMatrixCompactPoisson>::enterFor() {
    // --- Error check ---------------------------------------------------------------------
    if (_filter != NULL)
      _pair->error->all(FLERR, "FunctorIncompNavierStokesCompactPoissonBoundary:: Filter is not allowed");
    
    if (_pair->boundary_particle == PairIsph::Solid) 
      _pair->error->all(FLERR, "FunctorIncompNavierStokesCompactPoissonBoundary:: SPH is not allowed");
    
    if (_pair->A.is_filled == 1) 
      _pair->error->all(FLERR, "FunctorIncompNavierStokesCompactPoissonBoundary:: A is already filled"); 
    // --- Error check ---------------------------------------------------------------------
    
    _pair->A.crs->PutScalar(0.0);

    {
      const int nsize = (_pair->atom->nlocal + _pair->atom->nghost);
      Epetra_SerialDenseVector material(nsize);
      for (int i=0;i<nsize;++i)
        material[i] = 1.0/_rho[i];

      FilterBinary filter_op_laplace;
      filter_op_laplace.setPairYes(PairIsph::All, PairIsph::All);
      
      LaplacianMatrixCompactPoisson<PairIsph> op_laplace(_pair, 
                                                         -1.0, 
                                                         _f, _normal, _g, 
                                                         _b, 
                                                         material.Values());
      op_laplace.setFilter(&filter_op_laplace);
      
      PairFor(op_laplace, op_laplace.getNumberOfWork());
    }

    _pair->A.crs->ExtractDiagonalCopy(*_pair->A.diagonal); 

    _is_once = false;
  }
  
  template<class PairIsph,
           template<typename> class LaplacianMatrixCompactPoisson>
  inline void 
  FunctorOuterAleIncompNavierStokesCompactPoissonBoundary<PairIsph,
                                                          LaplacianMatrixCompactPoisson>::operator()(const int ii) {
    const int i = _ilist[ii];
    const int itype = _type[i];
    const int ikind = _pair->getParticleKind(itype);
    
    switch (ikind) {
    case PairIsph::Boundary: {
      _b[i] += _f[i];
      break;
    }
    case PairIsph::Fluid: {
      _b[i] += _f[i];

      if (_pair->comm->me == 0 && !_is_once) {
        double &diag = (*_pair->A.diagonal)[i];
        _pair->modifySingularMatrix(_tag[i], diag, _b[i]); 
        _is_once = true;
      }
      break;
    }
    default:
      _pair->error->all(FLERR, "FunctorIncompNavierStokesCompactPoissonBoundary:: Particle types are not supported");
    }
  }

  template<class PairIsph,
           template<typename> class LaplacianMatrixCompactPoisson>
  inline void 
  FunctorOuterAleIncompNavierStokesCompactPoissonBoundary<PairIsph,
                                                          LaplacianMatrixCompactPoisson>::exitFor() {
    _pair->A.crs->ReplaceDiagonalValues(*_pair->A.diagonal);
  }      
  
}

#endif

