#pragma once
#ifndef __FUNCTOR_APPLIED_ELECTRIC_POTENTIAL_H__
#define __FUNCTOR_APPLIED_ELECTRIC_POTENTIAL_H__

#include <iostream>
#include <string>
#include "utils.h"
#include "pair_for.h"
#include "functor.h"

namespace LAMMPS_NS {

  template<class PairIsph,
           template<typename> class LaplacianMatrix>
  class FunctorOuterAppliedElectricPotential : public FunctorOuter<PairIsph> {
  public:

    FunctorOuterAppliedElectricPotential(PairIsph *isph,
                                         double *sigma, double *phi, double **normal,
                                         double *b)
      : FunctorOuter<PairIsph>(isph),
        _sigma(sigma), _phi(phi), _normal(normal),
        _b(b) { }
    
    void enterFor();
    void operator()(const int ii);
    void exitFor();

  protected:
    double *_sigma, *_phi, **_normal, *_b;
    FUNCTOR_REMOVE_THIS_POINTER;
  };
  
  template<class PairIsph,
           template<typename> class LaplacianMatrix>
  inline void
  FunctorOuterAppliedElectricPotential<PairIsph,LaplacianMatrix>::enterFor() {
    // --- Error check --------------------------------------------------------------------- 
    if (_filter != NULL)
      _pair->error->all(FLERR, "FunctorAppliedElectricPotential:: Filter is not allowed");

    if (_pair->A.is_filled == 1) 
      _pair->error->all(FLERR, "FunctorAppliedElectricPotential:: A is already filled");      
    // --- Error check --------------------------------------------------------------------- 

    _pair->A.crs->PutScalar(0.0);

    // compute laplacin matrix
    {
      // buffer plays like a dirichlet condition
      FilterMatchBinary filter;
      filter.setPairYes(PairIsph::Fluid, PairIsph::Fluid);

      LaplacianMatrix<PairIsph> functor(_pair, -1.0, _sigma);
      functor.setFilter(&filter);
    
      PairFor(functor, functor.getNumberOfWork());
    }
    _pair->A.crs->ExtractDiagonalCopy(*_pair->A.diagonal);
  }

  template<class PairIsph,
           template<typename> class LaplacianMatrix>
  inline void 
    FunctorOuterAppliedElectricPotential<PairIsph,LaplacianMatrix>::operator()(const int ii) {
    const int i = _ilist[ii];
    const int itype = _type[i];
    const int ikind = _pair->getParticleKind(itype);

    // ** modify diagonal and right hand side
    double &diag = (*_pair->A.diagonal)[i];

    _b[i] = 0.0; 
    switch (ikind) {
    case PairIsph::Fluid: 
      break;
    case PairIsph::Solid:
      diag = 1.0;
      break;
    case PairIsph::BufferNeumann: 
    case PairIsph::BufferDirichlet: 
      diag = 1.0; _b[i] = _phi[i];
      break;
    default:
      break;
    }
  }
  
  template<class PairIsph,
           template<typename> class LaplacianMatrix>
  inline void 
  FunctorOuterAppliedElectricPotential<PairIsph,LaplacianMatrix>::exitFor() {
    _pair->A.crs->ReplaceDiagonalValues(*_pair->A.diagonal);
  }
}

#endif

