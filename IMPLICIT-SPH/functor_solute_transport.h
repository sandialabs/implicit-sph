#pragma once
#ifndef __FUNCTOR_SOLUTE_TRANSPORT_H__
#define __FUNCTOR_SOLUTE_TRANSPORT_H__

#include <iostream>
#include <string>
#include "utils.h"
#include "pair_for.h"
#include "functor.h"

namespace LAMMPS_NS {

  template<class PairIsph,
           template<typename> class LaplacianMatrix,
           template<typename> class GradientOperator>
  class FunctorOuterSoluteTransport : public FunctorOuter<PairIsph> {
  public:

    FunctorOuterSoluteTransport(PairIsph *isph, 
                                double **v, bool is_particle_fixed,
                                double dt, double theta, 
                                double dcoeff, 
                                double *b, double *w)
      : FunctorOuter<PairIsph>(isph),
        _v(v),_is_particle_fixed(is_particle_fixed),
        _dt(dt),_theta(theta),
        _dcoeff(dcoeff),
        _b(b),  // b is assumed to have v^{n}
        _w(w) { }
    
    void enterFor();
    void operator()(const int ii);
    void exitFor();

  protected:
    double **_v;
    bool _is_particle_fixed;

    double _dt,_theta,_dcoeff;
    double *_b,*_w;

    FUNCTOR_REMOVE_THIS_POINTER;
  };
  
  template<class PairIsph,
           template<typename> class LaplacianMatrix,
           template<typename> class GradientOperator>
  inline void
  FunctorOuterSoluteTransport<PairIsph,LaplacianMatrix,GradientOperator>::enterFor() {
    // --- Error check --------------------------------------------------------------------- 
    if (_filter != NULL)
      _pair->error->all(FLERR, "FunctorSoluteTransport:: Filter is not allowed");

    if (_pair->A.is_filled == 1) 
      _pair->error->all(FLERR, "FunctorSoluteTransport:: A is already filled");      
    // --- Error check --------------------------------------------------------------------- 

    _pair->A.crs->PutScalar(0.0);

    // compute laplacin matrix
    {
      FilterMatchBinary filter;
      filter.setPairYes(PairIsph::Fluid, PairIsph::Fluid - PairIsph::BufferNeumann);

      LaplacianMatrix<PairIsph> functor(_pair, _dt*_dcoeff);
      functor.setFilter(&filter);
    
      PairFor(functor, functor.getNumberOfWork());
    }

    // do not think about the case where fluid particles are fixed
    // if (_is_particle_fixed) {
    //   FilterMatchBinary filter;
    //   filter.setPairYes(PairIsph::Fluid, PairIsph::Fluid);
      
    //   FunctorOuterGradientDotOperatorMatrix<PairIsph,GradientOperator> functor(_pair, _v, -_dt);
    //   functor.setFilter(&filter);

    //   PairFor(functor, functor.getNumberOfWork());
    // }
    
    // b and workspace
    Epetra_Vector b(View, *_pair->nodalmap, _b);
    Epetra_Vector w(View, *_pair->nodalmap, _w);
    
    // w = A u, where u is be and A is laplacian dt nu
    _pair->A.crs->Multiply(false, b, w);
    w.Scale(1.0 - _theta);
    
    // complete left hand side
    _pair->A.crs->Scale(-_theta);
    _pair->A.crs->ExtractDiagonalCopy(*_pair->A.scaled_laplace_diagonal);   
  }

  template<class PairIsph,
           template<typename> class LaplacianMatrix,
           template<typename> class GradientOperator>
  inline void 
  FunctorOuterSoluteTransport<PairIsph,LaplacianMatrix,GradientOperator>::operator()(const int ii) {
    int i = _ilist[ii];
    int itype = _type[i];
    int ikind = _pair->getParticleKind(itype);

    // ** modify diagonal
    double &diag = (*_pair->A.diagonal)[i];

    switch (ikind) {
    case PairIsph::BufferDirichlet: 
    case PairIsph::BufferNeumann: 
    case PairIsph::Solid: {
      diag = 1.0;
      break;
    }
    case PairIsph::Fluid: {
      // compute diagonal (I + (-theta dt dcoeff \laplace))
      diag = 1.0 + (*_pair->A.scaled_laplace_diagonal)[i];
      
      // b += w
      _b[i] += _w[i];
      break;      
    }
    default:
      _pair->error->all(FLERR, "FunctorSoluteTransport:: Particle types are not supported");     
    }
  }
  
  template<class PairIsph,
           template<typename> class LaplacianMatrix,
           template<typename> class GradientOperator>
  inline void 
  FunctorOuterSoluteTransport<PairIsph,LaplacianMatrix,GradientOperator>::exitFor() {
    _pair->A.crs->ReplaceDiagonalValues(*_pair->A.diagonal);
  }
}

#endif

