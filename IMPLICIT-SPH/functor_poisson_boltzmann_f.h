#pragma once
#ifndef __FUNCTOR_POISSON_BOLTZMAN_F_H__
#define __FUNCTOR_POISSON_BOLTZMAN_F_H__

#include <iostream>
#include <string>
#include "utils.h"

namespace LAMMPS_NS {

  template<class PairIsph,
           template<typename> class Laplacian>
  class FunctorOuterPoissonBoltzmannF : public FunctorOuter<PairIsph> {
  public:
    FunctorOuterPoissonBoltzmannF(PairIsph *isph, 
                                  double *psi, double *psi0, double *eps,
                                  Epetra_Vector *f) 
      : FunctorOuter<PairIsph>(isph),
        _op_laplace(isph, psi, -1.0, eps),
        _psi(psi),_psi0(psi0), _f(f) { };
    void enterFor();
    void operator()(const int ii);
    
  protected:
    Laplacian<PairIsph> _op_laplace;
    FilterBinary _filter_op_laplace;

    double *_psi,*_psi0;
    Epetra_Vector *_f; 

    bool _is_linearized_pb;
    double _kappasq,_gamma;

    FUNCTOR_REMOVE_THIS_POINTER;
  };

  template<class PairIsph,
           template<typename> class Laplacian>
  inline void
  FunctorOuterPoissonBoltzmannF<PairIsph,Laplacian>::enterFor() {
    if (_filter != NULL)
      _pair->error->all(FLERR, "FunctorPoissonBoltzmannF:: Filter is not allowed");
    _is_linearized_pb = _pair->pb.is_linearized;

    _kappasq = 2.0*_pair->pb.ezcb/_pair->pb.psiref;

    _gamma = _pair->pb.gamma;

    // default boundary
    _filter_op_laplace.setPairYes(PairIsph::Fluid,
                                  PairIsph::All);
    _op_laplace.setFilter(&_filter_op_laplace);
  }        

  template<class PairIsph,
           template<typename> class Laplacian>
  inline void 
  FunctorOuterPoissonBoltzmannF<PairIsph,Laplacian>::operator()(const int ii) {
    int i = _ilist[ii];
    int itype = _type[i];
    int ikind = _pair->getParticleKind(itype);

    switch (ikind) {
    case PairIsph::Solid:
    case PairIsph::Boundary: {
      // psi0 should be given as a normalized value
      (*_f)[i] = (-_psi[i] + _psi0[i]);
      break;
    }
    case PairIsph::BufferDirichlet:
    case PairIsph::BufferNeumann:
    case PairIsph::Fluid: {
      _op_laplace(ii);
      double minus_laplace_at_i = _op_laplace.getScalarLaplacian();
      
      // Poisson Boltzmann F = \laplacian{psi} + sinh{psi} 
      (*_f)[i] = minus_laplace_at_i;
      if (_is_linearized_pb)
        (*_f)[i] += _kappasq*(_psi[i]/(1.0 + 2.0*_gamma*pow(_psi[i]/2,2)));
      else
        (*_f)[i] += _kappasq*(sinh(_psi[i])/(1.0 + 2.0*_gamma*pow(sinh(_psi[i]/2.0),2)));
      break;
    }
    }
  }
  
}

#endif

