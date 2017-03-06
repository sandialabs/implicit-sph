#pragma once
#ifndef __FUNCTOR_POISSON_BOLTZMAN_JACOBIAN_H__
#define __FUNCTOR_POISSON_BOLTZMAN_JACOBIAN_H__

#include <iostream>
#include <string>
#include "utils.h"
#include "pair_for.h"
#include "functor.h"
#include "functor_laplacian_matrix.h"

namespace LAMMPS_NS {

  template<class PairIsph,
           template<typename> class LaplacianMatrix>
  class FunctorOuterPoissonBoltzmannJacobian : public FunctorOuter<PairIsph> {
  public:
    FunctorOuterPoissonBoltzmannJacobian(PairIsph *isph, 
                                         double *psi, double *eps)
      : FunctorOuter<PairIsph>(isph),
        _psi(psi), _eps(eps) { }
    void enterFor();
    void operator()(const int ii);
    void exitFor();

  protected:
    double *_pnd, *_bd_coord, **_normal,*_psi, *_eps;
 
    bool _is_linearized_pb;   
    double _kappasq,_gamma;

    FUNCTOR_REMOVE_THIS_POINTER;
  };
  
  template<class PairIsph,
           template<typename> class LaplacianMatrix>
  inline void
  FunctorOuterPoissonBoltzmannJacobian<PairIsph,LaplacianMatrix>::enterFor() {
    if (_filter != NULL)
      _pair->error->all(FLERR, "FunctorPoissonBoltzmannJacobian:: Filter is not allowed");

    _is_linearized_pb = _pair->pb.is_linearized;

    _kappasq = 2.0*_pair->pb.ezcb/_pair->pb.psiref;

	_gamma = _pair->pb.gamma;

    // compute laplace jacobian and store the diagonals separately and reuse jacobian
    // for non-linear iteration
    if (_pair->A.is_filled == 0) {
      _pair->A.crs->PutScalar(0.0);

      FilterBinary filter_op_laplace;
      filter_op_laplace.setPairYes(PairIsph::Fluid, PairIsph::All);    

      // construct operator with -1.0 scaling
      LaplacianMatrix<PairIsph> op_laplace(_pair, -1.0, _eps);
      op_laplace.setFilter(&filter_op_laplace);
      
      // construct matrix
      PairFor(op_laplace, op_laplace.getNumberOfWork());
      
      _pair->A.crs->ExtractDiagonalCopy(*_pair->A.scaled_laplace_diagonal);   
      _pair->A.is_filled = 1;
    }
  }

  template<class PairIsph,
           template<typename> class LaplacianMatrix>
  inline void
  FunctorOuterPoissonBoltzmannJacobian<PairIsph,LaplacianMatrix>::operator()(const int ii) {
    const int i = _ilist[ii];
    const int itype = _type[i];
    const int ikind = _pair->getParticleKind(itype);

    // minus unit matrix
    switch (ikind) {
    case PairIsph::Solid:
    case PairIsph::Boundary:
      (*_pair->A.diagonal)[i] = -1.0;
      break;
    case PairIsph::BufferDirichlet:
    case PairIsph::BufferNeumann:
    case PairIsph::Fluid:
      // Poisson Boltzmann J = \minu-laplacian-jac{psi} + cosh{psi} 
      (*_pair->A.diagonal)[i] = (*_pair->A.scaled_laplace_diagonal)[i];
      if (_is_linearized_pb) {
        const double psi = _psi[i];
        const double numerator   = 4.0 - 2.0*_gamma*pow(psi, 2);
        const double denominator = pow(_gamma,2)*pow(psi,4) + 4.0*_gamma*pow(psi,2) + 4.0;
        (*_pair->A.diagonal)[i] += _kappasq*(numerator/denominator);
      } else {
        const double psi = _psi[i];
        const double numerator   = 2.0*_gamma*cosh(0.5*psi)*sinh(0.5*psi)*sinh(psi);
        const double denominator = 2.0*_gamma*pow(sinh(0.5*psi),2) + 1.0;
        (*_pair->A.diagonal)[i] += _kappasq*(cosh(psi)/denominator - numerator/pow(denominator,2));
      }
      break;
    }
  }

  template<class PairIsph,
           template<typename> class LaplacianMatrix>
  inline void
  FunctorOuterPoissonBoltzmannJacobian<PairIsph,LaplacianMatrix>::exitFor() {
    _pair->A.crs->ReplaceDiagonalValues(*_pair->A.diagonal);
  }
}

#endif

