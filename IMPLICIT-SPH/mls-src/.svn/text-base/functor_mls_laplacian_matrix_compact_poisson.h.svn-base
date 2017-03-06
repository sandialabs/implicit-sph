#pragma once
#ifndef __FUNCTOR_MLS_LAPLACIAN_MATRIX_COMPACT_POISSON_H__
#define __FUNCTOR_MLS_LAPLACIAN_MATRIX_COMPACT_POISSON_H__

#include <iostream>
#include <string>
#include "utils.h"
#include "functor.h"

namespace LAMMPS_NS {

  namespace MLS {

    template<class PairIsph>
    class FunctorOuterLaplacianMatrixCompactPoisson : public FunctorOuter<PairIsph> {
    public:
      FunctorOuterLaplacianMatrixCompactPoisson(PairIsph *isph,
                                                double alpha, 
                                                double *f,
                                                double **normal, double *g,
                                                double *b,
                                                double *material = NULL)
        : FunctorOuter<PairIsph>(isph),
          _op_laplace_a(isph, alpha, material),
          _op_laplace_b(isph, NULL, alpha, f, normal, g, material, b)
      { }
      
      void enterFor();
      void operator()(const int ii);
      void exitFor();

    protected:
      FUNCTOR_REMOVE_THIS_POINTER;

      FunctorOuterLaplacianMatrix<PairIsph> _op_laplace_a;
      FunctorOuterLaplacianCompactPoisson<PairIsph> _op_laplace_b;
    };

    template<class PairIsph> inline void
    FunctorOuterLaplacianMatrixCompactPoisson<PairIsph>::enterFor() {
      if (_filter == NULL)
        _pair->error->all(FLERR, "FunctorLaplacianMatrixCompactPoisson:: Filter is required to consider dirichlet-like boundary");

      _op_laplace_a.setFilter(_filter);
      _op_laplace_b.setFilter(_filter);

      _op_laplace_a.enterFor();
      _op_laplace_b.enterFor();
    }

    template<class PairIsph> inline void
    FunctorOuterLaplacianMatrixCompactPoisson<PairIsph>::operator()(const int ii) {
      _op_laplace_a(ii);
      _op_laplace_b(ii);
    }

    template<class PairIsph> inline
    void FunctorOuterLaplacianMatrixCompactPoisson<PairIsph>::exitFor() {
      _op_laplace_a.exitFor();
      _op_laplace_b.exitFor();
    }

  }
}
#endif

