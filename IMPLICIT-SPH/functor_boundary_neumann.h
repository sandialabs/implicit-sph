#pragma once
#ifndef __FUNCTOR_BOUNDARY_NEUMANN_H__
#define __FUNCTOR_BOUNDARY_NEUMANN_H__

#include <iostream>
#include <string>
#include <vector>
#include "utils.h"
#include "functor.h"

namespace LAMMPS_NS {

  namespace Corrected {
    
    template<class PairIsph>
    class FunctorOuterBoundaryNeumann : public FunctorOuterLaplacianHelper<PairIsph> {
    public:
      FunctorOuterBoundaryNeumann(PairIsph *isph, 
                                  double **normal,
                                  double **w, double **vstar,
                                  double gamma = 1.0) 
        : FunctorOuterLaplacianHelper<PairIsph>(isph, gamma, NULL, true, false),
          _normal(normal),
          _w(w),_vstar(vstar),
          _gamma(gamma) { }
      
      // overload the diff function
      double diff(const int i, const int j);    

      void enterFor();

      double getNeumannBoundary();
      
    protected:
      FUNCTOR_REMOVE_THIS_POINTER;
      
      double **_normal,**_w,**_vstar,_gamma;
    };

    template<class PairIsph> inline double
    FunctorOuterBoundaryNeumann<PairIsph>::getNeumannBoundary() {
      return this->getScalarLaplacian();
    }

    template<class PairIsph> inline void
    FunctorOuterBoundaryNeumann<PairIsph>::enterFor() {
      if (_pair->boundary_particle == PairIsph::Boundary)
        _pair->error->all(FLERR, "FunctorBoundaryNeumann:: Boundary particles are not allowed");
    }
    
    template<class PairIsph> inline double
    FunctorOuterBoundaryNeumann<PairIsph>::diff(const int i, const int j) {
      // compute normal distance : rn = rij dot n
      // n and rij are defined from solid to fluid
      double rn = 0.0;
      for (int k=0;k<_dim;++k)
        rn += (_x[i][k] - _x[j][k])*_normal[i][k];

      // make sure distance is positive however normal is defined
      rn = abs(rn);

      // check rn is safe to compute
      if (rn < ISPH_EPSILON) 
        return 0.0;

      // compute difference between boundary value and extrapolated velocity
      // val = (w - vstar) dot n
      double val = 0.0;    
      for (int k=0;k<_dim;++k)
        val += (_w[i][k] - _vstar[i][k])*_normal[i][k];

      // - dt /rho (p_i - p_j) = gamma_0 * rn * val
      return (_gamma*rn*val);
    }

  }
}
#endif

