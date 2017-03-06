#pragma once
#ifndef __FUNCTOR_MLS_BOUNDARY_NEUMANN_H__
#define __FUNCTOR_MLS_BOUNDARY_NEUMANN_H__

#include <iostream>
#include <string>
#include <vector>
#include "utils.h"
#include "functor.h"

namespace LAMMPS_NS {

  namespace MLS {

    template<class PairIsph>
    class FunctorOuterBoundaryNeumann : public FunctorOuter<PairIsph> {
    public:
      FunctorOuterBoundaryNeumann(PairIsph *isph, 
                                  double **normal,
                                  double **w, double **vstar,
                                  double gamma,
                                  double *neumann = NULL) 
        : FunctorOuter<PairIsph>(isph),
          _normal(normal),
          _w(w),_vstar(vstar),
          _gamma(gamma),
          _neumann(neumann) { }
      
      double getNeumannBoundary();
      
      void enterFor();
      void operator()(const int ii);

    protected:
      FUNCTOR_REMOVE_THIS_POINTER;

      double **_normal,**_w,**_vstar,_gamma;
      double _neumann_at_i, *_neumann;
    };

    template<class PairIsph> inline double
    FunctorOuterBoundaryNeumann<PairIsph>::getNeumannBoundary() {
      return _neumann_at_i;
    }

    template<class PairIsph> inline void
    FunctorOuterBoundaryNeumann<PairIsph>::enterFor() {
      if (_pair->boundary_particle == PairIsph::Solid)
        _pair->error->all(FLERR, "FunctorBoundaryNeumann:: Solid particles are not allowed");
    }

    template<class PairIsph> inline void
    FunctorOuterBoundaryNeumann<PairIsph>::operator()(const int ii) {
      int i = _ilist[ii];
      int itype = _type[i];
      int ikind = _pair->getParticleKind(itype);

      _neumann_at_i = 0.0;

      if (_filter != NULL && !_filter->yes(ikind))
        return;

      for (int k=0;k<_dim;++k)
        _neumann_at_i += (_w[i][k] - _vstar[i][k])*_normal[i][k];
      _neumann_at_i *= _gamma;
      
      if (_neumann != NULL)
        _neumann[i] = _neumann_at_i;
    }

  }
}
#endif

