#pragma once
#ifndef __FUNCTOR_BOUNDARY_MORRIS_NORMAL_H__
#define __FUNCTOR_BOUNDARY_MORRIS_NORMAL_H__

#include <iostream>
#include <string>
#include "utils.h"
#include "filter.h"
#include "functor.h"

#include "mirror_morris_normal.h"

namespace LAMMPS_NS {

  namespace Corrected {

    template<class PairIsph>
    class FunctorOuterLaplacian_MorrisNormal : public FunctorOuterLaplacian<PairIsph> {
    public:
      FunctorOuterLaplacian_MorrisNormal(PairIsph *isph,
                                         double *f,
                                         double alpha, double *scal = NULL)
        : FunctorOuterLaplacian<PairIsph>(isph, f, alpha, scal),
          _bd_coord(isph->bd_coord), 
          _n(isph->normal) { }
      void createMirrorFunction();
    protected:
      double *_bd_coord, **_n;
      FUNCTOR_REMOVE_THIS_POINTER;    
    };
    template<class PairIsph> inline void
    FunctorOuterLaplacian_MorrisNormal<PairIsph>::createMirrorFunction() {
      this->_mirror = new MirrorMorrisNormal(_dim, _type, _pair->h, _pair->morris_safe_coeff, _x, _bd_coord, _n);
    }
  
    template<class PairIsph>
    class FunctorOuterLaplacianMatrix_MorrisNormal : public FunctorOuterLaplacianMatrix<PairIsph> {
    public:
      FunctorOuterLaplacianMatrix_MorrisNormal(PairIsph *isph,
                                               double alpha, double *scal = NULL, int iblock = -1, double** normal=NULL)
        : FunctorOuterLaplacianMatrix<PairIsph>(isph, alpha, scal, iblock, normal),
          _bd_coord(isph->bd_coord), 
          _n(isph->normal) { }
      void createMirrorFunction();  
    protected:
      double *_bd_coord, **_n ;
      FUNCTOR_REMOVE_THIS_POINTER;  
    };
    template<class PairIsph> inline void
    FunctorOuterLaplacianMatrix_MorrisNormal<PairIsph>::createMirrorFunction() {
      this->_mirror = new MirrorMorrisNormal(_dim, _type, _pair->h, _pair->morris_safe_coeff, _x, _bd_coord, _n);
    }

  }  
}

#endif

