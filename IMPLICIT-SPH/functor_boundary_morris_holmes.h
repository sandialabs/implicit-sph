#pragma once
#ifndef __FUNCTOR_BOUNDARY_MORRIS_HOLMES_H__
#define __FUNCTOR_BOUNDARY_MORRIS_HOLMES_H__

#include <iostream>
#include <string>
#include "utils.h"
#include "filter.h"
#include "functor.h"

#include "mirror_morris_holmes.h"

namespace LAMMPS_NS {

  namespace Corrected {

    // implementation is due to the paper titled with 
    // "Smooth particle hydrodynamics simulation of low Reynolds number flows through porous media"
    // -- Holmes et al., Int, J. Nuer. Anal. Meth. Geomech. 2000

    // ============================================================================
    template<class PairIsph>
    class FunctorOuterLaplacian_MorrisHolmes : public FunctorOuterLaplacian<PairIsph> {
    public:
      FunctorOuterLaplacian_MorrisHolmes(PairIsph *isph, 
                                         double *f,
                                         double alpha = 1.0, double *scal = NULL)
        : FunctorOuterLaplacian<PairIsph>(isph, f, alpha, scal),
          _pnd(isph->pnd) { }

      FunctorOuterLaplacian_MorrisHolmes(PairIsph *isph, 
                                         double **f3,
                                         double alpha = 1.0, double *scal = NULL)
        : FunctorOuterLaplacian<PairIsph>(isph, f3, alpha, scal),
          _pnd(isph->pnd) { }

      void createMirrorFunction();

    protected:
      FUNCTOR_REMOVE_THIS_POINTER;    
      double *_pnd;
    };
    template<class PairIsph> inline void
    FunctorOuterLaplacian_MorrisHolmes<PairIsph>::createMirrorFunction() {
      this->_mirror = new MirrorMorrisHolmes(_dim,_type,_pair->h,_pair->morris_safe_coeff,_pnd,_vfrac);
    }

    // ============================================================================
    template<class PairIsph, bool AntiSymmetric = false>
    class FunctorOuterLaplacianMatrix_MorrisHolmes : public FunctorOuterLaplacianMatrix<PairIsph,AntiSymmetric> {
    public:
      FunctorOuterLaplacianMatrix_MorrisHolmes(PairIsph *isph, 
                                               double alpha, double *scal = NULL, int iblock = -1, double** normal = NULL)
        : FunctorOuterLaplacianMatrix<PairIsph,AntiSymmetric>(isph, alpha, scal, iblock, normal),
          _pnd(isph->pnd) { }
      void createMirrorFunction();  
    protected:
      FUNCTOR_REMOVE_THIS_POINTER;  
      double *_pnd;
    };
    template<class PairIsph, bool AntiSymmetric> inline void
    FunctorOuterLaplacianMatrix_MorrisHolmes<PairIsph,AntiSymmetric>::createMirrorFunction() {
      this->_mirror = new MirrorMorrisHolmes(_dim,_type,_pair->h,_pair->morris_safe_coeff,_pnd,_vfrac);
    }

    template<class PairIsph> using FunctorOuterLaplacianMatrixSymmetric_MorrisHolmes     = FunctorOuterLaplacianMatrix_MorrisHolmes<PairIsph,false>;
    template<class PairIsph> using FunctorOuterLaplacianMatrixAntiSymmetric_MorrisHolmes = FunctorOuterLaplacianMatrix_MorrisHolmes<PairIsph,true>;

    // ============================================================================
    template<class PairIsph>
    class FunctorOuterGradient_MorrisHolmes : public FunctorOuterGradient<PairIsph> {
    public:
      FunctorOuterGradient_MorrisHolmes(PairIsph *isph,
                                        double *f,
                                        double alpha = 1.0,
                                        double **grad = NULL)
        : FunctorOuterGradient<PairIsph>(isph, f, alpha, grad),
          _pnd(isph->pnd) { }

      FunctorOuterGradient_MorrisHolmes(PairIsph *isph,
                                        double **f3,
                                        double alpha = 1.0)
        : FunctorOuterGradient<PairIsph>(isph, f3, alpha),
          _pnd(isph->pnd) { }

      FunctorOuterGradient_MorrisHolmes(PairIsph *isph,
                                        double *f,
                                        double **f3,
                                        double alpha = 1.0)
        : FunctorOuterGradient<PairIsph>(isph, f, f3, alpha),
          _pnd(isph->pnd) { }

      void createMirrorFunction();

    protected:
      FUNCTOR_REMOVE_THIS_POINTER;  
      double *_pnd;
    };
    template<class PairIsph> inline void
    FunctorOuterGradient_MorrisHolmes<PairIsph>::createMirrorFunction() {
      this->_mirror = new MirrorMorrisHolmes(_dim,_type,_pair->h,_pair->morris_safe_coeff,_pnd,_vfrac);
    }
  
    // ============================================================================
    template<class PairIsph>
    class FunctorOuterGradientOperator_MorrisHolmes : public FunctorOuterGradientOperator<PairIsph> {
    public:
      FunctorOuterGradientOperator_MorrisHolmes(PairIsph *isph, double alpha = 1.0)
        : FunctorOuterGradientOperator<PairIsph>(isph, alpha),
          _pnd(isph->pnd) { }

      void createMirrorFunction();

    protected:
      FUNCTOR_REMOVE_THIS_POINTER;  
      double *_pnd;
    };
    template<class PairIsph> inline void
    FunctorOuterGradientOperator_MorrisHolmes<PairIsph>::createMirrorFunction() {
      this->_mirror = new MirrorMorrisHolmes(_dim,_type,_pair->h,_pair->morris_safe_coeff,_pnd,_vfrac);
    }

    // ============================================================================
    template<class PairIsph, bool AntiSymmetric = false>
    class FunctorOuterDivergence_MorrisHolmes : public FunctorOuterDivergence<PairIsph,AntiSymmetric> {
    public:
      FunctorOuterDivergence_MorrisHolmes(PairIsph *isph, 
                                          double **f,
                                          double alpha = 1.0, double *div = NULL) 
        : FunctorOuterDivergence<PairIsph,AntiSymmetric>(isph, f, alpha, div), 
          _pnd(isph->pnd) { }
      void createMirrorFunction();
    protected:
      FUNCTOR_REMOVE_THIS_POINTER;    
      double *_pnd;
    };
    template<class PairIsph, bool AntiSymmetric> inline void
    FunctorOuterDivergence_MorrisHolmes<PairIsph,AntiSymmetric>::createMirrorFunction() {
      this->_mirror = new MirrorMorrisHolmes(_dim,_type,_pair->h,_pair->morris_safe_coeff,_pnd,_vfrac);
    }

    template<class PairIsph> using FunctorOuterDivergenceSymmetric_MorrisHolmes     = FunctorOuterDivergence_MorrisHolmes<PairIsph,false>;
    template<class PairIsph> using FunctorOuterDivergenceAntiSymmetric_MorrisHolmes = FunctorOuterDivergence_MorrisHolmes<PairIsph,true>;
  }
}

#endif

