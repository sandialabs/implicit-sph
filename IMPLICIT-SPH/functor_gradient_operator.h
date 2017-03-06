#pragma once
#ifndef __FUNCTOR_GRADIENT_OPERATOR_H__
#define __FUNCTOR_GRADIENT_OPERATOR_H__

#include <iostream>
#include <string>
#include "utils.h"
#include "mirror.h"
#include "functor.h"

namespace LAMMPS_NS {

  namespace Corrected {
  
    template<class PairIsph>
    class FunctorOuterGradientOperator : public FunctorOuter<PairIsph> {
    public:
      FunctorOuterGradientOperator(PairIsph *isph, double alpha = 1.0)
        : FunctorOuter<PairIsph>(isph),
          _mirror(NULL),
          _alpha(alpha),
          _cnt(0),
          _lda(0),
          _is_entered(false) { }
      virtual~FunctorOuterGradientOperator() { 
        if (_mirror != NULL)
          delete _mirror; 
      }

      void enterFor();
      void operator()(const int ii);
      
      int getSize() const;
      int getNeighbor(const int jj);
      double getValue(const int jj, const int idx);
      
      virtual void createMirrorFunction(); 
    protected:
      FUNCTOR_REMOVE_THIS_POINTER;

      MirrorFunction *_mirror; 

      double _alpha;

      int _cnt, _lda;
      Epetra_IntSerialDenseVector _idx;
      Epetra_SerialDenseVector _val;

      bool _is_entered;
    };

    template<class PairIsph> 
    inline int
    FunctorOuterGradientOperator<PairIsph>::getSize() const {
      return _cnt;
    }
    template<class PairIsph> 
    inline int
    FunctorOuterGradientOperator<PairIsph>::getNeighbor(const int jj) {
      return _idx[jj];
    }
    template<class PairIsph> 
    inline double
    FunctorOuterGradientOperator<PairIsph>::getValue(const int jj, const int k) {
      return _val[jj + k*_lda];
    }
    template<class PairIsph> 
    inline void
    FunctorOuterGradientOperator<PairIsph>::createMirrorFunction() {
      _mirror = new MirrorNothing(_dim, _type, _pair->h);
    }

    template<class PairIsph> 
    inline void
    FunctorOuterGradientOperator<PairIsph>::enterFor() {
      if (!_is_entered) {
        if (_mirror == NULL) 
          createMirrorFunction();

        _lda = _pair->tags_in_cut->MaxNumIndices();

        _idx.Size(_lda);
        _val.Size(_lda*3);

        _is_entered = true;
      }
    }

    template<class PairIsph> 
    inline void
    FunctorOuterGradientOperator<PairIsph>::operator()(const int ii) {
      if (!_is_entered) 
        enterFor();

      int i = _ilist[ii];
      int itype = _type[i];

      // extract pointers
      double *val = _val.Values();
      int *idx = _idx.Values();

      // clean up the work arrays
      memset(idx, 0, sizeof(int)*_lda);
      memset(val, 0, sizeof(double)*_lda*3);

      _cnt = 0;

      if (_filter != NULL && !_filter->yes(_pair->getParticleKind(itype)))
        return;

      // neighbor around i
      int *jlist = _firstneigh[i];
      int jnum = _numneigh[i];
      
      // mirroring for particle i
      _mirror->setFreeParticle(i); 

      // initialize the correction matrix
      double *G = _pair->Gc[i];

      // reserve self contribution at idx = 0
      idx[_cnt++] = i; 

      for (int jj=0;jj<jnum;++jj) {
        int j = (jlist[jj] & NEIGHMASK);
        int jtype = _type[j];

        if (_filter != NULL && !_filter->yes(_pair->getParticleKind(itype), 
                                             _pair->getParticleKind(jtype)))
          continue;
      
        // compute distance between partile i and j
        double rsq = 0.0, rij[3] = {};
        for (int k=0;k<_dim;++k) {
          rij[k] = _x[i][k] - _x[j][k];
          rsq += (rij[k]*rij[k]);
        }

        double cutsq = LOOKUP(_pair->cutsq, itype, jtype); 
        if (rsq < cutsq) {
          _mirror->setMirrorParticle(j);
          double coeff = 1.0;
          if(_pair->getParticleKind(itype) != PairIsph::Solid &&
             _pair->getParticleKind(jtype) == PairIsph::Solid)
            coeff = _mirror->computeMirrorCoefficient(sqrt(cutsq));

          double r = sqrt(rsq) + ISPH_EPSILON;
          double dwdr = _pair->kernel->dval(r, LOOKUP(_pair->h, itype, jtype));
        
          // 1st order consistent corrected gradient
          // \grad_{1} psi_{ij}^{k} = (psi_{j} - psi_{i}) G_i^{kq} dwdr (r_{ij}^{q} / r) V_{j}
          //                        = (psi_{j} - psi_{i}) ( G_i^{kq} r_{ij}^{q} ) dwdr / r V_{j} 
          double vjtmp = dwdr/r*_vfrac[j]*coeff;          
          for (int k2=0;k2<_dim;++k2) {   // k
            double gitmp = 0.0;
            for (int k1=0;k1<_dim;++k1)   // q
              gitmp += VIEW2(G, _dim, k1, k2)*rij[k1];

            double ijtmp = gitmp*vjtmp;
            
            VIEW2(val, _lda, _cnt, k2)  = ijtmp;
            VIEW2(val, _lda,    0, k2) -= ijtmp;
          }
          idx[_cnt++] = j; 
        }
      }
    
      // scale with alpha
      _val.Scale(_alpha);
    }

  }
}

#endif

