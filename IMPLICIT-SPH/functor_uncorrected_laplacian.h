#pragma once
#ifndef __FUNCTOR_UNCORRECTED_LAPLACIAN_H__
#define __FUNCTOR_UNCORRECTED_LAPLACIAN_H__

#include <iostream>
#include <string>
#include "utils.h"
#include "mirror.h"
#include "functor.h"

namespace LAMMPS_NS {

  namespace Corrected {

    template<class PairIsph>
    class FunctorOuterLaplacian : public FunctorOuter<PairIsph> {
    public:
      FunctorOuterLaplacian(PairIsph *isph, double *f, 
                            double alpha, double *laplace = NULL) 
        : FunctorOuter<PairIsph>(isph),_mirror(NULL),_f(f),_alpha(alpha),_laplace(laplace) { }
      virtual~FunctorOuterLaplacian() {
        if (_mirror != NULL)
          delete _mirror;
      }

      void operator()(const int ii);
      double getScalarLaplacian();

      virtual void createMirrorFunction(); 
    protected:
      MirrorFunction *_mirror;

      double *_f, _alpha, *_laplace, _laplace_at_i;

      FUNCTOR_REMOVE_THIS_POINTER;
    };

    template<class PairIsph> inline double
    FunctorOuterLaplacian<PairIsph>::getScalarLaplacian() {
      return _laplace_at_i;
    }
    template<class PairIsph> inline void
    FunctorOuterLaplacian<PairIsph>::createMirrorFunction() {
      _mirror = new MirrorNothing(_dim, _type, _pair->h);
    }
  
    template<class PairIsph> inline void
    FunctorOuterLaplacian<PairIsph>::operator()(const int ii) {
      if (_mirror == NULL) 
        createMirrorFunction();

      int i = _ilist[ii];
      int itype = _type[i];

      // neighbor around i
      int *jlist = _firstneigh[i];
      int jnum = _numneigh[i];

      _laplace_at_i = 0.0;

      if (_filter != NULL && !_filter->yes(_pair->getParticleKind(itype)))
        return; 

      // mirroring for particle i
      _mirror->setFreeParticle(i);
    
      {
        // compute laplacian 
        for (int jj=0;jj<jnum;++jj) {
          // j th particle; always put neighmask two bits are reserved for LAMMPS
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
            if(_pair->getParticleKind(jtype) == PairIsph::Solid)
              coeff=_mirror->computeMirrorCoefficient(sqrt(cutsq));
          
            double r = sqrt(rsq) + ISPH_EPSILON;
            double dwdr = _pair->kernel->dval(r, LOOKUP(_pair->h, itype, jtype));
          
            // compute unit vector e_{ij}
            double eij[3];
            for (int k=0;k<_dim;++k)
              eij[k] = rij[k]/r;
          

            _laplace_at_i += 2.0*dwdr*_vfrac[j]*coeff*(_f[i] - _f[j])/r;
          }
        }
      }
      _laplace_at_i *= _alpha;
    
      if (_laplace != NULL) 
        _laplace[i] = _laplace_at_i;
    }

  }
}

#endif

