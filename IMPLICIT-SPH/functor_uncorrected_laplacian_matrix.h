#pragma once
#ifndef __FUNCTOR_UNCORRECTED_LAPLACIAN_MATRIX_H__
#define __FUNCTOR_UNCORRECTED_LAPLACIAN_MATRIX_H__

#include <iostream>
#include <string>
#include "utils.h"
#include "mirror.h"
#include "functor.h"

namespace LAMMPS_NS {

  namespace Corrected {

    template<class PairIsph>
    class FunctorOuterLaplacianMatrix : public FunctorOuter<PairIsph> {
    public:
      FunctorOuterLaplacianMatrix(PairIsph *isph, 
                                  double alpha, double *scal = NULL) 
        : FunctorOuter<PairIsph>(isph),_mirror(NULL),_alpha(alpha),_scal(scal) { }

      virtual void enterFor();
      virtual void operator()(const int ii);
      virtual void exitFor();

      virtual void createMirrorFunction();

    protected:
      Epetra_IntSerialDenseVector _idx;
      Epetra_SerialDenseVector _val;

      MirrorFunction *_mirror;
      double _alpha, *_scal;

      FUNCTOR_REMOVE_THIS_POINTER;
    };

    template<class PairIsph> inline void
    FunctorOuterLaplacianMatrix<PairIsph>::createMirrorFunction() {
      _mirror = new MirrorNothing(_dim, _type, _pair->h);
    }

    template<class PairIsph> inline void
    FunctorOuterLaplacianMatrix<PairIsph>::enterFor() {
      if (_filter == NULL)
        _pair->error->all(FLERR, "FunctorLaplacianMatrix:: Filter is required to consider dirichlet-like boundary");

      int n = _pair->tags_in_cut->MaxNumIndices();
      _idx.Size(n);
      _val.Size(n);

      createMirrorFunction();
    }

    template<class PairIsph> inline void
    FunctorOuterLaplacianMatrix<PairIsph>::operator()(const int ii) {
      int i = _ilist[ii];
      int itype = _type[i];

      // extract pre-allocated row
      double *val = _val.Values();
      int *idx = _idx.Values();

      // for solid particles, diagonal terms are unit
      // scaling is still applied to them (alpha and scal)
      if (!_filter->yes(_pair->getParticleKind(itype))) {
        idx[0] = _tag[i];

        // when we do not compute for i, leave the diagonal entree zero
        // and later need to change diagonal
        val[0] = 0.0; //_alpha;
        _pair->A.crs->ReplaceGlobalValues(_tag[i], 1, val, idx);  
        return;
      }

      // neighbor around i
      int *jlist = _firstneigh[i];
      int jnum = _numneigh[i];

      // mirroring for particle i
      _mirror->setFreeParticle(i);

      // initialize the correction matrix
      double diag = 0.0;
      int cnt = 0;

      // compute laplacian

      for (int jj=0;jj<jnum;++jj) {
        // j th particle; always put neighmask two bits are reserved for LAMMPS
        int j = (jlist[jj] & NEIGHMASK);
        int jtype = _type[j];
      
        // compute distance between partile i and j
        double rsq = 0.0, rij[3] = {};
        for (int k=0;k<_dim;++k) {
          rij[k] = _x[i][k] - _x[j][k];
          rsq += (rij[k]*rij[k]);
        }
      
        double cutsq = LOOKUP(_pair->cutsq, itype, jtype);
        if (rsq < cutsq) {
          if (_filter->yes(_pair->getParticleKind(itype),
                           _pair->getParticleKind(jtype))) {
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

            double aij = 2.0*dwdr*_vfrac[j]/r * coeff;

            val[cnt] = -aij;
            diag    +=  aij;
            idx[cnt] = _tag[j];

            ++cnt;
          }
        }
      }
      val[cnt] = diag;
      idx[cnt] = _tag[i];
      
      ++cnt;

      _val.Scale(_alpha*(_scal == NULL ? 1.0 : _scal[i]));
      _pair->A.crs->SumIntoGlobalValues(_tag[i], cnt, val, idx);
    }

    template<class PairIsph> inline
    void FunctorOuterLaplacianMatrix<PairIsph>::exitFor() {

      // finalize connectivity construction
      if (_pair->A.is_filled == 0) {
        _pair->A.crs->FillComplete();
        _pair->A.crs->OptimizeStorage();
        _pair->A.is_filled = 1;
      }
      delete _mirror;
    }

  }
}

#endif

