#pragma once
#ifndef __FUNCTOR_MLS_MASS_MATRIX_H__
#define __FUNCTOR_MLS_MASS_MATRIX_H__

#include <iostream>
#include <string>
#include "utils.h"
#include "functor.h"

namespace LAMMPS_NS {

  namespace MLS {

    template<class PairIsph>
    class FunctorOuterMassMatrix : public FunctorOuter<PairIsph> {
    public:
      FunctorOuterMassMatrix(PairIsph *isph) 
        : FunctorOuter<PairIsph>(isph),
          _is_pseudo_inverse_used(false) { }

      void usePseudoInverse(const bool pseudo_inverse);

      void enterFor();
      void operator()(const int ii);

    protected:
      bool _is_pseudo_inverse_used;
      int _lw;

      Epetra_SerialDenseVector _M, _Mfac, _q, _w;

      FUNCTOR_REMOVE_THIS_POINTER;
    };

    template<class PairIsph> inline void 
    FunctorOuterMassMatrix<PairIsph>::usePseudoInverse(const bool pseudo_inverse) {
      _is_pseudo_inverse_used = pseudo_inverse;
    }

    template<class PairIsph> inline void 
    FunctorOuterMassMatrix<PairIsph>::enterFor() {
      // create a local workspace (add one more dof for safety)
      const int ndof = _pair->basis->ndof(ISPH_MLS_MAX_ORDER) + 1;

      // mass matrix 
      _M.Size(ndof*ndof);
      _Mfac.Size(ndof*ndof);

      // base;
      _q.Size(ndof);

      // work space array
      _lw = ndof*max(6, ndof);
      _w.Size(_lw);
    }

    template<class PairIsph> inline void 
    FunctorOuterMassMatrix<PairIsph>::operator()(const int ii) {
      const int i = _ilist[ii];
      const int itype = _type[i];
      const int ikind = _pair->getParticleKind(itype);

      // neighbor around i
      int *jlist = _firstneigh[i];
      int jnum = _numneigh[i];

      // approx order at particle i
      const int np = _pair->np;

      // initialize the mass matrix
      const int ndof  = _pair->basis->ndof(np);
      const int msize = ndof*ndof;

      memset(_pair->M[i], 0, sizeof(double)*msize);
    
      if (_filter != NULL && !_filter->yes(ikind))
        return;

      double *M = _M.Values();
      double *Mfac = _Mfac.Values();
      double *q = _q.Values();
      double *ww = _w.Values();
    
      // clean up local array
      memset(M, 0, sizeof(double)*msize);
      memset(q, 0, sizeof(double)*ndof);

      for (int jj=0;jj<jnum;++jj) {
        const int j = (jlist[jj] & NEIGHMASK);
        const int jtype = _type[j];
        const int jkind = _pair->getParticleKind(jtype);

        if (_filter != NULL && !_filter->yes(ikind,jkind))
          continue;

        // compute distance between partile i and j
        double rsq = 0.0, rij[3] = {};
        for (int k=0;k<_dim;++k) {
          rij[k] = _x[i][k] - _x[j][k];
          rsq += (rij[k]*rij[k]);
        }

        const double r   = sqrt(rsq) + ISPH_EPSILON;
        const double rth = LOOKUP(_pair->h, itype, jtype);
        if (r < rth) {
          // compute kernel
          const double w = _pair->kernel->val(r, rth);

          // ** standard mass matrix
          {
            // compute basis upto order np
            _pair->basis->val(np, rij, rth, q);
            
            // update the symmetric matrix M (upper triangular only)
            // _M += w*q*q^t
            util.updateSymRankDenseMatrix(w,                         // scalar scale
                                          ndof, q,                   // updating vector
                                          Utils::UpperTriangular,    // uplo
                                          ndof, M);                  // matrix
          }
        }
      } // end jj

      // self contribution
      {
        const double r   = 0.0;
        const double rij[3] = {};
        const double rth = LOOKUP(_pair->h, itype, itype);
        
        // compute kernel
        const double w = _pair->kernel->val(r, rth);
        
        // ** standard mass matrix
        {
          // compute basis upto order np
          _pair->basis->val(np, rij, rth, q);
          
          // update the symmetric matrix M (upper triangular only)
          // _M += w*q*q^t
          util.updateSymRankDenseMatrix(w,                           // scalar scale
                                        ndof, q,                     // updating vector
                                        Utils::UpperTriangular,      // uplo
                                        ndof, M);                    // matrix
        }
      }

      // invert M and store it to Ma (M = Ma)
    
      // make a full matrix then use LU
      util.symmetrize(Utils::UpperTriangular, ndof, M);

      // keep M for a debugging purpose
      memcpy(Mfac, M, sizeof(double)*ndof*ndof);

      if (_is_pseudo_inverse_used) 
        util.invertDenseMatrixSVD(ndof, Mfac, _pair->M[i], _lw, ww);
      else
        util.invertDenseMatrix(ndof, Mfac, _pair->M[i], _lw, ww);

      // util.printDenseMatrix("Minv", ndof, ndof, _pair->M[i]);
      // util.printDenseMatrix("M", ndof, ndof, M);
    }

  }
}

#endif

