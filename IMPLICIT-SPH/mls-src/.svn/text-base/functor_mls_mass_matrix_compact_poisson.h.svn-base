#pragma once
#ifndef __FUNCTOR_MLS_MASS_MATRIX_COMPACT_POISSON_H__
#define __FUNCTOR_MLS_MASS_MATRIX_COMPACT_POISSON_H__

#include <iostream>
#include <string>
#include "utils.h"
#include "functor.h"

namespace LAMMPS_NS {

  namespace MLS {
    
    template<class PairIsph>
    class FunctorOuterMassMatrixCompactPoisson : public FunctorOuter<PairIsph> {
    public:
      FunctorOuterMassMatrixCompactPoisson(PairIsph *isph, 
                                           double **normal) 
        : FunctorOuter<PairIsph>(isph), 
          _normal(normal),
          _is_pseudo_inverse_used(false) { }

      void usePseudoInverse(const bool pseudo_inverse);

      void enterFor();
      void operator()(const int ii);

    protected:
      bool _is_pseudo_inverse_used;
      int _lw;
      double **_normal;
      Epetra_SerialDenseVector _M, _Mfac, _q, _dq, _w;

      FUNCTOR_REMOVE_THIS_POINTER;
    };

    template<class PairIsph> inline void 
    FunctorOuterMassMatrixCompactPoisson<PairIsph>::usePseudoInverse(const bool pseudo_inverse) {
      _is_pseudo_inverse_used = pseudo_inverse;
    }

    template<class PairIsph> inline void 
    FunctorOuterMassMatrixCompactPoisson<PairIsph>::enterFor() {
      // create a local workspace that includes the lagrange 
      // multipler and an additional dof for safety 
      const int ndofpl = _pair->basis->ndof(ISPH_MLS_MAX_ORDER) + 2;

      // mass matrix
      _M.Size(ndofpl*ndofpl);
      _Mfac.Size(ndofpl*ndofpl);

      // base;
      _q.Size(ndofpl);
      _dq.Size(ndofpl);

      _lw = ndofpl*max(6, ndofpl);
      _w.Size(_lw);
    }

    template<class PairIsph> inline void 
    FunctorOuterMassMatrixCompactPoisson<PairIsph>::operator()(const int ii) {
      const int i = _ilist[ii];
      const int itype = _type[i];
      const int ikind = _pair->getParticleKind(itype);

      // neighbor around i
      int *jlist = _firstneigh[i];
      int jnum = _numneigh[i];

      // approx order at particle i
      const int np = _pair->np;

      // initialize the mass matrix
      const int ndof   = _pair->basis->ndof(np);
      const int ndofpl = ndof + (ikind == PairIsph::Boundary); // ndof plus lagrange multiplier
      const int msize  = ndofpl*ndofpl;

      memset(_pair->M[i], 0, sizeof(double)*msize);
      
      if (_filter != NULL && !_filter->yes(ikind))
        return;

      double *M = _M.Values();
      double *Mfac = _Mfac.Values();

      // for standard mass matrix
      double *q = _q.Values();

      // diffuse derivative: laplacian and gradient
      double *dq = _dq.Values();

      // work array
      double *ww = _w.Values();

      // for interior penalty: laplacian
      const int p2[3][3] = { {2,0,0},    // dx2
                             {0,2,0},    // dy2
                             {0,0,2} };  // dz2

      // for boundary penalty: gradient normal
      const int p1[3][3] = { {1,0,0},    // dx
                             {0,1,0},    // dy
                             {0,0,1} };  // dz

      // clean up local array
      memset(M,   0, sizeof(double)*msize);

      for (int jj=0;jj<jnum;++jj) {
        const int j = (jlist[jj] & NEIGHMASK);
        const int jtype = _type[j];
        const int jkind = _pair->getParticleKind(jtype);

        if (_filter != NULL && !_filter->yes(ikind, jkind))
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
          
          // construction of mass matrix (mass + laplacian penalty + boundary penalty)
          // - note that M size ndofpl and updateSymRank will update top left part
          //   for interior particles

          // ** standard mass matrix
          {
            // compute basis upto order np
            memset(q, 0, sizeof(double)*ndofpl);
            _pair->basis->val(np, rij, rth, q);
            
            // update the symmetric matrix M: M += w*q*q^t
            util.updateSymRankDenseMatrix(w,                           // scalar scale
                                          ndof, q,                   // updating vector
                                          Utils::UpperTriangular,      // uplo
                                          ndofpl, M);                  // matrix
          }

          // ** penalize laplacian for all nodes
          {
            const double c_interior = _pair->cp.interior_penalty_coeff(rth);
            memset(dq, 0, sizeof(double)*ndofpl);
            if (c_interior > ISPH_EPSILON) {
              // need to accumulate contribution of the multi-index representation
              const bool update = true;

              // compute derivatives to construct laplacian penalty upto order np
              for (int k=0;k<_dim;++k)
                _pair->basis->dval(np, rij, rth,   p2[k][0], p2[k][1], p2[k][2], 1.0, dq, update);
           
              // update the symmetric matrix M: M = c_interior*w*dq*dq^t
              util.updateSymRankDenseMatrix(c_interior*w,                // scalar scale 
                                            ndof, dq,                    // updating vector
                                            Utils::UpperTriangular,      // uplo
                                            ndofpl, M);                  // matrix
            }
          }
          
          // ** penalize poisson boundary condition
          if (jkind == PairIsph::Boundary) {
            const double c_boundary = _pair->cp.boundary_penalty_coeff(rth);
            memset(dq, 0, sizeof(double)*ndofpl);
            if (c_boundary > ISPH_EPSILON) {
              // need to accumulate contribution of the multi-index representation
              const bool update = true;
              
              // compute derivative to account for the boundary penalty upto order np
              for (int k=0;k<_dim;++k) 
                _pair->basis->dval(np, rij, rth,   p1[k][0], p1[k][1], p1[k][2], _normal[j][k], dq, update);
              
              // update the symmetric matrix M: M = c_boundary*w*dq*dq^t
              util.updateSymRankDenseMatrix(c_boundary*w,                // scalar scale 
                                            ndof, dq,                    // updating vector
                                            Utils::UpperTriangular,      // uplo
                                            ndofpl, M);                  // matrix
            }
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
        
        int idx_at_i[3];
        double dq_at_i[3];

        // ** standard mass matrix
        {
          // compute basis upto order np
          memset(q, 0, sizeof(double)*ndofpl);
          _pair->basis->val(np, rij, rth, q);
          
          // update the symmetric matrix M: M += w*q*q^t
          util.updateSymRankDenseMatrix(w,                           // scalar scale
                                        ndof, q,                     // updating vector
                                        Utils::UpperTriangular,      // uplo
                                        ndofpl, M);                  // matrix
        }
       
        // ** penalize laplacian for all nodes
        {
          const double c_interior = _pair->cp.interior_penalty_coeff(rth);
          if (c_interior > ISPH_EPSILON) {
            for (int k=0;k<_dim;++k)
              _pair->basis->dval(np, rth,   p2[k][0], p2[k][1], p2[k][2],   idx_at_i[k], dq_at_i[k]);
            
            // update the symmetric matrix M: M = c_interior*w*dq*dq^t
            util.updateSymRankDenseMatrix(c_interior*w,                // scalar scale
                                          _dim, idx_at_i, dq_at_i,     // updating sparse vector
                                          Utils::UpperTriangular,      // uplo
                                          ndofpl, M);                  // matrix
          }
        }
        
        // ** penalize poisson boundary condition
        if (ikind == PairIsph::Boundary) {
          for (int k=0;k<_dim;++k) {
            _pair->basis->dval(np, rth,   p1[k][0], p1[k][1], p1[k][2],   idx_at_i[k], dq_at_i[k]);
            
            // compute the normal component
            dq_at_i[k] *= _normal[i][k];

            // apply the constraint on the upper triangular
            VIEW2(M, ndofpl, idx_at_i[k], ndof) = dq_at_i[k];
          }

          // update the symmetric matrix M: M = c_boundary*w*dq*dq^t          
          const double c_boundary = _pair->cp.boundary_penalty_coeff(rth);
          if (c_boundary > ISPH_EPSILON) 
            util.updateSymRankDenseMatrix(c_boundary*w,                // scalar scale
                                          _dim, idx_at_i, dq_at_i,     // updating sparse vector
                                          Utils::UpperTriangular,      // uplo
                                          ndofpl, M);                  // matrix
        } 
      }

      // invert M and store it to Mb (M = Mb)
    
      // make a full matrix
      util.symmetrize(Utils::UpperTriangular, ndofpl, M);

      memcpy(Mfac, M, sizeof(double)*ndofpl*ndofpl);

      if (_is_pseudo_inverse_used)
        util.invertDenseMatrixSVD(ndofpl, Mfac, _pair->M[i], _lw, ww);
      else
        util.invertDenseMatrix(ndofpl, Mfac, _pair->M[i], _lw, ww);        

      // util.printDenseMatrix("Minv", ndof, ndof, _pair->M[i]);
      // util.printDenseMatrix("M", ndofpl, ndofpl, M);
    }

  }
}

#endif

