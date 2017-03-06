#pragma once
#ifndef __FUNCTOR_MLS_HELPER_COMPACT_POISSON_H__
#define __FUNCTOR_MLS_HELPER_COMPACT_POISSON_H__

#include <iostream>
#include <string>
#include "utils.h"
#include "functor.h"

namespace LAMMPS_NS {

  namespace MLS {

    // this compute temporary vector for evluating \sum P(x_j) W_ij f_j
    template<class PairIsph>
    class FunctorOuterHelperCompactPoisson : public FunctorOuter<PairIsph> {
    public:
      FunctorOuterHelperCompactPoisson(PairIsph *isph, 
                                       double *f, 
                                       double **normal, double *g)
        : FunctorOuter<PairIsph>(isph),
          _f(f),
          _f3(NULL),
          _normal(normal),
          _g(g),
          _g3(NULL),
          _is_entered(false) { }

      FunctorOuterHelperCompactPoisson(PairIsph *isph, 
                                       double **f3, 
                                       double **normal, double **g3)
        : FunctorOuter<PairIsph>(isph),
          _f(NULL),
          _f3(f3),
          _normal(normal),
          _g(NULL),
          _g3(g3),
          _is_entered(false) { }

      FunctorOuterHelperCompactPoisson(PairIsph *isph, 
                                       double *f, double **f3,
                                       double **normal, double *g, double **g3) 
        : FunctorOuter<PairIsph>(isph),
          _f(f),
          _f3(f3),
          _normal(normal),
          _g(g),
          _g3(g3),
          _is_entered(false) { }

      virtual void enterFor();
      virtual void operator()(const int ii);

      double* getScalarHelper();
      double* getVectorHelper(const int k);

      bool hasScalarField() const;
      bool hasVectorField() const;
      
    protected:
      FUNCTOR_REMOVE_THIS_POINTER;

      // function value
      double *_f, **_f3;

      // boundary value
      double **_normal, *_g, **_g3;

      int _lda;
      Epetra_SerialDenseVector _q, _dq, _q_at_i, _q3_at_i;

      bool _is_entered;
    };

    template<class PairIsph> inline double*
    FunctorOuterHelperCompactPoisson<PairIsph>::getScalarHelper() {
      return &_q_at_i[0];
    }
    template<class PairIsph> inline double*
    FunctorOuterHelperCompactPoisson<PairIsph>::getVectorHelper(const int k) {
      return &_q3_at_i[k*_lda];
    }
    template<class PairIsph> inline bool
    FunctorOuterHelperCompactPoisson<PairIsph>::hasScalarField() const {
      return (_f != NULL);
    }
    template<class PairIsph> inline bool
    FunctorOuterHelperCompactPoisson<PairIsph>::hasVectorField() const {
      return (_f3 != NULL);
    }

    template<class PairIsph> inline void
    FunctorOuterHelperCompactPoisson<PairIsph>::enterFor() {
      if (!_is_entered) {
        // ndof of basis functions plus the dof for the lagrange multiplier
        _lda = _pair->basis->ndof(ISPH_MLS_MAX_ORDER) + 1;
        
        // workspace for j
        _q.Size(_lda);
        _dq.Size(_lda);
        
        // scalar field workspace for i
        if (_f != NULL)
          _q_at_i.Size(_lda);
        
        // vector field workspace for i
        if (_f3 != NULL) 
          _q3_at_i.Size(_lda*_dim);
        
        _is_entered = true;
      }
    }

    template<class PairIsph> inline void
    FunctorOuterHelperCompactPoisson<PairIsph>::operator()(const int ii) {
      if (!_is_entered)
        enterFor();

      const int i = _ilist[ii];
      const int itype = _type[i];
      const int ikind = _pair->getParticleKind(itype);

      // basis dof
      const int np   = _pair->np;
      const int ndof = _pair->basis->ndof(np);
      
      // neighbor around i
      int *jlist = _firstneigh[i];
      int jnum = _numneigh[i];

      // for interior penalty: laplacian
      const int p2[3][3] = { {2,0,0},    // dx2
                             {0,2,0},    // dy2
                             {0,0,2} };  // dz2

      // for boundary penalty: gradient normal
      const int p1[3][3] = { {1,0,0},    // dx
                             {0,1,0},    // dy
                             {0,0,1} };  // dz

      // pointers to workspace q and dq
      double *q  = _q.Values();
      double *dq = _dq.Values();

      // basis
      double *q_at_i = _q_at_i.Values();
      double *q3_at_i = _q3_at_i.Values();

      // initialization for early return
      if (_f != NULL) 
        memset(q_at_i, 0, sizeof(double)*_lda);

      if (_f3 != NULL) 
        memset(q3_at_i, 0, sizeof(double)*_lda*_dim);
      
      if (_filter != NULL && !_filter->yes(ikind))
        return;
      
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

          // ** compute \laplace P W f_j
          {
            const double c_interior = _pair->cp.interior_penalty_coeff(rth);
            memset(dq, 0, sizeof(double)*ndof);
            if (c_interior > ISPH_EPSILON) {
              const bool update = true;
              for (int k=0;k<_dim;++k)
                _pair->basis->dval(np, rij, rth,   p2[k][0], p2[k][1], p2[k][2], 1.0, dq, update);
            }

            // scalar field 
            if (_f != NULL) 
              for (int k=0;k<ndof;++k)
                q_at_i[k] += c_interior*dq[k]*w*_f[j];
            
            // vector field; component-wise
            if (_f3 != NULL) 
              for (int k2=0;k2<_dim;++k2) 
                for (int k1=0;k1<ndof;++k1) 
                  VIEW2(q3_at_i, _lda, k1, k2) += c_interior*dq[k1]*w*_f3[j][k2];
          }

          // compute gradient basis
          if (jkind == PairIsph::Boundary) {
            const double c_boundary = _pair->cp.boundary_penalty_coeff(rth);
            memset(dq, 0, sizeof(double)*ndof);
            if (c_boundary > ISPH_EPSILON) {
              const bool update = true;              
              for (int k=0;k<_dim;++k)
                _pair->basis->dval(np, rij, rth,   p1[k][0], p1[k][1], p1[k][2], _normal[j][k], dq, update);
              
              // scalar field 
              if (_g != NULL) 
                for (int k=0;k<ndof;++k)
                  q_at_i[k] += c_boundary*dq[k]*w*_g[j];
              
              // vector field; component-wise
              if (_g3 != NULL) 
                for (int k2=0;k2<_dim;++k2) 
                  for (int k1=0;k1<ndof;++k1) 
                    VIEW2(q3_at_i, _lda, k1, k2) += c_boundary*dq[k1]*w*_g3[j][k2];
            }
          }
        }
      }

      // penalty
      {
        const double r   = 0.0;
        const double rth = LOOKUP(_pair->h, itype, itype); 

        // compute kernel
        const double w = _pair->kernel->val(r, rth);

        int idx_at_i[3];
        double dq_at_i[3];
        
        {
          const double c_interior = _pair->cp.interior_penalty_coeff(rth);
          if (c_interior > ISPH_EPSILON) {
            for (int k=0;k<_dim;++k)
              _pair->basis->dval(np, rth,   p2[k][0], p2[k][1], p2[k][2],   idx_at_i[k], dq_at_i[k]);

            // scalar field 
            if (_f != NULL) 
              for (int k=0;k<_dim;++k)
                q_at_i[idx_at_i[k]] += c_interior*dq_at_i[k]*w*_f[i];
            
            // vector field; component-wise
            if (_f3 != NULL) 
              for (int k2=0;k2<_dim;++k2) 
                for (int k1=0;k1<_dim;++k1) 
                  VIEW2(q3_at_i, _lda, idx_at_i[k1], k2) += c_interior*dq[k1]*w*_f3[i][k2];
          }
        }

        if (ikind == PairIsph::Boundary) {
          const double c_boundary = _pair->cp.boundary_penalty_coeff(rth);
          if (c_boundary > ISPH_EPSILON) {
            for (int k=0;k<_dim;++k) {
              _pair->basis->dval(np, rth,   p1[k][0], p1[k][1], p1[k][2],   idx_at_i[k], dq_at_i[k]);
              dq_at_i[k] *= _normal[i][k];
            }
            
            // scalar field
            if (_g != NULL)
              for (int k=0;k<_dim;++k)
                q_at_i[idx_at_i[k]] += c_boundary*dq_at_i[k]*w*_g[i];
            
            // vector field; component-wise
            if (_g3 != NULL)
              for (int k2=0;k2<_dim;++k2)
                for (int k1=0;k1<_dim;++k1)
                  VIEW2(q3_at_i, _lda, idx_at_i[k1], k2) += c_boundary*dq_at_i[k1]*w*_g3[i][k2];
          }
          
          // apply the constraint
          if (_g != NULL)
            q_at_i[ndof] = _g[i];
          
          if (_g3 != NULL) 
            for (int k=0;k<_dim;++k) 
              VIEW2(q3_at_i, _lda, ndof, k) = _g3[i][k];
        }
      }
    }
  }
}
#endif

