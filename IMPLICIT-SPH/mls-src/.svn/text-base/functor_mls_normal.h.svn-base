#pragma once
#ifndef __FUNCTOR_MLS_NORMAL_H__
#define __FUNCTOR_MLS_NORMAL_H__

#include <iostream>
#include <string>
#include "utils.h"
#include "functor.h"

namespace LAMMPS_NS {

  namespace MLS {
    
    template<class PairIsph>
    class FunctorOuterNormal : public FunctorOuter<PairIsph> {
    public:
      FunctorOuterNormal(PairIsph *isph, 
                         double *orient, double **normal)
        : FunctorOuter<PairIsph>(isph),
          _lda(0),
          _orient(orient),
          _normal(normal) { }
      
      void enterFor();
      void operator()(const int ii);
      
      double* getNormal();
      
    protected:
      FUNCTOR_REMOVE_THIS_POINTER;
      
      int _lda;
      Epetra_SerialDenseVector _q, _q_at_i;
      
      double *_orient,**_normal;
      double _normal_at_i[3];
    };
    
    template<class PairIsph> inline double*
    FunctorOuterNormal<PairIsph>::getNormal() {
      return _normal_at_i;
    }
    
    template<class PairIsph> inline void
    FunctorOuterNormal<PairIsph>::enterFor() {
      if (_filter == NULL)
        _pair->error->all(FLERR, "FunctorNormal:: Filter is required which defines interface of different particles");
      
      _lda = _pair->basis->ndof(ISPH_MLS_MAX_ORDER);
      
      _q.Size(_lda);
      _q_at_i.Size(_lda);      
    }

    template<class PairIsph> inline void
    FunctorOuterNormal<PairIsph>::operator()(const int ii) {
      int i = _ilist[ii];
      int itype = _type[i];
      int ikind = _pair->getParticleKind(itype);

      // initialize with zero
      memset(_normal_at_i, 0, sizeof(double)*3);
      if (!_filter->yes(ikind))
        return;
      
      // basis dof
      int np   = _pair->np;
      int ndof = _pair->basis->ndof(np);

      // neighbor around i
      int *jlist = _firstneigh[i];
      int jnum = _numneigh[i];

      // pointers to q
      double *q = _q.Values();
      double *q_at_i = _q_at_i.Values();

      // initialization for early return
      memset(q_at_i, 0, sizeof(double)*_lda);

      for (int jj=0;jj<jnum;++jj) {
        int j = (jlist[jj] & NEIGHMASK);
        int jtype = _type[j];
        int jkind = _pair->getParticleKind(jtype);

        if (!_filter->yes(ikind, jkind)) 
          continue;

        // compute distance between partile i and j
        double rsq = 0.0, rij[3] = {};
        for (int k=0;k<_dim;++k) {
          rij[k] = _x[i][k] - _x[j][k];
          rsq += (rij[k]*rij[k]);
        }

        double r   = sqrt(rsq) + ISPH_EPSILON;
        double rth = LOOKUP(_pair->h, itype, jtype); 
        if (r < rth) {
          // compute kernel
          double w = _pair->kernel->val(r, rth);

          // compute basis
          _pair->basis->val(np, rij, rth, q); 

          // scalar field 
          for (int k=0;k<ndof;++k)
            q_at_i[k] += q[k]*w;
        }
      }

      // basis arguments
      double rth = LOOKUP(_pair->h, itype, itype);

      const int p[3][3] = { {1,0,0},    // dx
                            {0,1,0},    // dy
                            {0,0,1} };  // dz

      // mass matrix
      double *M = _pair->M[i];

      // orientation
      double orient = _orient[ikind];

      // compute gradient
      for (int k2=0;k2<_dim;++k2) {
        int idx = 0;
        double dq = 0.0, *row = NULL;

        // d/d{x,y,z}, pick a row correspond to the combination of p
        _pair->basis->dval(np, rth,   p[k2][0], p[k2][1], p[k2][2],   idx, dq); 
        row = &VIEW2(M, ndof, idx, 0);
        
        _normal_at_i[k2] = dq*util.dotVectors(ndof, row, ndof, q_at_i, 1)*orient;
      }

      double alpha = 0.0;
      for (int k=0;k<_dim;++k)
        alpha += _normal_at_i[k]*_normal_at_i[k];
      
      alpha = sqrt(alpha);

      if (alpha != 0.0)
        for (int k = 0; k < _dim; ++k)
          _normal_at_i[k] /= alpha;

      if (_normal != NULL)
        memcpy(_normal[i], _normal_at_i, sizeof(double)*_dim);
    }

  }
}
#endif

