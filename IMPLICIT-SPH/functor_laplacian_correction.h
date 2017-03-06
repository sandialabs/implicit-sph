#pragma once
#ifndef __FUNCTOR_LAPLACIAN_CORRECTION_H__
#define __FUNCTOR_LAPLACIAN_CORRECTION_H__

#include <iostream>
#include <string>
#include "utils.h"
#include "functor.h"

namespace LAMMPS_NS {

  namespace Corrected {

    template<class PairIsph>
    class FunctorOuterLaplacianCorrection : public FunctorOuter<PairIsph> {
    public:
      FunctorOuterLaplacianCorrection(PairIsph *isph) : FunctorOuter<PairIsph>(isph) { }
      void operator()(const int ii);
    protected:
      FUNCTOR_REMOVE_THIS_POINTER;
    };

    // Reference version
    template<class PairIsph> inline void 
    FunctorOuterLaplacianCorrection<PairIsph>::operator()(const int ii) {
      int i = _ilist[ii];
      int itype = _type[i];

      // neighbor around i
      int *jlist = _firstneigh[i];
      int jnum = _numneigh[i];

      memset(_pair->Lc[i], 0, sizeof(double)*_dimL);

      if (_filter != NULL && !_filter->yes(_pair->getParticleKind(itype)))
        return; 
    
      // initialize the correction matrix
      double A[27] = {};
      double L[36] = {};
      double *G = _pair->Gc[i];

      // compute intermediate tensor A
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

        if (rsq < LOOKUP(_pair->cutsq, itype, jtype)) {
          double r = sqrt(rsq) + ISPH_EPSILON;
          double dwdr = _pair->kernel->dval(r, LOOKUP(_pair->h, itype, jtype));

          // compute intermediate array: 
          // a_{ij}^{k} = \sum_{q} { G_i^{kq} dwdr (r_{ij}^{q}/r) V_{j} }
          //            = \sum_{q} { G_i^{kq} r_{ij}^{q} (dwdr / r V_{j}) }
          double aij[3] = {};
          for (int k2=0;k2<_dim;++k2) {  // k
            for (int k1=0;k1<_dim;++k1)  // q
              aij[k2] += VIEW2(G, _dim, k1, k2)*rij[k1];
            aij[k2] *= dwdr/r*_vfrac[j]; 
          }

          // construct 3rd order tensor A:
          // A_{ij}^{kmn} = \sum_{q} { G_i^{kq} r_{ij}^{m} r_{ij}^{n} dwdr (r_{ij}^{q} / r) V_{j} }
          //              = \sum_{q} { G_i^{kq} dwdr (r_{ij}^{q} / r) V_{j} } r_{ij}^{m} r_{ij}^{n}
          //              = a_{ij}^{k} r_{ij}^{m} r_{ij}^{n}
          for (int k3=0;k3<_dim;++k3) {        // k
            // slice the tensor to view in 2D matrix
            double *slice = &A[k3*_dimsq];
            for (int k2=0;k2<_dim;++k2)        // n
              for (int k1=0;k1<(k2+1);++k1)   // m
                VIEW2(slice, _dim, k1, k2) += aij[k3]*rij[k1]*rij[k2];
          }
        }
      }
      // compute correction tensor 
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
      
        // integrate 
        if (rsq < LOOKUP(_pair->cutsq, itype, jtype)) {
          double r = sqrt(rsq) + ISPH_EPSILON;
          double dwdr = _pair->kernel->dval(r, LOOKUP(_pair->h, itype, jtype));

          // compute unit vector e_{ij}
          double eij[3];
          for (int k=0;k<_dim;++k)
            eij[k] = (rij[k]/r);

          // compute intermediate value C; use half
          // C_{ij}^{mn} = \sum_{k} { A_{i}^{kmn} e_{ij}^{k} + r_{ij}^{m} e_{ij}^{n} } dwdr V_{j}
          double C[9] = {}; 
          for (int k3=0;k3<_dim;++k3) {         // k
            double *slice = &A[k3*_dimsq];          
            for (int k2=0;k2<_dim;++k2)         // n
              for (int k1=0;k1<(k2+1);++k1)    // m
                VIEW2(C, _dim, k1, k2) += VIEW2(slice, _dim, k1, k2)*eij[k3];
          }
          // rank one update with scaling
          for (int k2=0;k2<_dim;++k2)          // n
            for (int k1=0;k1<(k2+1);++k1) {   // m
              VIEW2(C, _dim, k1, k2) += rij[k1]*eij[k2];
              VIEW2(C, _dim, k1, k2) *= dwdr*_vfrac[j];
            }
        
          // compute a linear system L
          // L_{ij}^{mn,op} = \sum_{k} { A_{i}^{kmn} e_{ij}^{k} + r_{ij}^{m} e_{ij}^{n} }
          //                  e_{ij}^{o} dwdr (r_{ij}^{p} / r) V_{j}
          //                = C_{ij}^{mn} e_{ij}^{o} e_{ij}^{p} 
          //
          // consider symmetry in the summation
          // if o == p, L_{ij}^{oo,mn} = C_{ij}^{mn}  e_{ij}^{o} e_{ij}^{o}
          // if o != p, L_{ij}^{op,mn} = C_{ij}^{mn} (e_{ij}^{o} e_{ij}^{p}) * 2.0
          const double scale[2] = {2.0, 1.0};
          for (int k4=0, op=0;k4<_dim;++k4)           // p
            for (int k3=0;k3<(k4+1);++k3,++op)        // o
              for (int k2=0, mn=0;k2<_dim;++k2)       // n
                for (int k1=0;k1<(k2+1);++k1,++mn)   // m
                  VIEW2(L, _dimL, mn, op) += VIEW2(C, _dim, k1, k2)*eij[k3]*eij[k4]*scale[k3==k4];
        }
      }
    
      // set rhs
      for (int k2=0, op=0;k2<_dim;++k2)  
        for (int k1=0;k1<(k2+1);++k1,++op) 
          _pair->Lc[i][op] = -double(k1 == k2);

      // solve L for Lc (rhs) and store the solution to Lc
      double ipiv[6];
      util.solveDenseMatrix(Utils::Full,
                            _dimL, _dimL, L, 
                            _dimL,     1, _pair->Lc[i],
                            6, &ipiv[0]); 
    }

  }
}

#endif

