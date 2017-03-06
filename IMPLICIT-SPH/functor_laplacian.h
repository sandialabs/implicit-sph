#pragma once
#ifndef __FUNCTOR_LAPLACIAN_H__
#define __FUNCTOR_LAPLACIAN_H__

#include <iostream>
#include <string>
#include "utils.h"
#include "mirror.h"
#include "functor.h"

namespace LAMMPS_NS {

  namespace Corrected {

    // compute laplacian with given difference operator
    // =============================================================================
    template<class PairIsph>
    class FunctorOuterLaplacianHelper : public FunctorOuter<PairIsph> {
    public:
      FunctorOuterLaplacianHelper(PairIsph *isph, 
                                  double alpha, double *material,
                                  bool scalar_diff, bool vector_diff)
        : FunctorOuter<PairIsph>(isph),
          _mirror(NULL),
          _alpha(alpha),
          _material(material),
          _scalar_diff(scalar_diff),
          _vector_diff(vector_diff) { }

      virtual~FunctorOuterLaplacianHelper() {
        if (_mirror != NULL)
          delete _mirror;
      }
      
      void operator()(const int ii);
      double getScalarLaplacian();
      double* getVectorLaplacian();
      
      virtual void createMirrorFunction();

      virtual double diff(const int i, const int j);
      virtual double diff(const int i, const int j, const int k);

    protected:
      FUNCTOR_REMOVE_THIS_POINTER;

      MirrorFunction *_mirror;
      
      double _alpha, *_material;
      double _laplace_at_i,_laplace3_at_i[3];
      
      bool _scalar_diff,_vector_diff;

      Epetra_SerialDenseVector _aij;
    };

    template<class PairIsph> inline double
    FunctorOuterLaplacianHelper<PairIsph>::getScalarLaplacian() {
      return _laplace_at_i;
    }
    template<class PairIsph> inline double*
    FunctorOuterLaplacianHelper<PairIsph>::getVectorLaplacian() {
      return _laplace3_at_i;
    }
    template<class PairIsph> inline void
    FunctorOuterLaplacianHelper<PairIsph>::createMirrorFunction() {
      _mirror = new MirrorNothing(_dim, _type, _pair->h);
    }
    template<class PairIsph> inline double
    FunctorOuterLaplacianHelper<PairIsph>::diff(const int i, const int j) {
      return 0.0;
    }
    template<class PairIsph> inline double
    FunctorOuterLaplacianHelper<PairIsph>::diff(const int i, const int j, const int k) {
      return 0.0;
    }
  
    template<class PairIsph> inline void
    FunctorOuterLaplacianHelper<PairIsph>::operator()(const int ii) {
      if (_mirror == NULL) 
        createMirrorFunction();
      
      const int i = _ilist[ii];
      const int itype = _type[i];
      const int ikind = _pair->getParticleKind(itype);
      const double material_at_i = (_material == NULL ? 1.0 : _material[i]);

      // neighbor around i
      const int *jlist = _firstneigh[i];
      const int jnum = _numneigh[i];

      if (_aij.Length() < jnum) 
        _aij.Resize(jnum);

      // initialization
      _laplace_at_i = 0.0;
      memset(_laplace3_at_i, 0, sizeof(double)*3);

      if (_filter != NULL && !_filter->yes(ikind))
        return; 

      // mirroring for particle i
      _mirror->setFreeParticle(i);

      double grad_material_at_i[3] = {};
      double grad_at_i[3] = {};
      double grad3_at_i[3][3] = {};
      {
        // initialize the correction matrix
        const double *G = _pair->Gc[i];
      
        for (int jj=0;jj<jnum;++jj) {
          const int j = (jlist[jj] & NEIGHMASK);
          const int jtype = _type[j];
          const int jkind = _pair->getParticleKind(jtype);
          const double material_at_j = (_material == NULL ? 1.0 : _material[j]);
        
          if (_filter != NULL && !_filter->yes(ikind, jkind))
            continue;
        
          // compute distance between partile i and j
          double rsq = 0.0, rij[3] = {};
          for (int k=0;k<_dim;++k) {
            rij[k] = _x[i][k] - _x[j][k];
            rsq += (rij[k]*rij[k]);
          }
        
          const double cutsq = LOOKUP(_pair->cutsq, itype, jtype); 
          if (rsq < cutsq) {
            _mirror->setMirrorParticle(j);
            double coeff = 1.0;
            if (jkind & PairIsph::Solid)
              coeff=_mirror->computeMirrorCoefficient(sqrt(cutsq));
            
            const double r = sqrt(rsq) + ISPH_EPSILON;
            const double dwdr = _pair->kernel->dval(r, LOOKUP(_pair->h, itype, jtype));
          
            // 1st order consistent corrected gradient
            // \grad_{1} psi_{ij}^{k} = (psi_{j} - psi_{i}) G_i^{kq} dwdr (r_{ij}^{q} / r) V_{j}
            //                        = (psi_{j} - psi_{i}) ( G_i^{kq} r_{ij}^{q} ) dwdr / r V_{j} 
            const double vjtmp = dwdr/r*_vfrac[j]*coeff;          
            for (int k2=0;k2<_dim;++k2) {   // k
              double gitmp = 0.0;
              for (int k1=0;k1<_dim;++k1)   // q
                gitmp += VIEW2(G, _dim, k1, k2)*rij[k1];

              const double ijtmp = gitmp*vjtmp;

              if (_filter->yes(jkind))
                grad_material_at_i[k2] += ijtmp*(material_at_j - material_at_i);

              if (_scalar_diff)
                grad_at_i[k2] += ijtmp*(diff(j, i));
              
              if (_vector_diff)
                for (int k=0;k<_dim;++k)
                  grad3_at_i[k][k2] += ijtmp*(diff(j, i, k));
            }
          }
        }
      }

      {    
        // initialize the correction matrix
        const double *L = _pair->Lc[i];
      
        // compute laplacian 
        for (int jj=0;jj<jnum;++jj) {
          // j th particle; always put neighmask two bits are reserved for LAMMPS
          const int j = (jlist[jj] & NEIGHMASK);
          const int jtype = _type[j];
        
          // compute distance between partile i and j
          double rsq = 0.0, rij[3] = {};
          for (int k=0;k<_dim;++k) {
            rij[k] = _x[i][k] - _x[j][k];
            rsq += (rij[k]*rij[k]);
          }
        
          const double cutsq = LOOKUP(_pair->cutsq, itype, jtype);      
          if (rsq < cutsq) {
            const double r = sqrt(rsq) + ISPH_EPSILON;
            const double dwdr = _pair->kernel->dval(r, LOOKUP(_pair->h, itype, jtype));
          
            // compute unit vector e_{ij}
            double eij[3];
            for (int k=0;k<_dim;++k)
              eij[k] = rij[k]/r;
          
            // compute an intermediate scalar: aij
            // a_{ij} = 2.0 L \dot (e_{ij} \tprod dwdr e_{ij}) V_{j}
            //        = 2.0 L_{i}^{op} e_{ij}^{o} e_{ij}^{p} dwdr  V_{j}
            // if o == p, a_{ij}
            // if o != p, a_{ij} *= 2.0
            double aij = 0.0;
            const double scale[2] = {2.0, 1.0};
            for (int k2=0, op=0;k2<_dim;++k2)           // p
              for (int k1=0;k1<(k2+1);++k1,++op)        // o
                aij += L[op]*eij[k1]*eij[k2]*scale[k1==k2];
            aij *= 2.0*dwdr*_vfrac[j];
            
            if (_scalar_diff) {
              double bij = -1.0*util.dotVectors(_dim, eij, grad_at_i);
              _laplace_at_i += aij*bij;
            }
            
            if (_vector_diff)
              for (int k=0;k<_dim;++k) {
                double bij = -1.0*util.dotVectors(_dim, eij, grad3_at_i[k]);
                _laplace3_at_i[k] += aij*bij;
              }

            // store aij
            _aij[jj] = aij;
          }
        }
      }

      {    
        // compute laplacian 
        for (int jj=0;jj<jnum;++jj) {
          // j th particle; always put neighmask two bits are reserved for LAMMPS
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
        
          const double cutsq = LOOKUP(_pair->cutsq, itype, jtype);      
          if (rsq < cutsq) {
            _mirror->setMirrorParticle(j);
            double coeff = 1.0;
            if (jkind & PairIsph::Solid)
              coeff=_mirror->computeMirrorCoefficient(sqrt(cutsq));
          
            const double r = sqrt(rsq) + ISPH_EPSILON;
            //const double dwdr = _pair->kernel->dval(r, LOOKUP(_pair->h, itype, jtype));
            const double aij = _aij[jj];

            // compute an intermediate scalar: bij
            // b_{ij} = ( (psi_{i} - psi_{j}) / r - e_{ij} \dot \grad_{1} psi_{ij} ) 
            if (_scalar_diff) {
              double bij = coeff*(diff(i, j))/r;
              
              // compute f_{i} = \laplacian_{ij} = a_{ij} b_{ij}
              _laplace_at_i += aij*bij;
            }

            if (_vector_diff) {
              for (int k=0;k<_dim;++k) {
                // compute an intermediate scalar: bij
                // b_{ij} = ( (psi_{i} - psi_{j}) / r - e_{ij} \dot \grad_{1} psi_{ij} ) 
                double bij = coeff*(diff(i, j, k))/r;
                
                // compute f_{i} = \laplacian_{ij} = a_{ij} b_{ij}
                _laplace3_at_i[k] += aij*bij;
              }
            }
          }
        }
      }

      if (_scalar_diff) 
        _laplace_at_i = _alpha*(material_at_i*_laplace_at_i + 
                                util.dotVectors(_dim, grad_material_at_i, grad_at_i));
      
      if (_vector_diff) 
        for (int k=0;k<_dim;++k) 
          _laplace3_at_i[k] = _alpha*(material_at_i*_laplace3_at_i[k] + 
                                      util.dotVectors(_dim, grad_material_at_i, grad3_at_i[k]));
    }

    // compute laplacian for given varialbe f
    // =============================================================================
    template<class PairIsph>
    class FunctorOuterLaplacian : public FunctorOuterLaplacianHelper<PairIsph> {
    public:
      // laplacian for scalar field
      FunctorOuterLaplacian(PairIsph *isph, 
                            double *f, 
                            double alpha = 1.0,
                            double *material = NULL)
        : FunctorOuterLaplacianHelper<PairIsph>(isph, alpha, material, true, false),
          _f(f),
          _f3(NULL) { }

      // component wise laplacian
      FunctorOuterLaplacian(PairIsph *isph, 
                            double **f3, 
                            double alpha = 1.0,
                            double *material = NULL)
        : FunctorOuterLaplacianHelper<PairIsph>(isph, alpha, material, false, true),
          _f(NULL),
          _f3(f3) { }

      double diff(const int i, const int j);
      double diff(const int i, const int j, const int k);

    protected:
      FUNCTOR_REMOVE_THIS_POINTER;
      double *_f,**_f3;
    };

    template<class PairIsph> inline double
    FunctorOuterLaplacian<PairIsph>::diff(const int i, const int j) {
      return (_f[i] - _f[j]);
    }      

    template<class PairIsph> inline double
    FunctorOuterLaplacian<PairIsph>::diff(const int i, const int j, const int k) {
      return (_f3[i][k] - _f3[j][k]);
    }      

  }
}

#endif

