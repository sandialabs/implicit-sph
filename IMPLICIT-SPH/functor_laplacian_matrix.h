#pragma once
#ifndef __FUNCTOR_LAPLACIAN_MATRIX_H__
#define __FUNCTOR_LAPLACIAN_MATRIX_H__

#include <iostream>
#include <string>
#include "utils.h"
#include "mirror.h"
#include "functor.h"

namespace LAMMPS_NS {

  namespace Corrected {

    template<class PairIsph, bool AntiSymmetric = false>
    class FunctorOuterLaplacianMatrix : public FunctorOuter<PairIsph> {
    public:
      FunctorOuterLaplacianMatrix(PairIsph *isph,
                                  double alpha,
                                  double *material = NULL,
                                  int iblock = -1,
                                  double** normal = NULL)
        : FunctorOuter<PairIsph>(isph),_mirror(NULL),
          _alpha(alpha),
          _material(material),
          _normal(normal),
          _iblock(iblock) {
        if (_iblock < 0)
          _A = this->_pair->A.crs;
        else
          _A = this->_pair->A_blk(_iblock,_iblock);
      }

      virtual void enterFor();
      virtual void operator()(const int ii);
      virtual void exitFor();

      virtual void createMirrorFunction();

    protected:
      Epetra_CrsMatrix * _A;
      Epetra_IntSerialDenseVector _idx, _lidx;
      Epetra_SerialDenseVector _val, _valmn, _Qi, _Qj;
      MirrorFunction *_mirror;
      double _alpha, *_material, **_normal;
      int _iblock;

      FUNCTOR_REMOVE_THIS_POINTER;
    };

    template<class PairIsph, bool AntiSymmetric> inline void
    FunctorOuterLaplacianMatrix<PairIsph,AntiSymmetric>::createMirrorFunction() {
      _mirror = new MirrorNothing(_dim, _type, _pair->h);
    }

    template<class PairIsph, bool AntiSymmetric> inline void
    FunctorOuterLaplacianMatrix<PairIsph,AntiSymmetric>::enterFor() {
      if (_filter == NULL)
        _pair->error->all(FLERR, "FunctorLaplacianMatrix:: Filter is required to consider dirichlet-like boundary");

      int n = _pair->tags_in_cut->MaxNumIndices();
      _idx.Size(n);
      _val.Size(n);
      _valmn.Size(n);
      _lidx.Size(n);
      _Qi.Size(_dim*_dim);
      _Qj.Size(_dim*_dim);

      createMirrorFunction();
    }

    template<class PairIsph, bool AntiSymmetric> inline void
    FunctorOuterLaplacianMatrix<PairIsph,AntiSymmetric>::operator()(const int ii) {
      const int i = _ilist[ii];
      const int itype = _type[i];
      const int ikind = _pair->getParticleKind(itype);
      const double material_at_i = (_material == NULL ? 1.0 : _material[i]);

      // extract pre-allocated row
      double *val = _val.Values();
      double *valmn = _valmn.Values();

      int *idx = _idx.Values();
      int *lidx = _lidx.Values();

      // for solid particles, diagonal terms are unit
      // scaling is still applied to them (alpha and material)
      if (!_filter->yes(ikind)) {
        idx[0] = _tag[i];

        // when we do not compute for i, leave the diagonal entree zero
        // and later need to change diagonal
        val[0] = 0.0; //_alpha;
        _A->ReplaceGlobalValues(_tag[i], 1, val, idx);
        return;
      }

      const bool has_normal = (_normal != NULL) && (util.dotVectors(_dim,_normal[i],_normal[i]) > 0.5);
//      if (_normal != NULL) { //return if we are not near the boundary and the interaction of particle i with particles of the same type.
//        if (!has_normal && !_filter->yes(ikind, ikind))
//          return;
//      }

      // neighbor around i
      const int *jlist = _firstneigh[i];
      const int jnum = _numneigh[i];

      double grad_material_at_i[3] = {};

      // c_{i}^{k} = \sum_{j} aij e_{ij}^{k}
      double ci[3] = {};

      // mirroring for particle i
      _mirror->setFreeParticle(i);

      int cnt;
      {
        // initialize the correction matrix
        const double *L = (AntiSymmetric ? _pair->Li : _pair->Lc[i]);
        const double *G = (AntiSymmetric ? _pair->Gi : _pair->Gc[i]);

        double diag = 0.0;
        cnt = 0;

        // compute laplacian

        for (int jj=0;jj<jnum;++jj) {
          // j th particle; always put neighmask two bits are reserved for LAMMPS
          const int j = (jlist[jj] & NEIGHMASK);
          const int jtype = _type[j];
          const int jkind = _pair->getParticleKind(jtype);
          const double material_at_j = (_material == NULL ? 1.0 : _material[j]);

          // compute distance between partile i and j
          double rsq = 0.0, rij[3] = {};
          for (int k=0;k<_dim;++k) {
            rij[k] = _x[i][k] - _x[j][k];
            rsq += (rij[k]*rij[k]);
          }

          const double cutsq = LOOKUP(_pair->cutsq, itype, jtype);
          if (rsq < cutsq) {
            _mirror->setMirrorParticle(j);
            double coeff = _filter->yes(ikind,ikind);
            if (!(ikind & PairIsph::Solid) && (jkind & PairIsph::Solid)) 
              coeff = (_filter->yes(ikind,jkind) ? _mirror->computeMirrorCoefficient(sqrt(cutsq)) : 0.0);

            const double r = sqrt(rsq) + ISPH_EPSILON;
            const double dwdr = _pair->kernel->dval(r, LOOKUP(_pair->h, itype, jtype));

            // compute unit vector e_{ij}
            double eij[3];
            for (int k=0;k<_dim;++k)
              eij[k] = rij[k]/r;

            const double vfrac = (AntiSymmetric ? sqrt(_vfrac[i]*_vfrac[j]) : _vfrac[j]);
            const double vjtmp = dwdr*vfrac;
            for (int k2=0;k2<_dim;++k2) {   // k
              double gitmp = 0.0;
              for (int k1=0;k1<_dim;++k1)   // q
                gitmp += VIEW2(G, _dim, k1, k2)*eij[k1];

              double ijtmp = gitmp*vjtmp;
              if (ikind & jkind)
                grad_material_at_i[k2] += ijtmp*(sphOperator<AntiSymmetric>(material_at_i, material_at_j));
            }

            // a_{ij} = 2.0 L \dot (e_{ij} \tprod dwdr e_{ij}) V_{j}
            //        = 2.0 L_{i}^{op} e_{ij}^{o} e_{ij}^{p} dwdr V_{j}
            // if o == p, a_{ij}
            // if o != p, a_{ij} *= 2.0
            double aij = 0.0;
            const double scale_a[2] = {2.0, 1.0};
            for (int k2=0, op=0;k2<_dim;++k2)    // p
              for (int k1=0;k1<(k2+1);++k1,++op)       // o
                aij += L[op]*eij[k1]*eij[k2]*scale_a[k1==k2];
            aij *= 2.0*dwdr*vfrac;

            // precompute ci
            if (!AntiSymmetric) {
              for (int k=0;k<_dim;++k)
                ci[k] += aij*eij[k];
            }

            // redefine aij
            aij *= material_at_i*coeff/r;

            val[cnt] = -aij;
            diag    +=  aij;

            idx[cnt] = _tag[j];
            lidx[cnt] = j;
            ++cnt;
          }
        }
        val[cnt] = diag;
        idx[cnt] = _tag[i];
        lidx[cnt] = i;

        ++cnt;
      }

      {
        // initialize the correction matrix
        const double *G = (AntiSymmetric ? _pair->Gi : _pair->Gc[i]);
        double diag = 0.0;
        cnt = 0;

        // compute gradient

        for (int jj=0;jj<jnum;++jj) {
          const int j = (jlist[jj] & NEIGHMASK);
          const int jtype = _type[j];
          const int jkind = _pair->getParticleKind(jtype);

          // compute distance between partile i and j
          double rsq = 0.0, rij[3] = {};
          for (int k=0;k<_dim;++k) {
            rij[k] = _x[i][k] - _x[j][k];
            rsq += (rij[k]*rij[k]);
          }

          const double cutsq = LOOKUP(_pair->cutsq, itype, jtype);
          if (rsq < cutsq) {
            double coeff = _filter->yes(ikind,ikind);
            if (!(ikind & PairIsph::Solid) && (jkind & PairIsph::Solid)) 
              coeff = _filter->yes(ikind,jkind);

            const double r = sqrt(rsq) + ISPH_EPSILON;
            const double dwdr = _pair->kernel->dval(r, LOOKUP(_pair->h, itype, jtype));

            const double vfrac = (AntiSymmetric ? sqrt(_vfrac[i]*_vfrac[j]) : _vfrac[j]);
            const double vjtmp = dwdr*vfrac;

            // unit vector
            double eij[3];
            for (int k=0;k<_dim;++k)
              eij[k] = rij[k]/r;

            // 1st order consistent corrected gradient
            // a_{ij} = e_{ij} \dot ( G_{i} e_{ij} ) dwdr V_{j}
            //        = e_{ij}^{k} \dot ( G_{i}^{kq} e_{ij}^{q} ) dwdr V_{j}
            double bij[3] = {};
            for (int k2=0;k2<_dim;++k2)    // k
              for (int k1=0;k1<_dim;++k1)   // q
                bij[k2] += VIEW2(G, _dim, k1, k2)*eij[k1];

            const double tmp = coeff*(material_at_i*util.dotVectors(_dim, bij, ci)*vjtmp -
                                      util.dotVectors(_dim, bij, grad_material_at_i)*vjtmp);

            val[cnt] -= tmp;
            diag     += tmp;

            if (idx[cnt] != _tag[j])
              _pair->error->all(FLERR, "FunctorLaplacianMatrix:: Index does not match");

            ++cnt;
          }
        }
        val[cnt] += diag;
        idx[cnt]  = _tag[i];
        lidx[cnt] = i;
        ++cnt;
      }

      {
        _val.Scale(_alpha);

        if (_iblock < 0 || !has_normal )
          _A->SumIntoGlobalValues(_tag[i], cnt, val, idx);
        else  {
          double n_i[3] = {};/// = _normal[i];
          for(int k=0; k<cnt; ++k) {
            const int j = lidx[k];
            for(int d=0; d<_dim; ++d)
              n_i[d] += _normal[j][d];
          }
          double norm = sqrt(util.dotVectors(_dim, n_i, n_i));
          for(int k=0; k<_dim; ++k)
            n_i[k] /= norm;


   //*
       //   int ib[3];
       int ib=0;
        //  util.computeRotationMatrix(_dim, n_i, ib, _Qi.Values());
          for(; ib < _dim-1 && (n_i[ib]*n_i[ib] < 1.0/_dim); ++ib);  //pick component ib of normal which is big enough
     //     for(int ib=0; ib < _dim; ++ib) {
            for(int jb=0; jb<_dim; ++jb) {
              for(int k=0; k<cnt; ++k) {
                const int j = lidx[k];
                valmn[k] = val[k]*n_i[jb]*n_i[ib];///fabs(n_i[ib[0]]);
              }
              this->_pair->A_blk(ib,jb)->SumIntoGlobalValues(_tag[i], cnt, valmn, idx);
            }
       //   }
/*/
          util.computeRotationMatrix(_dim, n_i, _Qi.Values());
          for(int m=0; m<_dim; ++m) {  //Not efficient
            for(int n=0; n<_dim; ++n) {
              for(int k=0; k<cnt; ++k) {
                const int j = lidx[k];
                const double *n_j = _normal[j];
                util.computeRotationMatrix(_dim, n_j, _Qj.Values());
                if(util.dotVectors(_dim,n_j,n_j)<0.5)
                  valmn[k] = val[k]*(m==_iblock)*(n==_iblock);
                else 
                  valmn[k] = val[k]*VIEW2(_Qi, _dim, _iblock, m)*VIEW2(_Qj, _dim, _iblock, n);
              }
              this->_pair->A_blk(m,n)->SumIntoGlobalValues(_tag[i], cnt, valmn, idx);
            }
          }
         // */
        }
      }
    }

    template<class PairIsph, bool AntiSymmetric> inline
    void FunctorOuterLaplacianMatrix<PairIsph,AntiSymmetric>::exitFor() {

      // finalize connectivity construction
      if (_pair->A.is_filled == 0) {
        _pair->A.crs->FillComplete();
        _pair->A.crs->OptimizeStorage();
        _pair->A.is_filled = 1;
      }
      delete _mirror;
    }

    template<class PairIsph> using FunctorOuterLaplacianMatrixSymmetric     = FunctorOuterLaplacianMatrix<PairIsph,false>;
    template<class PairIsph> using FunctorOuterLaplacianMatrixAntiSymmetric = FunctorOuterLaplacianMatrix<PairIsph,true>;
  }
}

#endif

