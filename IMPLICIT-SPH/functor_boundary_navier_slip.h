#pragma once
#ifndef __FUNCTOR_BOUNDARY_NAVIER_SLIP_H__
#define __FUNCTOR_BOUNDARY_NAVIER_SLIP_H__

#include <iostream>
#include <string>
#include "utils.h"
#include "functor.h"

namespace LAMMPS_NS {

  namespace Corrected {
    // implement robin boundary condition
    // if add_neumann_term == true,  - \nabla u n = beta u - u_s
    // if add_neumann_term == false, - \nabla u n = beta u
    template<class PairIsph>
    class FunctorOuterBoundaryNavierSlip : public FunctorOuter<PairIsph> {
    public:
      FunctorOuterBoundaryNavierSlip(PairIsph *isph, 
                                     double beta, double **n, 
                                     double *rho,
                                     bool add_neumann_term = true,
                                     int iblock = -1)
        : FunctorOuter<PairIsph>(isph),
          _beta(beta), _n(n), _rho(rho), _add_neumann_term(add_neumann_term), _iblock(iblock) {}

      void enterFor();
      void operator()(const int ii);

    protected:
      Epetra_IntSerialDenseVector _idx, _lidx;
      Epetra_SerialDenseVector _val, _valmn, _Qi, _Qj;
      double _beta,**_n,*_rho;
      int _add_neumann_term, _iblock;

      FUNCTOR_REMOVE_THIS_POINTER;
    };

    template<class PairIsph> inline void
    FunctorOuterBoundaryNavierSlip<PairIsph>::enterFor() {
      if (_filter != NULL)
        _pair->error->all(FLERR, "FunctorBoundaryNavierSlip:: Filter is not allowed");

      int n = _pair->tags_in_cut->MaxNumIndices();
      _idx.Size(n);
      _lidx.Size(n);
      _val.Size(n);
      _valmn.Size(n);
      _Qi.Size(_dim*_dim);
      _Qj.Size(_dim*_dim);
    }

    template<class PairIsph> inline void
    FunctorOuterBoundaryNavierSlip<PairIsph>::operator()(const int ii) {
      int i = _ilist[ii];
      int itype = _type[i];
      int ikind = _pair->getParticleKind(itype);

      switch (ikind) {
      case PairIsph::Solid: {
        // do nothing;
        break;
      }
      case PairIsph::BufferDirichlet:
      case PairIsph::BufferNeumann:
      case PairIsph::Fluid: {
        // neighbor around i
        int *jlist = _firstneigh[i];
        int jnum = _numneigh[i];

        // extract pre-allocated row
        double *val = _val.Values();
        double *valmn = _valmn.Values();
        int *idx = _idx.Values();
        int *lidx = _lidx.Values();

        // if solid particles are stationary but has non-zero velocity
        double robin_at_i = 0.0;
        int cnt = 0;

        // compute robin boundary
        double *G = _pair->Gc[i];
        for (int jj=0;jj<jnum;++jj) {
          // j th particle; always put neighmask two bits are reserved for LAMMPS
          int j = (jlist[jj] & NEIGHMASK);
          int jtype = _type[j];

          // compute for only solid particles
          if (_pair->getParticleKind(jtype) == PairIsph::Solid) {
            // compute distance between partile i and j
            double rsq = 0.0, rij[3] = {};
            for (int k=0;k<_dim;++k) {
              rij[k] = _x[i][k] - _x[j][k];
              rsq += (rij[k]*rij[k]);
            }
            if (rsq < LOOKUP(_pair->cutsq, itype, jtype)) {
              double r = sqrt(rsq) + ISPH_EPSILON;
              double dwdr = _pair->kernel->dval(r, LOOKUP(_pair->h, itype, jtype));
            
              // a_{ij} = (n_{ij} + n_{ji}) \dot ( G_{i} e_{ij} ) dwdr V_{j}
              //        = (n_{ij}^{k} +n_{ji}^{k}) \dot ( G_{i}^{kq} r_{ij}^{q} ) dwdr/r V_{j}
              double aij[3] = {};
              for (int k2=0;k2<_dim;++k2)     // k
                for (int k1=0;k1<_dim;++k1)   // q
                  aij[k2] += VIEW2(G, _dim, k1, k2)*rij[k1];
            
              double tmp = 0.0;
              for (int k=0;k<_dim;++k)
                tmp += (_n[i][k] + _n[j][k])*aij[k];

              double beta_coeff = _beta;

//              double w = 1.0;
//              beta_coeff = _beta*(1.0+ sin(w*_x[i][0]) * sin(w*_x[i][1]));

              double robin_at_j = beta_coeff*dwdr/r*_vfrac[j]/_rho[i]*tmp;

              if(_add_neumann_term) {
                val[cnt] =  robin_at_j;
                idx[cnt] = _tag[j];
                lidx[cnt] = j;
                ++cnt;
              }

              robin_at_i -= robin_at_j;
            }
          }
        }

        val[cnt] = robin_at_i;
        idx[cnt] = _tag[i];
        lidx[cnt] = i;
        ++cnt;

        if(_iblock >= 0) {
//*

          double n_i[3] ={};/// = _normal[i];
          for(int k=0; k<cnt; ++k) {
            const int j = lidx[k];
            for(int d=0; d<_dim; ++d)
              n_i[d] += _n[j][d];
          }
          double norm = sqrt(util.dotVectors(_dim, n_i, n_i));
          if(norm !=0 )
            for(int k=0; k<_dim; ++k)
              n_i[k] /= norm;


       //   for(; iib < _dim-1 && (_n[i][iib]*_n[i][iib] < 1.0/_dim); ++iib);
          for(int ib=0; ib < _dim; ++ib) {
            for(int jb=0; jb <_dim; ++jb) {
              for(int k=0; k<cnt; ++k) {
                const int j = lidx[k];
                valmn[k] = val[k]*(double(ib==jb)-n_i[jb]*n_i[ib]);///fabs(n_i[ib[0]]);
              }
              this->_pair->A_blk(ib,jb)->SumIntoGlobalValues(_tag[i], cnt, valmn, idx);
            }
          }
          /*/
          int ib[3];
          util.computeRotationMatrix(_dim, _n[i], ib, _Qi.Values());

          for(int m=0; m<_dim; ++m) {  //Not efficient
            for(int n=0; n<_dim; ++n) {
              for(int k=0; k<cnt; ++k) {
                int j = lidx[k];
                util.computeRotationMatrix(_dim, _n[j],ib, _Qj.Values());
                valmn[k] = val[k] * VIEW2(_Qi, _dim, _iblock,  m) * VIEW2(_Qj, _dim, _iblock, n); //Q' B Q
              }
              this->_pair->A_blk(m,n)->SumIntoGlobalValues(_tag[i], cnt, valmn, idx);
            }
          } //*/
        }
        else {
        // diagonal entry will separately updated
        _pair->A.crs->SumIntoGlobalValues(_tag[i], cnt, val, idx);
        }
        break;
      }
      default:
        _pair->error->all(FLERR, "FunctorBoundaryNavierSlip:: Particle types are not supported"); 
      }
    }

  }
}

#endif

