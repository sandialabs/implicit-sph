#pragma once
#ifndef __FUNCTOR_BOUNDARY_DIRICHLET_H__
#define __FUNCTOR_BOUNDARY_DIRICHLET_H__

#include <iostream>
#include <string>
#include <vector>
#include "utils.h"
#include "functor.h"

namespace LAMMPS_NS {

  namespace Corrected {

    template<class PairIsph>
    class FunctorOuterBoundaryDirichlet : public FunctorOuter<PairIsph> {
    public:
      FunctorOuterBoundaryDirichlet(PairIsph *isph,
                                    double **n, double *b, int lda, double scale=1)
        : FunctorOuter<PairIsph>(isph), _n(n), _b(b), _lda(lda), _scale(scale) { }

      void enterFor();
      void operator()(const int ii);

    protected:
      Epetra_IntSerialDenseVector _idx;
      Epetra_SerialDenseVector _val;
      double **_n, *_b;
      int _lda;
      double _scale;

      FUNCTOR_REMOVE_THIS_POINTER;
    };

    template<class PairIsph> inline void
    FunctorOuterBoundaryDirichlet<PairIsph>::enterFor() {
      if (_filter != NULL)
        _pair->error->all(FLERR, "FunctorBoundaryDirichlet:: Filter is not allowed");

      int n = _pair->tags_in_cut->MaxNumIndices();
      _idx.Size(n);
      _val.Size(n);
    }

    template<class PairIsph> inline void
    FunctorOuterBoundaryDirichlet<PairIsph>::operator()(const int ii) {
      int i = _ilist[ii];
      int itype = _type[i];
      int ikind = _pair->getParticleKind(itype);

      // neighbor around i
      int *jlist = _firstneigh[i];
      int jnum = _numneigh[i];

      // extract pre-allocated row
      double *val = _val.Values();
      int *idx = _idx.Values();

      switch (ikind) {
      case PairIsph::Solid: {
        // do nothing;
        break;
      }
      case PairIsph::Fluid: {
        int cnt = 0, n_solid = 0;
        for (int jj=0;jj<jnum;++jj) {
          int j = (jlist[jj] & NEIGHMASK);
          int jtype = _type[j];
          double rsq = 0.0, rij[3] = {};
          for (int k=0;k<_dim;++k) {
            rij[k] = _x[i][k] - _x[j][k];
            rsq += (rij[k]*rij[k]);
          }

          if (rsq < pow(LOOKUP(_pair->h, itype, jtype),2))
            n_solid += int(_pair->getParticleKind(jtype) == PairIsph::Solid);
        }


        double xn_i(0), xn_av(0);

        if(n_solid > 0) {
          //compute normal coordinate of boundary
          map<double, int> normal_coord;

          for (int k=0;k<_dim;++k)
            xn_i += _x[i][k] * _n[i][k];

          vector<double> xn(jnum), tag(jnum);

          int natoms_cut = 0;
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

            // we use h^2 instead of cut to precisely define the interface layer
            if (rsq < pow(LOOKUP(_pair->h, itype, jtype),2)) {

              double xn_j(0.);
              for (int k=0;k<_dim;++k)
                xn_j += _x[j][k]*_n[i][k];

              xn_av += xn_j;
              tag[natoms_cut] = _tag[j];
              xn[natoms_cut++] = xn_j;
            }
            else if(rsq < LOOKUP(_pair->cutsq, itype, jtype)) {
              val[cnt]=0;
              idx[cnt++]=_tag[j];
            }
          }
          xn_av /= natoms_cut;

          double xterm(0.);
          for(int j=0; j<natoms_cut; j++)
            xterm += xn[j] * ( xn[j] - xn_av );

          for(int j=0; j<natoms_cut; j++) {
            val[cnt] = -(xn_i - xn_av)*(xn[j]-xn_av)/xterm - 1.0/natoms_cut;//
            idx[cnt++]  = tag[j];
          }
          val[cnt] = 1.0;
          idx[cnt]  = _tag[i];

          _pair->A.crs->ReplaceGlobalValues(_tag[i], ++cnt, val, idx);

          for (int k=0;k<_dim;++k)
            VIEW2(_b, _lda, i, k) = 0;
        }

        break;
      }
      default:
        _pair->error->all(FLERR, "FunctorBoundaryDirichlet:: Particle types are not supported"); 
      }
    }

  }
}
#endif

