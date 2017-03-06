#pragma once
#ifndef __FUNCTOR_COMPUTE_NORMAL_H__
#define __FUNCTOR_COMPUTE_NORMAL_H__

#include <iostream>
#include <string>
#include <map>
#include "math.h"
#include "float.h"
#include "utils.h"
#include "functor.h"

namespace LAMMPS_NS {

  namespace Corrected { 

    template<class PairIsph>
    class FunctorOuterNormal: public FunctorOuter<PairIsph> {
    public:
      FunctorOuterNormal(PairIsph *isph, 
                         int *part, double *orient, 
                         double **normal,
                         double *pnd, double* bd_coord) 
        : FunctorOuter<PairIsph>(isph), 
          _part(part),
          _orient(orient), _normal(normal), 
          _pnd(pnd), _bd_coord(bd_coord) { }

      void enterFor();
      void operator()(const int ii);

      double* getNormal();
      double getParticleNumberDensity();

    protected:
      int *_part;
      double *_orient, **_normal, *_pnd, *_bd_coord;
      double _normal_at_i[3], _pnd_at_i;

      FUNCTOR_REMOVE_THIS_POINTER;
    };

    template<class PairIsph> inline double*
    FunctorOuterNormal<PairIsph>::getNormal() {
      return _normal_at_i;
    }
    template<class PairIsph> inline double 
    FunctorOuterNormal<PairIsph>::getParticleNumberDensity() {
      return _pnd_at_i;
    }
    template<class PairIsph> inline void 
    FunctorOuterNormal<PairIsph>::enterFor() {
      if (_filter == NULL)
        _pair->error->all(FLERR, "FunctorNormal:: Filter is required which defines interface of different particles");
    }

    template<class PairIsph> inline void 
    FunctorOuterNormal<PairIsph>::operator()(const int ii) {
      const int i = _ilist[ii];
      const int itype = _type[i];
      const int ikind = _pair->getParticleKind(itype);

      const bool use_part = (_part != NULL);
      const int ipart = (use_part ? _part[i] : ikind);

      // neighbor around i
      int *jlist = _firstneigh[i];
      int jnum = _numneigh[i];

      // normal is computed on the given storage or member variable
      memset(_normal_at_i, 0, sizeof(double) * _dim);
      _pnd_at_i = 0.0;

      if (!_filter->yes(ikind))
        return;

      // initialize the correction matrix
      double *G = _pair->Gc[i];

      for (int jj = 0; jj < jnum; ++jj) {
        const int j = (jlist[jj] & NEIGHMASK);
        const int jtype = _type[j];
        const int jkind = _pair->getParticleKind(jtype);
        const int jpart = (use_part ? _part[j] : jkind);

        // compute distance between partile i and j
        double rsq = 0.0, rij[3] = { };
        for (int k = 0; k < _dim; ++k) {
          rij[k] = _x[i][k] - _x[j][k];
          rsq += (rij[k] * rij[k]);
        }
        const double r = sqrt(rsq) + ISPH_EPSILON;

        if (rsq < LOOKUP(_pair->cutsq, itype, jtype)) {
          if (_filter->yes(ipart, jpart)) {
            const double orient = _orient[ikind];
            const double dwdr = _pair->kernel->dval(r, LOOKUP(_pair->h, itype, jtype));

            // 1st order consistent corrected normal
            // \normal_{1} psi_{ij}^{k} =  orient*G_i^{kq} dwdr (r_{ij}^{q} / r V_{j}
            //                          =  orient*( G_i^{kq} r_{ij}^{q} ) dwdr / r V_{j}
            for (int k2=0;k2<_dim;++k2) {   // k
              double gitmp = 0.0;
              for (int k1=0;k1<_dim;++k1)// q
                gitmp += VIEW2(G, _dim, k1, k2)*rij[k1];
              _normal_at_i[k2] += gitmp*orient*dwdr/r*_vfrac[j];
            }
          } else {
            _pnd_at_i += _pair->kernel->val(r, LOOKUP(_pair->h, itype, jtype));
          }
        }
      }

      // account for the self contribution
      _pnd_at_i += _pair->kernel->val(0.0, LOOKUP(_pair->h, itype, itype));

      // scale with alpha
      double alpha = 0.0;
      for (int k = 0; k < _dim; ++k)
        alpha += _normal_at_i[k] * _normal_at_i[k];

      alpha = sqrt(alpha);

      if (alpha != 0.0)
        for (int k = 0; k < _dim; ++k)
          _normal_at_i[k] /= alpha;

      // store normal at i if normal is not null
      if (_normal != NULL)
        memcpy(_normal[i], _normal_at_i, sizeof(double) * _dim);

      if (_pnd != NULL)
        _pnd[i] = _pnd_at_i;

      // ------------------------------------------------------------------------


      if (_bd_coord != NULL) {
        int n_solid(0);
        for (int jj = 0; jj < jnum; ++jj) {
          int j = (jlist[jj] & NEIGHMASK);
          int jtype = _type[j];
          n_solid += int(_pair->getParticleKind(jtype) == PairIsph::Solid);
        }

        if (n_solid > 0) {
          //compute normal coordinate of boundary
          map<double, int> normal_coord;
          int n_fluid(0);
        
          //_n[i][0]=0; _n[i][1] = (_x[i][1]<0) - (_x[i][1]>0);
          for (int jj = 0; jj < jnum; ++jj) {
            // j th particle; always put neighmask two bits are reserved for LAMMPS
            int j = (jlist[jj] & NEIGHMASK);
            int jtype = _type[j];
          
            // compute distance between partile i and j
            double rsq = 0.0, rij[3] = { };
            for (int k = 0; k < _dim; ++k) {
              rij[k] = _x[i][k] - _x[j][k];
              rsq += (rij[k] * rij[k]);
            }
            if (rsq < LOOKUP(_pair->cutsq, itype, jtype)) {
              double ncoord(0);
              for (int k=0;k<_dim;++k)
                ncoord += _x[j][k]*_normal_at_i[k];
            
              normal_coord.insert(make_pair(ncoord, jtype));
            }
          }
        
          //self contribution
          double ncoord(0), bd_ncoord(0);
          for (int k = 0; k < _dim; ++k)
            ncoord += _x[i][k] * _normal_at_i[k];
          normal_coord.insert(make_pair(ncoord, itype));
        
          int misclassification = max(n_solid, n_fluid);
          for (auto it=normal_coord.begin();it!=normal_coord.end();++it) {
            n_solid -= int(_pair->getParticleKind(it->second) == PairIsph::Solid);
            n_fluid += int(_pair->getParticleKind(it->second) == PairIsph::Fluid);
            int misclass = max(n_solid, n_fluid);
            if (misclass > misclassification) {
              _bd_coord[ii] = 0.5 * (bd_ncoord + it->first);
              break;
            }
          
            misclassification = misclass;
            bd_ncoord = (it)->first;
          }
        }
      }
    }

  }
}
#endif

