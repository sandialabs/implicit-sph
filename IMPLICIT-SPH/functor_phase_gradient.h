#pragma once
#ifndef __FUNCTOR_PHASE_GRADINET_H__
#define __FUNCTOR_PHASE_GRADIENT_H__

#include <iostream>
#include <string>

#include "utils.h"
#include "functor.h"
#include "color.h"

namespace LAMMPS_NS {

  namespace Corrected { 

    template<class PairIsph>
    class FunctorOuterPhaseGradient: public FunctorOuter<PairIsph> {
    public:
      FunctorOuterPhaseGradient(PairIsph *isph, double *rho,
                                double **grad)
        : FunctorOuter<PairIsph>(isph),
        _rho(rho),
        _grad(grad) { }

      double* getGradient();
      
      void enterFor();
      void operator()(const int ii);

    protected:
      double *_rho, **_grad;
      double _grad_at_i[3];

      FUNCTOR_REMOVE_THIS_POINTER;
    };

    template<class PairIsph> inline double*
    FunctorOuterPhaseGradient<PairIsph>::getGradient() {
      return _grad_at_i;
    }
    
    template<class PairIsph> inline void 
    FunctorOuterPhaseGradient<PairIsph>::enterFor() {
      if (_filter == NULL)
        _pair->error->all(FLERR, "FunctorPhaseGradient:: Filter is required");
    }

    template<class PairIsph> inline void 
    FunctorOuterPhaseGradient<PairIsph>::operator()(const int ii) {
      const int i = _ilist[ii];
      const int itype = _type[i];

      const int ikind  = _pair->getParticleKind(itype);
      const int iphase = _pair->getParticlePhase(itype);

      const double irho = (_rho != NULL ? _rho[i] : 1.0);
      
      // neighbor around i
      const int *jlist = _firstneigh[i];
      const int jnum = _numneigh[i];

      // grad is computed on the given storage or member variable
      memset(_grad_at_i, 0, sizeof(double)*_dim);

      if (!_filter->yes(ikind))
        return;

      // initialize the correction matrix
      double *G = _pair->Gc[i];

      double vol_in_phase = _vfrac[i], vol_out_phase = 0.0;

      for (int jj=0;jj<jnum;++jj) {
        const int j = (jlist[jj] & NEIGHMASK);
        const int jtype = _type[j];

        const int jkind  = _pair->getParticleKind(jtype);
        const int jphase = _pair->getParticlePhase(jtype);

        const double jrho = (_rho != NULL ? _rho[j] : 1.0);

        bool ij_in_phase = true;
        if (_filter->yes(ikind, jkind)) {
          if (iphase != jphase) {
            // compute distance between partile i and j
            double rsq = 0.0, rij[3] = { };
            for (int k=0;k<_dim;++k) {
              rij[k] = _x[i][k] - _x[j][k];
              rsq += (rij[k] * rij[k]);
            }
            
            if (rsq < LOOKUP(_pair->cutsq, itype, jtype)) { // this should be replace by the csf model
              const double r = sqrt(rsq) + ISPH_EPSILON;              
              const double cij = _pair->st.csf.color->val(iphase, irho, jphase, jrho);
              const double dwdr = _pair->kernel->dval(r, LOOKUP(_pair->h, itype, jtype));

              ij_in_phase = false;
              vol_out_phase += _vfrac[j];

              // we do not use corrected gradients and use the standard average scheme
              switch (_pair->st.csf.color->gradtype()) {
              case ColorFunction::CorrectedGrad: {
                // ** 1st order consistent corrected grad
                for (int k2=0;k2<_dim;++k2) {   // k
                  double gitmp = 0.0;
                  for (int k1=0;k1<_dim;++k1)// q
                    gitmp += VIEW2(G, _dim, k1, k2)*rij[k1];
                  _grad_at_i[k2] += gitmp*cij*dwdr/r*_vfrac[j];
                }
                break;
              }
              case ColorFunction::AdamiGrad: {
                // ** standard gradient from Adami paper
                for (int k=0;k<_dim;++k) 
                  _grad_at_i[k] += (pow(_vfrac[i], 2) + pow(_vfrac[j], 2))*cij*dwdr*(rij[k]/r)/_vfrac[i];
                break;
              }
              }
            }
          }
        } 
        
        if (ij_in_phase)
          vol_in_phase += _vfrac[j]; 
      }

      // remove inaccurate normal
      // const double mag = sqrt(util.dotVectors(_dim, _grad_at_i, _grad_at_i));
      // const double eps = _pair->st.csf.epsilon*LOOKUP(_pair->h, itype, itype);
      // if (mag < eps)
      //   memset(_grad_at_i, 0, sizeof(double)*3);
      
      const double vol_ratio = vol_in_phase/(vol_in_phase + vol_out_phase);
      const double eps = _pair->st.csf.epsilon;

      if (vol_ratio < eps || vol_ratio > (1.0 - eps))
        memset(_grad_at_i, 0, sizeof(double)*3);

      if (_grad != NULL)
        memcpy(_grad[i], _grad_at_i, sizeof(double)*_dim);
    }

  }
}
#endif

