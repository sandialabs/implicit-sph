#pragma once
#ifndef __FUNCTOR_CORRECT_PHASE_NORMAL_H__
#define __FUNCTOR_CORRECT_PHASE_NORMAL_H__

#include <iostream>
#include <string>

#include "utils.h"
#include "functor.h"

namespace LAMMPS_NS {

  namespace Corrected {

    template<class PairIsph>
    class FunctorOuterCorrectPhaseNormal : public FunctorOuter<PairIsph> {
    public:
      FunctorOuterCorrectPhaseNormal(PairIsph *isph,
                                     double *pnd,
                                     double **knormal,
                                     double **pnormal)
        : FunctorOuter<PairIsph>(isph),
          _pnd(pnd),
          _knormal(knormal),
          _pnormal(pnormal) { }

      void enterFor();
      void operator()(const int ii);

    protected:
      double *_pnd, **_knormal, **_pnormal;

      FUNCTOR_REMOVE_THIS_POINTER;
    };

    template<class PairIsph> inline void
    FunctorOuterCorrectPhaseNormal<PairIsph>::enterFor() {
      if (_filter == NULL)
        _pair->error->all(FLERR, "FunctorCorrectPhaseNormal:: Filter is required");
    }

    template<class PairIsph> inline void
    FunctorOuterCorrectPhaseNormal<PairIsph>::operator()(const int ii) {
      const int i = _ilist[ii];
      const int itype = _type[i];

      const int ikind  = _pair->getParticleKind(itype);
      const int iphase = _pair->getParticlePhase(itype);

      const bool has_knormal = util.dotVectors(_dim, _knormal[i], _knormal[i]) > 0.5;
      const bool has_pnormal = util.dotVectors(_dim, _pnormal[i], _pnormal[i]) > 0.5;

      //const double cut = sqrt(LOOKUP(_pair->cutsq, itype, itype));
      //const double h   = LOOKUP(_pair->h, itype, itype);

      if (_filter->yes(ikind) && has_knormal && has_pnormal) {
        const double theta = (iphase == 1 ? _pair->st.csf.theta : M_PI - _pair->st.csf.theta);
        //const double theta = _pair->st.csf.theta;

        double *nw_at_i = _knormal[i];
        double *n_at_i  = _pnormal[i];

        double nt_at_i[3], ntl_at_i[3] = {};
        {

          for (int k=0;k<_dim;++k)
            nt_at_i[k] = n_at_i[k] - util.dotVectors(_dim, n_at_i, nw_at_i)*nw_at_i[k];
        }

        {
          const double mag = sqrt(util.dotVectors(_dim, nt_at_i, nt_at_i));
          if (mag != 0.0)
            for (int k=0;k<_dim;++k)
              nt_at_i[k] /= mag;
        }

        {
          for (int k=0;k<_dim;++k)
            ntl_at_i[k] = nt_at_i[k]*sin(theta) + nw_at_i[k]*cos(theta);
        }

        {
          const double d_at_i = 2.0*(_pnd[i]*_vfrac[i] - 0.5) - 0.5;     // -.5 < d < .5
          const double f_at_i = d_at_i < 0.0 ? 0.0 : 2.0*d_at_i;         // could use pow(2.0*d_at_i, 2);
          for (int k=0;k<_dim;++k)
            n_at_i[k] = (f_at_i)*n_at_i[k] + (1.0 - f_at_i)*ntl_at_i[k];
          
          const double mag = sqrt(util.dotVectors(_dim, n_at_i, n_at_i));
          if (mag != 0.0)
            for (int k=0;k<_dim;++k)
              n_at_i[k] /= mag;
        }
        
      } 
    }

  }
}
#endif

