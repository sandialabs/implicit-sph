#pragma once
#ifndef __FUNCTOR_VOLUME_H__
#define __FUNCTOR_VOLUME_H__

#include <iostream>
#include <string>
#include "utils.h"
#include "functor.h"

namespace LAMMPS_NS {

  namespace Corrected {

    template<class PairIsph>
    class FunctorOuterVolume : public FunctorOuter<PairIsph> {
    public:
      FunctorOuterVolume(PairIsph *isph) : FunctorOuter<PairIsph>(isph) { }
      void enterFor();
      void operator()(const int ii);
      void exitFor();

    protected:
      FUNCTOR_REMOVE_THIS_POINTER;
    };

    template<class PairIsph> inline 
    void FunctorOuterVolume<PairIsph>::enterFor() {
      // at this moment, filter is not allowed in computing volume
      // if we need separate set of quarature point remove the error check
      // and add a proper filter
      if (_filter != NULL)
        _pair->error->all(FLERR, "FunctorVolume:: Filter is not allowed - it computes volumes for all particles interacting with all neighbor particles)");      
    }  

    // compute volume fraction (quadrature) for all particles interacting 
    // with all particle types. when this quadrature is used in the interface,
    // this may result inaccuracy. however, the continuity on the normal direction. 
    // implies that the normal component at the interface is zero and two particle
    // would have the same values. this eventually neglect the contribution from 
    // the other type of particles. All in all, all is well and don't worry now.
    template<class PairIsph> inline 
    void FunctorOuterVolume<PairIsph>::operator()(const int ii) {
      int i = _ilist[ii];
      int itype = _type[i];

      // initialize
      _vfrac[i] = 0.0;

      // neighbor around i
      int *jlist = _firstneigh[i];
      int jnum = _numneigh[i];

      // initialize wtmp to include self contribution
      double wtmp = _pair->kernel->val(0.0, LOOKUP(_pair->h, itype, itype));

      for (int jj=0;jj<jnum;++jj) {
        // j th particle; always put neighmask two bits are reserved for LAMMPS
        int j = (jlist[jj] & NEIGHMASK);
        int jtype = _type[j];

        // compute distance between partile i and j
        double rsq = 0.0;
        for (int k=0;k<_dim;++k) {
          double r = _x[i][k] - _x[j][k];
          rsq += (r*r);
        }

        // integrate
        if (rsq < LOOKUP(_pair->cutsq, itype, jtype)) 
          wtmp += _pair->kernel->val(sqrt(rsq), LOOKUP(_pair->h, itype, jtype));
      }
      // compute volume fraction for ith particle
      _vfrac[i] = 1.0/wtmp;
    }
    template<class PairIsph> inline 
    void FunctorOuterVolume<PairIsph>::exitFor() {
      // communicate volumes
      _pair->comm_variable = PairIsph::Vfrac;
      _pair->comm_forward = 1;
      _pair->comm->forward_comm_pair(_pair);
    }

  }
}
#endif

