/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS

ComputeStyle(isph/velocity_divergence,ComputeISPH_VelocityDivergence)

#else

#ifndef LMP_COMPUTE_ISPH_VELOCITY_DIVERGENCE_H
#define LMP_COMPUTE_ISPH_VELOCITY_DIVERGENCE_H

#include "compute.h"

namespace LAMMPS_NS {
  
  class ComputeISPH_VelocityDivergence : public Compute {
  public:
    ComputeISPH_VelocityDivergence(class LAMMPS *, int, char **);
    virtual~ComputeISPH_VelocityDivergence();
    
    void init();
    void compute_peratom();
    
  private:
    int _nmax, _vflag;
    double *_div;
  };
  
}

#endif
#endif

