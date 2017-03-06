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

ComputeStyle(isph/applied-electric-potential/henry,ComputeISPH_AppliedElectricPotentialHenry)

#else

#ifndef LMP_COMPUTE_ISPH_APPLIED_ELECTRIC_POTENTIAL_HENRY_H
#define LMP_COMPUTE_ISPH_APPLIED_ELECTRIC_POTENTIAL_HENRY_H

#include "compute.h"

namespace LAMMPS_NS {
  
  class ComputeISPH_AppliedElectricPotentialHenry : public Compute {
  public:
    ComputeISPH_AppliedElectricPotentialHenry(class LAMMPS *, int, char **);
    virtual~ComputeISPH_AppliedElectricPotentialHenry();
    
    void init();

    double compute_scalar(); // L2 error
    void compute_peratom();  // put exact solution on particles
    void computeHenrySolution(const int dim, const double *dx,
                              double &phi_at_i, double *phigrad_at_i);
  private:
    double _c[3], _a, _sratio, _eapp;
  };
  
}

#endif
#endif

