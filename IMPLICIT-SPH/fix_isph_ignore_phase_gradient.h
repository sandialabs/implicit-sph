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

#ifdef FIX_CLASS

FixStyle(isph/ignore/phasegradient,FixISPH_IgnorePhaseGradient)

#else

#ifndef LMP_FIX_ISPH_IGNORE_PHASE_GRADIENT_H
#define LMP_FIX_ISPH_IGNORE_PHASE_GRADIENT_H

#include "fix.h"

namespace LAMMPS_NS {

  class FixISPH_IgnorePhaseGradient: public Fix {
  public:
    FixISPH_IgnorePhaseGradient(class LAMMPS *, int, char **);
    virtual~FixISPH_IgnorePhaseGradient();

    virtual int setmask();

    void ignorePhaseGradient(double **phase_grad, const double cut) const;
  protected:

    int _axis;
    double **_phase_grad, _pt_on_axis, _thres_over_cut;
  };

}

#endif
#endif
