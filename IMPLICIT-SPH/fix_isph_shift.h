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

FixStyle(isph/shift,FixISPH_Shift)

#else

#ifndef LMP_FIX_ISPH_SHIFT_H
#define LMP_FIX_ISPH_SHIFT_H

#include "fix.h"

namespace LAMMPS_NS {

  class FixISPH_Shift: public Fix {
  public:
    FixISPH_Shift(class LAMMPS *, int, char **);
    virtual~FixISPH_Shift();

    virtual int setmask();
    virtual void initial_integrate(int);
    virtual void final_integrate();

  protected:

    double _shift;
    double _shiftcut;
    double _nonfluidweight;
  };

}

#endif
#endif
