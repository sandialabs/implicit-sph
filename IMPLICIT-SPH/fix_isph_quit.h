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

FixStyle(isph/quit,FixISPH_Quit)

#else

#ifndef LMP_FIX_ISPH_QUIT_H
#define LMP_FIX_ISPH_QUIT_H

#include "fix.h"

namespace LAMMPS_NS {

  class FixISPH_Quit: public Fix {
    
  public:
    FixISPH_Quit(class LAMMPS *, int, char **);
    virtual~FixISPH_Quit();

    virtual int setmask();
    virtual void final_integrate();

  protected:
    class Compute *_stat;
    double _tfinal;
  };

}

#endif
#endif
