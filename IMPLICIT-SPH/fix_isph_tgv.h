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

FixStyle(isph/tgv,FixISPH_TGV)

#else

#ifndef LMP_FIX_ISPH_TGV_H
#define LMP_FIX_ISPH_TGV_H

#include "fix.h"
#include "utils.h"

namespace LAMMPS_NS {

  class FixISPH_TGV: public Fix {
  public:
    FixISPH_TGV(class LAMMPS *, int, char **);
    virtual ~FixISPH_TGV();

    int setmask();
    virtual void final_integrate();
    
    template<class ArgPair> void compute_error(ArgPair *pair);

  protected:
    class Compute *_stat;

  };


}

#endif
#endif

