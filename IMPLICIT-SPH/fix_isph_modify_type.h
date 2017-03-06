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

FixStyle(isph/modify/type,FixISPH_ModifyType)

#else

#ifndef LMP_FIX_ISPH_MODIFY_TYPE_H
#define LMP_FIX_ISPH_MODIFY_TYPE_H

#include "fix.h"

namespace LAMMPS_NS {

  class FixISPH_ModifyType: public Fix {
  public:
    FixISPH_ModifyType(class LAMMPS *, int, char **);
    virtual~FixISPH_ModifyType();

    virtual int setmask();
    virtual void initial_integrate(int);

  protected:

    int _type_from, _type_to;
    double _xmin, _xmax, _ymin, _ymax, _zmin, _zmax;
  };

}

#endif
#endif
