/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef BOND_CLASS

BondStyle(isph/harmonic,BondISPH_Harmonic)

#else

#ifndef LMP_BOND_ISPH_HARMONIC_H
#define LMP_BOND_ISPH_HARMONIC_H

#include "bond_harmonic.h"
#include "bond_isph.h"

namespace LAMMPS_NS {
  
  class BondISPH_Harmonic : public BondHarmonic,
                            public BondISPH {
  public:
    BondISPH_Harmonic(class LAMMPS *);
    virtual ~BondISPH_Harmonic();
    virtual void compute(int eflag, int vflag);
  };
  
}

#endif
#endif

/* ERROR/WARNING messages:

E: Incorrect args for bond coefficients

Self-explanatory.  Check the input script or data file.

*/
