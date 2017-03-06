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

#ifndef LMP_BOND_ISPH_H
#define LMP_BOND_ISPH_H

#include "stdio.h"

namespace LAMMPS_NS {
  
  class BondISPH {
  protected:
    bool is_compute_enabled;
    
  public:
    BondISPH() : is_compute_enabled(false) {}
    virtual ~BondISPH() {}
    void setComputeEnabled(const bool flag) {
      is_compute_enabled = flag;
    }
    bool isComputeEnabled() const {
      return is_compute_enabled;
    }
  };
  
}

#endif


