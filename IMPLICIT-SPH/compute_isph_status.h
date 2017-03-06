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

ComputeStyle(isph/status,ComputeISPH_Status)

#else

#ifndef LMP_COMPUTE_ISPH_STATUS_H
#define LMP_COMPUTE_ISPH_STATUS_H

#include "macrodef.h"
#include "compute.h"

namespace LAMMPS_NS {
  
  class ComputeISPH_Status : public Compute {
  public:
    ComputeISPH_Status(class LAMMPS *, int, char **);
    virtual~ComputeISPH_Status();
    
    void init();
    void compute_vector();
    
  private:
    std::ostream  *_os;

    bool _vflag;
    unsigned int _mask;

    double _r, _c[3];
    double _stat[ISPH_COMPUTE_STATUS_MAX_SIZE];
  };
  
}

#endif
#endif

