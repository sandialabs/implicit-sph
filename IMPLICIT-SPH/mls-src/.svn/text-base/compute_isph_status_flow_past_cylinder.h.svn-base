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

ComputeStyle(isph/status/flow_past_cylinder,ComputeISPH_StatusFlowPastCylinder)

#else

#ifndef LMP_COMPUTE_ISPH_STATUS_FLOW_PAST_CYLINDER_H
#define LMP_COMPUTE_ISPH_STATUS_FLOW_PAST_CYLINDER_H

#include "macrodef.h"
#include "compute.h"

namespace LAMMPS_NS {
  
  class ComputeISPH_StatusFlowPastCylinder : public Compute {
  public:
    ComputeISPH_StatusFlowPastCylinder(class LAMMPS *, int, char **);
    virtual~ComputeISPH_StatusFlowPastCylinder();
    
    void init();
    void compute_vector();
    
  private:
    // file io 
    std::ostream *_os;

    // cylinder surface integration 
    int _axis;
    double _cyl_radius, _cyl_lo, _cyl_hi, _cyl_center[3];

    // direction vector for drag and lift
    double _drag[3], _lift[3];

    // type mask for checking points on surface
    unsigned int _mask;

    // status
    double _stat[ISPH_COMPUTE_STATUS_MAX_SIZE];

    // temporary storage for traction vector
    int _nmax;
    double **_traction;
  };
  
}

#endif
#endif

