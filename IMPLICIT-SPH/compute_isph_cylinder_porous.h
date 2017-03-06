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

ComputeStyle(isph/cylinder/porous,ComputeISPH_CylinderPorous)

#else

#ifndef LMP_COMPUTE_ISPH_CYLINDER_POROUS_H
#define LMP_COMPUTE_ISPH_CYLINDER_POROUS_H

#include "math.h"
#include "compute.h"

namespace LAMMPS_NS {
  
  class ComputeISPH_CylinderPorous : public Compute {
  public:
    ComputeISPH_CylinderPorous(class LAMMPS *, int, char **);
    virtual~ComputeISPH_CylinderPorous();
    
    void init();

    void readBeadCoordinates(const char *filename);
    void setParticleGeometry();
    
  private:
    int _bead_type, _n_beads;
    double _bead_coords_scale, **_bead_coords, _bead_radius, _bead_inner_radius, _bead_coords_hi, _bead_coords_lo;
    
    int _axis;
    double _cyl_hi, _cyl_lo, _cyl_buffer_min, _cyl_buffer_max, _cyl_radius, _cyl_center[3];

    inline bool is_coords_in_beads(double *pt, int &pid) {
      double no_sq  = pow(_bead_inner_radius, 2);
      double yes_sq = pow(_bead_radius, 2);
      for (int i=0;i<_n_beads;++i) {
        double *bt = _bead_coords[i];

        if (_bead_coords_lo < bt[_axis] && bt[_axis] < _bead_coords_hi) {
          double dist_sq = (pow(pt[0] - bt[0], 2) +
                            pow(pt[1] - bt[1], 2) +
                            pow(pt[2] - bt[2], 2));
          if (dist_sq < yes_sq) {
            pid = (no_sq < dist_sq ? i+1 : 0); // positive part number
            return true;
          }
        }
      }
      pid = 0;
      return false;
    }

    inline bool is_coords_in_cylinder(double *pt) {
      double yes_sq = pow(_cyl_radius, 2);
      double dist_sq = 0.0;
      for (int i=0;i<3;++i) 
        dist_sq += (i == _axis ? 0.0 : pow(pt[i] - _cyl_center[i], 2));
      return (dist_sq < yes_sq);
    }

    inline bool is_coords_in_buffer(double *pt) {
      return (_cyl_buffer_min < pt[_axis] && pt[_axis] < _cyl_buffer_max); 
    }
  };
  
}

#endif
#endif

