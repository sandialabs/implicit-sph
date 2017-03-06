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

ComputeStyle(isph/sphere/porous,ComputeISPH_SpherePorous)

#else

#ifndef LMP_COMPUTE_ISPH_SPHERE_POROUS_H
#define LMP_COMPUTE_ISPH_SPHERE_POROUS_H

#include "math.h"
#include "compute.h"

namespace LAMMPS_NS {
  
  class ComputeISPH_SpherePorous : public Compute {
  public:
    ComputeISPH_SpherePorous(class LAMMPS *, int, char **);
    virtual~ComputeISPH_SpherePorous();
    
    void init();

    void readSphereCoordinates(const char *filename);
    void setParticleGeometry();
    
  private:
    int _n_spheres, _stype;
    double **_sphere_coords, *_sphere_radius, _sphere_hollow_gap;

    inline bool is_coords_in_spheres(double *pt, int &pid) {

      for (int i=0;i<_n_spheres;++i) {
        const double r_at_i = _sphere_radius[i];
        double diff = r_at_i - _sphere_hollow_gap;
        double no_sq  = (diff > 0.0 ? diff*diff : 0.0);
        double yes_sq = r_at_i*r_at_i;

        double *st = _sphere_coords[i];
        double dist_sq = (pow(pt[0] - st[0], 2) +
                          pow(pt[1] - st[1], 2) +
                          pow(pt[2] - st[2], 2));
        if (dist_sq < yes_sq) {
          pid = (no_sq < dist_sq ? i+1 : 0); // positive part number
          return true;
        }
      }
      pid = 0;
      return false;
    }

  };
  
}

#endif
#endif

