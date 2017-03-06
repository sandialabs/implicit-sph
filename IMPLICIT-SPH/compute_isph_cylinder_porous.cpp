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

#include "string.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "error.h"

#include <iostream>
#include <fstream>

#include "math.h"
#include "utils.h"

#include "pair_isph_corrected.h"
//#include "pair_isph_mls.h"

#include "compute_isph_cylinder_porous.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "ComputeISPH_CylinderPorous::ComputeISPH_CylinderPorous"
ComputeISPH_CylinderPorous::ComputeISPH_CylinderPorous(LAMMPS *lmp, int narg, char **arg) 
  : Compute(lmp, narg, arg) {
  FUNCT_ENTER(comm->me);

  if (atom->isph_flag != 1)
    error->all(FLERR, "Compute ipsh/status command requires a matching atom_style, i.e., isph");

  if (narg < 16)
    error->all(FLERR,"Illegal number of arguments for compute isph command");

  // arg[3] porous type : uniform_sphere only
  _bead_type = 0;
  
  // arg[4] scale
  _bead_coords_scale = force->numeric(FLERR, arg[4]);
  
  // arg[5] bead radius
  _bead_radius = force->numeric(FLERR, arg[5]);
  {
    PairISPH *pair = dynamic_cast<PairISPH*>(force->pair_match("isph", 0));
    _bead_inner_radius = max(_bead_radius - 1.0*sqrt(pair->getDefaultCutSq()), 0.0);
  }

  // arg[6,7] bead lower upper limits
  _bead_coords_lo = force->numeric(FLERR, arg[6]);
  _bead_coords_hi = force->numeric(FLERR, arg[7]);

  // arg[8] porous file : coordinates 
  _bead_coords = NULL;
  
  // arg[9,10,11,12] cylinder hi and row
  _cyl_lo = force->numeric(FLERR, arg[9]);
  _cyl_hi = force->numeric(FLERR, arg[10]);

  _cyl_buffer_min = force->numeric(FLERR, arg[11]);  
  _cyl_buffer_max = force->numeric(FLERR, arg[12]);  

  // arg[13] cylinder radius
  _cyl_radius = force->numeric(FLERR, arg[13]);

  // arg[14,15,16] cylinder center
  _cyl_center[0] = force->numeric(FLERR, arg[14]);
  _cyl_center[1] = force->numeric(FLERR, arg[15]);  
  _cyl_center[2] = force->numeric(FLERR, arg[16]);  

  // cylinder axis component
  if      (strcmp(arg[17],"x") == 0) _axis = 0;
  else if (strcmp(arg[17],"y") == 0) _axis = 1;
  else if (strcmp(arg[17],"z") == 0) _axis = 2;
  else
    error->all(FLERR,"Not supported keyword: x, y, z");

  if (comm->me == 0) {
    cout << " - ComputeISPH_CylinderPorous - " << endl
         << "   bead data file = " << arg[8] << endl
         << "   bead type = " << _bead_type << endl
         << "   bead coordinate scale = " << _bead_coords_scale << endl
         << "   bead radius = " << _bead_radius << endl
         << "   bead limits = " << _bead_coords_lo << ", " << _bead_coords_hi << endl
         << "  " << endl
         << "   cylinder radius = " << _cyl_radius << endl
         << "   cylinder low and high = " << _cyl_lo << ", " <<  _cyl_hi << endl
         << "   cylinder buffer = " << _cyl_buffer_min << ", " << _cyl_buffer_max << endl
         << "   cylinder center = "
         << _cyl_center[0] << ", " << _cyl_center[1] << ", " << _cyl_center[2] << endl
         << "   cylinder axis = " << arg[15] << ", " << _axis << endl;
  }

  readBeadCoordinates(arg[8]);
  setParticleGeometry();
  
  FUNCT_EXIT(comm->me);
}

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "ComputeISPH_CylinderPorous::~ComputeISPH_CylinderPorous"
ComputeISPH_CylinderPorous::~ComputeISPH_CylinderPorous() {
  FUNCT_ENTER(comm->me);

  memory->destroy(_bead_coords);
  
  FUNCT_EXIT(comm->me);  
}

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "ComputeISPH_CylinderPorous::init"
void ComputeISPH_CylinderPorous::init() {
  FUNCT_ENTER(comm->me);

  int count = 0;
  for (int i=0;i<modify->ncompute;++i)
    count += (strcmp(modify->compute[i]->style,"isph/cylinder/porous") == 0);

  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute isph/cylinder/porous");

  // compute geometry in the constructor so that we can delete some
  // setParticleGeometry();
  
  FUNCT_EXIT(comm->me);  
}

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "ComputeISPH_CylinderPorous::readBeadCoordinates"
void ComputeISPH_CylinderPorous::readBeadCoordinates(const char *filename) {
  FUNCT_ENTER(comm->me);

  double x[3];
  
  // count beads    
  _n_beads = 0;
  for (ifstream ifs(filename);!ifs.eof();++_n_beads) 
    ifs >> x[0] >> x[1] >> x[2];

  // allocate
  if (_n_beads) 
    memory->grow(_bead_coords, _n_beads, 3, "compute:cylinder:bead_coords");

  // read bead coordinates
  ifstream ifs(filename);
  int cnt = 0;
  for (int i=0;i<_n_beads;++i) {
    ifs >> x[0] >> x[1] >> x[2];
    
    for (int j=0;j<3;++j)
      x[j] *=_bead_coords_scale;

    double r = 0.0;
    for (int j=0;j<3;++j)
      if (j != _axis)
        r += pow(x[j] - _cyl_center[j], 2);
    r = sqrt(r);

    //if (r < (_cyl_radius - _bead_radius)) {
      for (int j=0;j<3;++j)
        _bead_coords[cnt][j] = x[j];
      ++cnt;
      //}
  }
  
  if (comm->me == 0) {
    cout << " - ComputeISPH_CylinderPorous - " << endl
         << "   # of beads in file = " << _n_beads << endl
         << "   # of beads counted = " << cnt << endl;
  }
  _n_beads = cnt;

  FUNCT_EXIT(comm->me);  
}

/* ---------------------------------------------------------------------- */
#undef  __FUNCT__
#define __FUNCT__ "ComputeISPH_CylinderPorous::setParticleGeometry"
void ComputeISPH_CylinderPorous::setParticleGeometry() {
  FUNCT_ENTER(comm->me);
  {
    // get accesses to data structures
    const int nsize = (atom->nlocal + atom->nghost);

    double **x = atom->x;
    int *type = atom->type;

    PairISPH *pair = dynamic_cast<PairISPH*>(force->pair_match("isph", 0));
    if (pair != NULL) {
      pair->setUsePartToComputeNormal(true);
      int *part = atom->part;

      // check all particles including ghost particles
      for (int i=0;i<nsize;++i) {
        int pid = 0;
        double *x_at_i = x[i];
        if       (!is_coords_in_cylinder(x_at_i))   { type[i] = 4;             part[i] = -1;  // out of cylinder
        } else if (is_coords_in_beads(x_at_i, pid)) { type[i] = (pid ? 3 : 0); part[i] = pid; // inside beads
        } else if (is_coords_in_buffer(x_at_i))     { type[i] = 2;             part[i] = 0;   // buffer fluids
        } else {                                      type[i] = 1;             part[i] = 0;   // fluids 
        }
      }

      // if type == 0 they should be erased
    } else {
      error->all(FLERR, "Fail to retrive a pointer to pair: PairISPH"); 
    }
  }

  {
    double part_local[5] = {};
    double part_global[5] = {};

    const int nlocal = atom->nlocal;
    double **x = atom->x;
    int *type = atom->type;

    // check all particles including ghost particles
    for (int i=0;i<nlocal;++i) {
      for (int j=0;j<5;++j)
        part_local[j] += (type[i] == j);
      if (type[i] == 1) { // fluid case
        double *x_at_i = x[i];
        if (x_at_i[_axis] < _cyl_lo || x_at_i[_axis] > _cyl_hi) {
          --part_local[1];
          ++part_local[2];
        }
      }
    }

    // compute porosity
    MPI_Allreduce(part_local, part_global, 5, MPI_DOUBLE, MPI_SUM, world);

    const int w = 15;
    if (comm->me == 0) {
      cout << "ComputeIsph_CylinderPorous::Parts = "
           << part_global[0] << " "
           << part_global[1] << " "
           << part_global[2] << " "
           << part_global[3] << " "
           << part_global[4] << " "
           << endl;

      cout << "ComputeIsph_CylinderPorous::Porosity = " 
           << scientific
           << setw(w) << ((part_global[1])/(part_global[0]+part_global[1]+part_global[3]))
           << fixed << endl;
    }

  }
  FUNCT_EXIT(comm->me);  
}

