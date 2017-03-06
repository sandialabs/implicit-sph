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

#ifdef ATOM_CLASS

AtomStyle(isph,AtomVecISPH)

#else

#ifndef LMP_ATOM_VEC_ISPH_H
#define LMP_ATOM_VEC_ISPH_H

#include "atom_vec.h"

namespace LAMMPS_NS {

class AtomVecISPH: public AtomVec {
public:
  AtomVecISPH(class LAMMPS *);
  void grow(int);
  void grow_reset();
  void copy(int, int, int);
  int pack_comm(int, int *, double *, int, int *);
  int pack_comm_vel(int, int *, double *, int, int *);
  void unpack_comm(int, int, double *);
  void unpack_comm_vel(int, int, double *);
  int pack_reverse(int, int, double *);
  void unpack_reverse(int, int *, double *);
  int pack_border(int, int *, double *, int, int *);
  int pack_border_vel(int, int *, double *, int, int *);
  void unpack_border(int, int, double *);
  void unpack_border_vel(int, int, double *);
  int pack_exchange(int, double *);
  int unpack_exchange(double *);
  int size_restart();
  int pack_restart(int, double *);
  int unpack_restart(double *);
  void create_atom(int, double *);
  void data_atom(double *, imageint, char **);
  bigint memory_usage();

  void pack_data(double **);
  void write_data(FILE *, int, double**);

private:
  tagint *tag;
  int  *type, *mask;
  imageint *image;

  // variables
  // - part id
  int *part;

  // - position, velocity at timestep n, velocity at timestep n+1, pressure
  double **x, **v, **f, *pressure;

  // - concentration
  double **concentration;

  // - ale particle velocity \dot{x}
  double **xdot;

  // - velocity history for high order time stepping
  double ***vprev, ***xprev;

  // - material parameters
  double *rmass, *viscosity, *vfrac, *density;

  // - electric potential and field
  double *eps, *sigma;
  double *psi, *psi0, **psigrad, *phi, **phigrad;

  // - bond structure
  tagint *molecule;
  int **nspecial;
  tagint **special;
  int *num_bond;
  int **bond_type;
  tagint **bond_atom;

};

}

#endif
#endif

/* ERROR/WARNING messages:

 E: Per-processor system is too big

 The number of owned atoms plus ghost atoms on a single
 processor must fit in 32-bit integer.

 E: Invalid atom ID in Atoms section of data file

 Atom IDs must be positive integers.

 E: Invalid atom type in Atoms section of data file

 Atom types must range from 1 to specified # of types.

 E: Invalid mass value

 Self-explanatory.

 */
