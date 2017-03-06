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

#ifdef PAIR_CLASS

PairStyle(isph,PairISPH)

#else

#ifndef LMP_PAIR_ISPH_H
#define LMP_PAIR_ISPH_H

#include "pair.h"
#include "random_mars.h"

// time derivative
#include "time_bdf.h"

// Trilinos Epetra, Thyra, and NOX interface
#include "Epetra_CrsMatrix.h"
#include "Thyra_DefaultBlockedLinearOp.hpp"

#include "NOX.H"
#include "NOX_Epetra.H"   

// preconditioner wrapper
#include "precond.h"

// solver linear
#include "solver_lin.h"
#include "solver_lin_belos.h"

// solver nox
#include "solver_nox.h"
#include "solver_nox_aztecOO.h"
#include "solver_nox_stratimikos.h"

// functor base
#include "filter.h"
#include "functor.h"
#include "functor_normalize_vector.h" 
#include "functor_exact_solution.h"

// setup
#include "functor_graph.h"
#include "functor_graph_boundary.h"
#include "functor_electrostatic_force.h"

// color function 
#include "color.h"

// pairwise force function
#include "pairwise_force.h"

// utility fix
#include "fix_isph_ignore_phase_gradient.h"

namespace LAMMPS_NS {
  
  class PairISPH : public Pair,
                   public NOX::Epetra::Interface::Required,
                   public NOX::Epetra::Interface::Jacobian,
                   public NOX::Epetra::Interface::Preconditioner {
  private:
    // ** short name for this class
    typedef class PairISPH                              Self;
    typedef class SolverLin_Belos                       SolverLinear;
    typedef class SolverNOX_Stratimikos<class PairISPH> SolverNonlinear;

    friend class FunctorOuter<class PairISPH>; 
    friend class FunctorOuterNormalizeVector<class PairISPH>;

    friend class FunctorOuterGraph<class PairISPH>;
    friend class FunctorOuterGraphBoundary<class PairISPH>; 
    friend class FunctorOuterElectrostaticForce<class PairISPH>;
    
    // debugging only
    friend class FunctorOuterAdvanceTimeExactSolutionBegin<class PairISPH>;
    friend class FunctorOuterAdvanceTimeExactSolutionEnd<class PairISPH>; 

    // ----------------------------------------------------------------

  public:
    enum DicreteType         { Corrected, MLS };

    enum CommType            { Vfrac,
                               NormalVector, 
                               Velocity, 
                               Pressure, 
                               Vstar, 
                               DeltaP, 
                               Psi,
                               WorkScalar, 
                               WorkVector,
                               TempScalar,
                               TempVector,
                               MaxCommType = 12 };
    
    enum SurfaceTensionModel { ContinuumSurfaceForce = 1,
                               PairwiseForce = 2,
                               MaxSurfaceTensionModel = 3};
    
    enum ParticleKind        { NoParticle = 0, 
                               FluidWithNormal = 1,
                               Fluid = 99,
                               SolidWithNormal = 4,
                               Solid = 12,
                               Boundary = 16,
                               BufferDirichlet = 32,
                               BufferNeumann = 64,
                               All = 127,
                               NormalFilter = (FluidWithNormal + SolidWithNormal + Boundary),
                               MaxParticleKind = 128 };
    
    enum BoundaryCond        { NoBoundaryCond = 0, 
                               ConstExtension = 1, 
                               NavierSlip = 2, 
                               Dirichlet = 3, 
                               MorrisNormal = 4, 
                               MorrisHolmes = 5, 
                               HomogeneousNeumann = 6, 
                               MaxBoundaryCond = 7 };
    
    enum SingularPoisson     { NotSingularPoisson = 0, 
                               NullSpace = 1, 
                               PinZero = 2, 
                               DoubleDiag = 3,
                               MaxSingularPoisson = 4 };
    
  public:
    PairISPH(class LAMMPS *);
    virtual~PairISPH();

    // interface is declared in base pair
    virtual void compute(int, int); // pure, so must be implemented

    // Question: what are these for ?
    //void compute_inner();
    //void compute_middle();
    //void compute_outer(int, int);

    virtual double single(int, int, int, int, double, double, double, double &);

    virtual void settings(int, char **);
    virtual void coeff(int, char **);

    virtual void init_style();
    virtual double init_one(int, int);
    //void init_list(int, class NeighList *);

    //void write_restart(FILE *);
    //void read_restart(FILE *);

    //void write_restart_settings(FILE *);
    //void read_restart_settings(FILE *);
 
    // overload comm routines to specialize them for this pairing
    virtual int pack_forward_comm(int, int *, double *, int, int *);
    virtual void unpack_forward_comm(int, int, double *);

    //int pack_reverse_comm(int, int, double *);
    //void unpack_reverse_comm(int, int *, double *);
  
    //double memory_usage();

    // - ISPH -
    int  getDiscretizationInfo() const;
    int  getParticleKind(const int itype) const; //, const int i = -1) const;
    void resetParticleKindMap(); 
    void mapParticleKind(const int from, const int to); 
    bool isParticleFixed(const int itype) const;
    int  getParticlePhase(const int itype) const;
    bool isAleEnabled() const;

    double getDefaultCutSq() const;
    double** getVstar() const;

    void setUsePartToComputeNormal(const bool flag);

    void clearCommArray(const int comm_type); 

    virtual void printParticles(ostream &s, string name);
    virtual void grow(double extend = 1.0);

    virtual void initializeSolvers();

    virtual void computeGraph();
    virtual void computePre(const int mask = 0);
    virtual void computeNormals();
    virtual double computeZeroMeanPressure(double *f, bool update = true);

    virtual void refreshParticles();

    virtual void modifySingularMatrix(const int row, double &diag, double &b);

    virtual void computeUnitTests();

    // PB - NOX interface
    virtual bool computeF(const Epetra_Vector &x, Epetra_Vector &f,
                          const FillType flag = Residual);
    virtual bool computeJacobian(const Epetra_Vector &x, Epetra_Operator &J);
    virtual bool computePreconditioner(const Epetra_Vector &x, Epetra_Operator &M,
                                       Teuchos::ParameterList *param = NULL);

    virtual void computePsiGradient();
    virtual void computePoissonBoltzmann();  /* main procedure of PB */

    virtual void computeAppliedElectricPotential(double *b);
    virtual void computePhiGradient();
    virtual void computeAppliedElectricField();

    virtual void computeElectrostaticForce();

    // ST - Continuum Surface Force (CSF) for surface tension
    virtual void computeSurfaceTension();

    // RS - Random tensor for fluctuating hydrodynamics
    virtual void computeRandomStressTensor(double **rstress_x, double **rstress_y, double **rstress_z);
    virtual void computeForceFromRandomStress();

    // RT - Transport equation
    virtual void computeSoluteTransport();
    virtual void computeSoluteTransportSpecies(const double dcoeff, double *b);

    // NS - Projection
    virtual void computeHelmholtz(double *b);
    virtual void computeBlockHelmholtz(double *b);
    virtual void computePoisson(double *b);
    virtual void correctVelocity();
    virtual void correctPressure();

    virtual void computeIncompressibleNavierStokes();  /* main procedure of NS */

    // NS - ALE
    virtual void predictAleVelocity();
    virtual void computeAlePoisson(double *b);
    virtual void correctAleVelocity();
    virtual void computeAleHelmholtz(double *b);

    virtual void computeAleIncompressibleNavierStokes();  /* main procedure of NS ALE */

    // NS - Time/Shift
    virtual void advanceTime();
    virtual void shiftParticles(const double shift,
                                const double shiftcut,
                                const double nonfluidweight);

    // NS - Post process
    virtual void computeVelocityCurl(double **v, double **curl);
    virtual void computeVelocityDivergence(double **v, double *div);
    virtual void computeTractionVector(double *p, double **v, double **traction);

  protected:

    virtual void allocate();

    // ----------------
    // user attributes
    // ----------------

    // test mode
    bool use_exact_solution;
    
    // discretization information: Corrected, MLS
    int discretization;

    // particle information
    int **pinfo, boundary_particle, particle_kind_map[MaxParticleKind];

    // lookup table for smoothing lengths: [itype][jtype]
    double **h;  

    // normal to fluid/structure interface and particle number density
    double **normal;

    // use part information to compute normal
    bool use_part;

    // intermediate vectors for the two-step scheme
    double **vstar, *dp;

    // Unit test
    struct {
      bool is_enabled;
      bool test_compact_poisson_gradient;
    } ut;

    // Applied Efield
    struct {
      bool is_enabled;
      bool is_solid2fluid;
      int smooth_phi;
      double e[3];
    } ae;

    // Poisson Boltzmann parameters
    struct {
      bool is_enabled, is_linearized;
      double psiref;
      double ezcb, e[3];
      double gamma;

      int boundary;

      bool has_extra_f;
      Teuchos::ParameterList extra_f;
    } pb;

    // Surface tension parameter
    struct {
      bool is_enabled;
      int model;

      struct {
        double alpha;  // surface tension coefficient
        double theta;  // solid contact angle 
        
        ColorFunction *color; // color function
        double epsilon;       // normal threshold
        double kappa;         // maximum curvature
      } csf;
      
      struct {
        PairwiseForceFunction *force; // interaction model
        double **s;
      } pf;
    } st;

    // Navier Stokes parameters
    struct {
      bool is_enabled;
      
      double theta;    // theta time stepping
      double beta;     // Wenxiao's BC parameter
      double g[3];     // constant body force field
      double kBT;      // Wenxiao's kinetic temperature

      int boundary;
      int singular_poisson;

      bool is_block_helmholtz_enabled;
      bool is_incremental_pressure_used;
      bool is_momentum_preserve_operator_used;

      bool is_fluctuation_enabled;
      bool is_bond_interaction_enabled;
    } ns;

    struct {
      bool is_enabled;
      TimeBDF *tdiff;
    } ale;

    // Solute Transport
    struct {
      bool is_enabled;

      bool mask[ISPH_MAX_CONCENTRATION];     // species mask
      double dcoeff[ISPH_MAX_CONCENTRATION]; // scalar diffusion coefficients
      double theta;                          // theta time stepping

      bool is_particle_fixed;
    } tr;

    // solver interface to NOX
    SolverNonlinear *nl_solver; 

    // preconditioner
    PrecondWrapper *prec;

    // solver interface to Belos or stratimikos
    SolverLinear *li_solver;

    // sparse matrix
    Epetra_Map *nodalmap;
    Epetra_CrsGraph *tags_in_cut, *tags_in_boundary;
      
    struct {
      Epetra_CrsMatrix *crs;
      Epetra_Vector *diagonal, *scaled_laplace_diagonal;
      int is_filled;
    } A;

    struct {
      Epetra_CrsMatrix *crs[9];
      Epetra_Vector *diagonal[3], *scaled_laplace_diagonal[3];
      void set(const int i, const int j, Epetra_CrsMatrix *aij) { crs[i + j*3] = aij; }
      Epetra_CrsMatrix* operator()(const int i, const int j) const { return crs[i + j*3]; }
    } A_blk;

    int is_once, nmax;

    // select communication variables, e.g., VFRAC and PSI
    CommType comm_variable;

    // LAMMPS Marsaglia random number generator
    class RanMars *random;

  public:

    // workspace that can be accessed from outside
    double *work, **work3;

    // temporary pointer holder for communication
    double *temp, **temp3;
  };

}

#endif
#endif
