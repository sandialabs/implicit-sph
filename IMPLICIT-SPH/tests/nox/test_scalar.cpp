#include "problem_scalar.h"
#include "solver_nox.h"
#include "solver_nox_impl.h"
#include "utils.h"

using namespace std;
using namespace LAMMPS_NS;

typedef class ProblemScalar Problem;
typedef class SolverNOX<ProblemScalar> SolverNonlinear;

int main(int argc, char **argv) {
  int r_val = LAMMPS_SUCCESS;

  MPI_Init(&argc,&argv);
  MPI_Comm comm = MPI_COMM_WORLD;

  Problem  problem;
  SolverNonlinear nl_solver(SolverNonlinear::MatrixFree, problem, comm);
  
  nl_solver.setParameters();
  nl_solver.setConvergenceTests();

  nl_solver.createLinearMap(1, 0);

  double x;
  nl_solver.createSolution(&x);
  nl_solver.setInitialSolution(SolverNonlinear::Value, 1.0);
  
  nl_solver.solveProblem();

  cout << "Solution (x) is : " << x << endl;

  MPI_Finalize();

  return r_val;
}
