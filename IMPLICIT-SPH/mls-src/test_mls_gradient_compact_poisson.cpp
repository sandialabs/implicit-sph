#include <iostream>
#include <string>
#include "utils.h"

#include "stdio.h"
#include "string.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "update.h"
#include "respa.h"
#include "error.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "pair.h"
#include "domain.h"
#include "math.h"
#include "memory.h"

#include "pair_isph_mls.h"
#include "test_mls.h"

namespace LAMMPS_NS {

  namespace MLS {

    using namespace std;

    typedef class PairISPH_MLS PairIsph;
    typedef class FunctorOuterGradientCompactPoisson<PairIsph> FunctorGradientCompactPoisson;
    
    int TestSuite::testMLS_GradientCompactPoisson() {
      int r_val = 0;

      // Given 2D periodic domain, this problem check gradient operator with 
      // manufactured solution u = sin(x) sin(y).

      // ** Initialize workspaces and setup necessary variables 
      double *u = NULL, **u_grad = NULL;
      double *f = NULL, *g = NULL;
      {
        const int nmax = _pair->atom->nmax;
        _pair->memory->grow(u,      nmax,    "test:mls:u");       // exact function
        _pair->memory->grow(f,      nmax,    "test:mls:f");       // \laplace u = f
        _pair->memory->grow(g,      nmax,    "test:mls:g");       // \grad u \normal = g

        _pair->memory->grow(u_grad, nmax, 3, "test:mls:u_grad");  // computed gradient of u
      }

      int     *type   = _pair->atom->type;      
      double **x      = _pair->atom->x;
      double **normal = _pair->normal;

      // ** Set up exact solution
      {
        const int nsize = _pair->atom->nlocal;
        for (int i=0;i<nsize;++i) {
          // u = sin(x)sin(y)
          u[i] = sin(x[i][0])*sin(x[i][1]);  

          // interior penalty f = \laplace u
          f[i] = -2.0*sin(x[i][0])*sin(x[i][1]);

          // boundary penalty g = \grad u \cdot normal
          if (_pair->getParticleKind(type[i]) == PairIsph::Boundary)
            g[i] = (cos(x[i][0])*sin(x[i][1])*normal[i][0] +
                    sin(x[i][0])*cos(x[i][1])*normal[i][1]);
          else
            g[i] = 0.0;
        }

        // communicate
        _pair->temp = u;
        _pair->comm_variable = PairIsph::TempScalar;
        _pair->comm_forward = 1;
        _pair->comm->forward_comm_pair(_pair);

        _pair->temp = f;
        _pair->comm_variable = PairIsph::TempScalar;
        _pair->comm_forward = 1;
        _pair->comm->forward_comm_pair(_pair);

        _pair->temp = g;
        _pair->comm_variable = PairIsph::TempScalar;
        _pair->comm_forward = 1;
        _pair->comm->forward_comm_pair(_pair);

        _pair->temp = NULL;
      }
      
      // ** Test body
      {
        _pair->enterCompactPoisson();
        FunctorGradientCompactPoisson functor(_pair, u, 1.0, f, normal, g, u_grad);
        PairFor(functor, functor.getNumberOfWork());
        _pair->exitCompactPoisson();
      }

      // ** Evaluation
      {
        const int nsize = _pair->atom->nlocal;

        double l_err[2] = {}, g_err[2] = {};
        for (int i=0;i<nsize;++i) {
          const double diff_x = cos(x[i][0])*sin(x[i][1]) - u_grad[i][0];
          const double diff_y = sin(x[i][0])*cos(x[i][1]) - u_grad[i][1];

          l_err[0] += (pow(diff_x,2) + pow(diff_y,2));

          // set error into u_grad for visualization
          u_grad[i][0] = abs(diff_x);
          u_grad[i][1] = abs(diff_y);
        }
        l_err[1] = nsize;

        // sum up over all processors
        MPI_Reduce(l_err, g_err, 2, MPI_DOUBLE, MPI_SUM, 0, _pair->world); 
        
        if (_pair->comm->me == 0) {
          const double err = sqrt(g_err[0]/g_err[1]);
          cout << "Error = " << err << endl;
          if (err > _threshold) 
            r_val = 1;
        }

        // dump to psigrad
        memcpy(&_pair->atom->psigrad[0][0], &u_grad[0][0], sizeof(double)*3*nsize);
      }

      // ** Finalize 
      {
        _pair->memory->destroy(u);
        _pair->memory->destroy(f);
        _pair->memory->destroy(g);

        _pair->memory->destroy(u_grad);
      }
      return r_val;
    }

  }
}

