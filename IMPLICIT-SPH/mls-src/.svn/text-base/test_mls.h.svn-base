#pragma once
#ifndef __TEST_MLS_H__
#define __TEST_MLS_H__

#include "atom.h"
#include "comm.h"

#undef __UNIT_TEST__
#define __UNIT_TEST__(test, r_value, count) {         \
    const int r_val_test = test;                      \
    r_value += r_val_test; ++count;                   \
                                                      \
    if (_pair->comm->me == 0) {                       \
      cout << setw(32) << #test ;                     \
      if (r_val_test != 0)                            \
        cout << ":: FAILED" << endl;                  \
      else                                            \
        cout << ":: PASSED" << endl;                  \
    }                                                 \
  }

namespace LAMMPS_NS {

  class PairISPH_MLS;

  namespace MLS {

    class TestSuite {
    private:
      PairISPH_MLS *_pair;
      double _threshold;

    public:
      TestSuite(PairISPH_MLS *pair, 
                double threshold = 1.0e-6) 
        : _pair(pair),
          _threshold(threshold) { } 
      
      int doUnitTests();

      int testMLS_GradientCompactPoisson();
      int testMLS_LaplacianCompactPoisson();
      int testMLS_LaplacianMatrixCompactPoisson();
    };

  }
}

#endif

