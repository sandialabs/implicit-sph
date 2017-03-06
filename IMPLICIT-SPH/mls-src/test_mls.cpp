#include <iostream>
#include <string>
#include "utils.h"

#include "pair_isph_mls.h"
#include "test_mls.h"

namespace LAMMPS_NS {

  namespace MLS {

    using namespace std;

    int TestSuite::doUnitTests() {
      int r_val = 0, cnt = 0;

      __UNIT_TEST__(testMLS_GradientCompactPoisson(),        r_val, cnt);
      __UNIT_TEST__(testMLS_LaplacianCompactPoisson(),       r_val, cnt);
      __UNIT_TEST__(testMLS_LaplacianMatrixCompactPoisson(), r_val, cnt);

      if (_pair->comm->me == 0) {
        cout << "MLS::TestSuite::"
             << "  # of tests = " << cnt 
             << ", # of failures = " << r_val
             << endl;
      }

      return r_val;
    }

  }
}

