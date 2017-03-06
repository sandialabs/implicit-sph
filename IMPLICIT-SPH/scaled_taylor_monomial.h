#pragma once
#ifndef __SCALED_TAYLOR_MONOMIAL_H__
#define __SCALED_TAYLOR_MONOMIAL_H__

#include <iostream>
#include "math.h"
#include "float.h"
#include "utils.h"

namespace LAMMPS_NS {

  using namespace std;

  class ScaledTaylorMonomial {
  protected:
    unsigned int _dim;
    double _fact[ISPH_MLS_MAX_ORDER+1];
    bool _is_interpolation_enabled;

  public:

    ScaledTaylorMonomial(const unsigned int dim, 
                         const bool is_interpolation_enabled)
      : _dim(dim), 
        _is_interpolation_enabled(is_interpolation_enabled) {
      util.factorial(ISPH_MLS_MAX_ORDER, _fact);       // factorials upto maxp
    }

    int ndof(const unsigned int np) {      
      int r_val = 0;
      switch (_dim) {
      case 2: r_val = (np+1)*(np+2)/2;        break;
      case 3: r_val = (np+1)*(np+2)*(np+3)/6; break;
      }
      return (r_val - _is_interpolation_enabled);
    }
    
    void show_dof_order(const unsigned int np, ostream &os) {
      int idx = 0;

      const int is_dim3 = (_dim == 3);
      for     (int k3=0;k3<=((int)np*is_dim3);++k3)          // z
        for   (int k2=0;k2<=((int)np-k3);     ++k2)          // y
          for (int k1=0;k1<=((int)np-k2-k3);  ++k1)  {       // x
            // skip the constant dof if interpolation property is enabled
            if (_is_interpolation_enabled && !(k1+k2+k3))
              continue;
            os << "idx = " << idx++ << ", order = " << k1 << ", " << k2 << ", " << k3 << endl;
          }
      IS_VALID(ndof(np) == idx, 0, ">> Size does not match");
    }


    void val(const unsigned int np, 
             const double *rij, 
             const double rth, 
             double *u) {
      // scale the distance 
      // - note that the input r_ij = x_i - x_j
      double sji[3];
      for (unsigned int k=0;k<_dim;++k)
        sji[k] = -rij[k]/rth;
      
      // loop over the monomials
      int idx = 0;

      const int is_dim3 = (_dim == 3);
      for     (int k3=0;k3<=((int)np*is_dim3);++k3)          // z
        for   (int k2=0;k2<=((int)np-k3);     ++k2)          // y
          for (int k1=0;k1<=((int)np-k2-k3);  ++k1)  {       // x
            // skip the constant dof if interpolation property is enabled
            if (_is_interpolation_enabled && !(k1+k2+k3))
              continue;

            // compute scaled monomial function
            u[idx++] = (pow(sji[0], k1)/_fact[k1] *
                        pow(sji[1], k2)/_fact[k2] * 
                        pow(sji[2], k3)/_fact[k3]);
          }
      IS_VALID(ndof(np) == idx, 0, ">> Size does not match");
    }

    // evaluate du at x = x_i (which is r = 0)
    // - compute non-zeo du and its location
    void dval(const unsigned int np, 
              const double rth, 
              const unsigned int p1,
              const unsigned int p2,
              const unsigned int p3,
              int &idx, double &du) {
      const unsigned int pp = p1 + p2 + p3;
      
      // loop over the monomials
      idx = 0; 
      du  = 0.0;

      const int is_dim3 = (_dim == 3);
      for     (int k3=0;k3<=((int)np*is_dim3);++k3)            // z
        for   (int k2=0;k2<=((int)np-k3);     ++k2)            // y
          for (int k1=0;k1<=((int)np-k2-k3);  ++k1,++idx)      // x
            if (k1 == (int)p1 && k2 == (int)p2 && k3 == (int)p3) {
              // du = 1/eps^|pp|
              du = 1.0/pow(rth, pp);

              // when the interpolation property is enabled,
              // the index location shifted since the const dof 
              // is removed
              idx -= _is_interpolation_enabled;
              return;
            }
    }

    // evaluate du at x = x_j (which r != 0)
    void dval(const unsigned int np,
              const double *rij, 
              const double rth,
              const unsigned int p1, 
              const unsigned int p2,
              const unsigned int p3,
              const double scale,
              double *du,
              const bool update = false) {
      // scale the distance 
      // - note that the input r_ij = x_i - x_j
      double sji[3];
      for (unsigned int k=0;k<_dim;++k)
        sji[k] = -rij[k]/rth;
      
      // compute eps^beta = (rth^p1, rth^p2, rth^p3)
      const double 
        rth_p1 = pow(rth, p1), 
        rth_p2 = pow(rth, p2), 
        rth_p3 = pow(rth, p3);
      
      // initialize if computing mode is not update
      if (!update)
        memset(du, 0, ndof(np)*sizeof(double));

      // loop over the monomials
      int idx = 0;

      const int is_dim3 = (_dim == 3);
      for     (int k3=0;k3<=((int)np*is_dim3);++k3) {  const int dp3 = (k3 - p3);    // z
        for   (int k2=0;k2<=((int)np-k3);     ++k2) {  const int dp2 = (k2 - p2);    // y
          for (int k1=0;k1<=((int)np-k2-k3);  ++k1) {  const int dp1 = (k1 - p1);    // x
            // skip the constant dof if interpolation property is enabled
            if (_is_interpolation_enabled && !(k1+k2+k3))
              continue;
            
            // compute derivatives of scaled monomial function
            // - du = 1/(alpha - beta)!*(sji^{alpha - beta)/eps^beta
            if (dp1 >= 0 && dp2 >= 0 && dp3 >= 0) 
              du[idx++] += scale*(pow(sji[0], dp1)/rth_p1/_fact[dp1] *
                                  pow(sji[1], dp2)/rth_p2/_fact[dp2] *
                                  pow(sji[2], dp3)/rth_p3/_fact[dp3]);
            else 
              ++idx; 
          }
        }
      }
      IS_VALID(ndof(np) == idx, 0, ">> Size does not match");
    }
  };

}

#endif
