#pragma once
#ifndef __TIME_BDF_H__
#define __TIME_BDF_H__

#include "utils.h"

namespace LAMMPS_NS {

  using namespace std;

  // solving du/dt + A u = f 
  // with the BDF scheme and unknown u_{n+1}
  // (gamma u_{n+1} - \sum alpha_{q} u_{n})/dt + A u_{n+1} = f_{n}

  class TimeBDF {
  public:
    
    // if LAMMPS restarts, nprev may not be zero
    TimeBDF(const int order)
      : _order(order), _nprev(0) { }
    virtual~TimeBDF() { }

    int getOrder() const;
    double getGamma() const;
    //double getCurrentTimestep() const;

    // setup new variables at the ntimestep of dt
    void updateTimestepBegin(const int ntimestep, 
                             const double dt);
    void updateTimestepEnd();

    // setting varialbes are placed between updateTimeBegin and End
    void setCurrentScalarVariable(const int npart,
                                  double *u, 
                                  double **uprev);
    void setCurrentVectorVariable(const int npart,
                                  const int dim,
                                  double **u, 
                                  double ***uprev);

    // setting varialbes are placed between updateTimeBegin and End
    void modifyCurrentScalarVariable(const char op,
                                     const int npart,
                                     double *u, 
                                     double **uprev);
    void modifyCurrentVectorVariable(const char op,
                                     const int npart,
                                     const int dim,
                                     double **u, 
                                     double ***uprev);

    // this is to compute position vector x
    void recoverRelativeVectorVariable(const int npart, 
                                       const int dim,
                                       double ***uprev);
    void updateRelativeVectorVariable(const int npart, 
                                      const int dim,
                                      double **unp1,
                                      double **un,
                                      double ***uprev);
    
    // extrapolation uhat = \sum beta uprev
    void extrapolateScalarVariable(const int npart,
                                   double **uprev,
                                   double *uhat);
    void extrapolateVectorVariable(const int npart,
                                   const int dim,
                                   double ***uprev,
                                   double **uhat);
    
    // diff udiff = \sum alpha uprev
    void diffScalarVariable(const int npart,
                            double **uprev,
                            double *udiff);
    void diffVectorVariable(const int npart,
                            const int dim,
                            double ***uprev,
                            double **udiff);
    
  protected:

    int _order, _nprev;

    double _dt[ISPH_BDF_MAX_ORDER], 
      _alpha[ISPH_BDF_MAX_ORDER], 
      _beta[ISPH_BDF_MAX_ORDER], 
      _gamma;
  };

  inline int 
  TimeBDF::getOrder() const {
    return (_order > _nprev ? _nprev : _order);
  }

  inline double
  TimeBDF::getGamma() const {
    return _gamma;
  }

  inline void
  TimeBDF::updateTimestepBegin(const int ntimestep,
                               const double dt) {
    // LAMMPS do initial run and repeat the same at ntimestep = 1
    _nprev = (ntimestep == 0 ? 0 : (ntimestep - 1));

    // shift previous solutions
    int nprev = min(_order-1,_nprev);

    // shift time 
    for (int k=0;k<nprev;++k) {
      int tgt = nprev - k;
      int src = tgt   - 1;

      // time
      _dt[tgt] = _dt[src];
    }

    // store the current solution in the first place
    _dt[0] = dt;
  }

  inline void
  TimeBDF::updateTimestepEnd() {
    // the current fields are stored between begin and end steps
    ++_nprev;

    // actual nprev is determined by min of nprev and history
    int order = getOrder();

    // compute gamma
    _gamma = 0.0;
    
    double cumsum_dt = 0.0, rho[ISPH_BDF_MAX_ORDER] = {};
    for (int i=0;i<order;++i) {
      cumsum_dt += _dt[i];
      rho[i] = _dt[0]/cumsum_dt;

      _gamma += rho[i];
    }

    // compute alpha and beta
    for (int i=0;i<order;++i) {
      double tmp = 1.0;
      for (int k=0;k<order;++k) 
        if (i != k) 
          tmp *= (1.0 - rho[k]/rho[i]);
      _beta[i] = 1.0/tmp;
      _alpha[i] = rho[i]/tmp;
    }
  }

  inline void
  TimeBDF::setCurrentScalarVariable(const int npart,
                                    double *u, 
                                    double **uprev) { 
    
    // shift previous solutions
    int nprev = min(_order-1,_nprev);

    // shift variables
    for (int k=0;k<nprev;++k) {
      int tgt = nprev - k;
      int src = tgt   - 1;

      for (int i=0;i<npart;++i) 
        uprev[i][tgt] = uprev[i][src];
    }

    // store the current solution in the first place
    for (int i=0;i<npart;++i)
      uprev[i][0] = u[i];
  }
  inline void
  TimeBDF::setCurrentVectorVariable(const int npart,
                                    const int dim, 
                                    double **u, 
                                    double ***uprev) { 
    
    // shift previous solutions
    int nprev = min(_order-1,_nprev);

    // shift variables
    for (int k=0;k<nprev;++k) {
      int tgt = nprev - k;
      int src = tgt   - 1;

      for (int i=0;i<npart;++i) 
        memcpy(uprev[i][tgt], uprev[i][src], sizeof(double)*dim); 
    }

    // store the current solution in the first place
    for (int i=0;i<npart;++i)
      memcpy(uprev[i][0], u[i], sizeof(double)*dim);
  }

  inline void
  TimeBDF::modifyCurrentScalarVariable(const char op,
                                       const int npart,
                                       double *u, 
                                       double **uprev) { 
    switch (op) {
    case '+': {
      for (int i=0;i<npart;++i)
        uprev[i][0] += u[i];
      break;
    }
    case '-': {
      for (int i=0;i<npart;++i)
        uprev[i][0] -= u[i];
      break;
    }
    case '=': {
      for (int i=0;i<npart;++i)
        uprev[i][0] = u[i];
      break;
    }
    }
  }
  inline void
  TimeBDF::modifyCurrentVectorVariable(const char op,
                                       const int npart,
                                       const int dim, 
                                       double **u, 
                                       double ***uprev) { 
    switch (op) {
    case '+': {
      for (int i=0;i<npart;++i)
        for (int k=0;k<dim;++k)
          uprev[i][0][k] += u[i][k];
      break;
    }
    case '-': {
      for (int i=0;i<npart;++i)
        for (int k=0;k<dim;++k)
          uprev[i][0][k] -= u[i][k];
      break;
    }
    case '=': {
      for (int i=0;i<npart;++i)
        for (int k=0;k<dim;++k)
          uprev[i][0][k] = u[i][k];
      break;
    }
    }
  }

  inline void
  TimeBDF::recoverRelativeVectorVariable(const int npart, 
                                        const int dim,
                                        double ***uprev) {
    int order = getOrder();
    for (int i=0;i<npart;++i) 
      for (int q=1;q<order;++q) 
        for (int k=0;k<dim;++k)
          uprev[i][q][k] = uprev[i][q-1][k] - uprev[i][q][k];
  }
  inline void
  TimeBDF::updateRelativeVectorVariable(const int npart, 
                                         const int dim,
                                         double **unp1,
                                         double **un,
                                         double ***uprev) {
    int order = getOrder();
    for (int i=0;i<npart;++i) 
      for (int q=(order-1);q>0;--q) 
        for (int k=0;k<dim;++k)
          uprev[i][q][k] = uprev[i][q-1][k] - uprev[i][q][k];

    for (int i=0;i<npart;++i) 
      for (int k=0;k<dim;++k)
        uprev[i][0][k] = unp1[i][k] - un[i][k];
  }

  inline void
  TimeBDF::extrapolateScalarVariable(const int npart,
                                     double **uprev,
                                     double *uhat) {
    int order = getOrder();
    for (int i=0;i<npart;++i) {
      uhat[i] = 0.0;
      for (int q=0;q<order;++q) 
        uhat[i] += _beta[q]*uprev[i][q];
    }
  }
  inline void
  TimeBDF::extrapolateVectorVariable(const int npart,
                                     const int dim,
                                     double ***uprev,
                                     double **uhat) {
    int order = getOrder();
    for (int i=0;i<npart;++i) {
      memset(uhat[i], 0, sizeof(double)*dim);
      for (int q=0;q<order;++q) 
        for (int k=0;k<dim;++k)
          uhat[i][k] += _beta[q]*uprev[i][q][k];
    }
  }

  inline void 
  TimeBDF::diffScalarVariable(const int npart,
                              double **uprev,
                              double *udiff) {
    int order = getOrder();
    for (int i=0;i<npart;++i) {
      udiff[i] = 0.0;
      for (int q=0;q<order;++q) 
        udiff[i] += _alpha[q]*uprev[i][q];
    }
  }
  inline void 
  TimeBDF::diffVectorVariable(const int npart,
                              const int dim,
                              double ***uprev,
                              double **udiff) {
    int order = getOrder();
    for (int i=0;i<npart;++i) {
      memset(udiff[i], 0, sizeof(double)*dim);
      for (int q=0;q<order;++q) 
        for (int k=0;k<dim;++k)
          udiff[i][k] += _alpha[q]*uprev[i][q][k];
    }
  }
  
}

#endif
