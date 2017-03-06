#pragma once
#ifndef __FUNCTOR_BOUNDARY_ROBIN_H__
#define __FUNCTOR_BOUNDARY_ROBIN_H__

#include <iostream>
#include <string>
#include <vector>
#include "utils.h"
#include "functor.h"

namespace LAMMPS_NS {

  namespace Corrected {

    // this assume matrix is 3x3 (or 2x2) block matrix form in A_blk
    
    template<class PairIsph>
    class FunctorOuterBoundaryRobin : public FunctorOuter<PairIsph> {
    public:
      FunctorOuterBoundaryRobin(PairIsph *isph, 
                                // fill up necessary arguments 
                                )
        : FunctorOuter<PairIsph>(isph)
      { }
      
      void enterFor();
      void operator()(const int ii);
      void exitFor();
      
    protected:
      FUNCTOR_REMOVE_THIS_POINTER;
    };

    template<class PairIsph> inline void
    FunctorOuterBoundaryRobin<PairIsph>::enterFor() {
      // Copy necessary things
      // - basic helmholtz matrix is computed on A.crs
      // - deep or shallow copy of A.crs into A_blk(k1,k2)
      // - no need to change the same structure on the right hand side multivectors
      //   as the multivector is decomposed internally; just modify the vectors in operator()
    }
    
    template<class PairIsph> inline void
    FunctorOuterBoundaryRobin<PairIsph>::operator(const int ii) {
      // do modification
      // blocks are accessible by _pair->A_blk(k1,k2)
      // A00 A01 A02
      // A10 A11 A12
      // A20 A21 A22
      // each row and column represent x,y,z component
      // do whatever you want....

    }
    template<class PairIsph> inline void
    FunctorOuterBoundaryRobin<PairIsph>::exitFor() {
      // FillComplete for all necessary blocks
    }
  }
}
#endif

