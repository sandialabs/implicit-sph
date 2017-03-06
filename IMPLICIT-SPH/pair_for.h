#pragma once
#ifndef __PAIR_FOR_H__
#define __PAIR_FOR_H__


namespace LAMMPS_NS {

  template<class Functor> inline void
  PairFor(Functor &functor, size_t nwork) {
    functor.enterFor();
    for (size_t iter=0;iter<nwork;++iter)
      functor(iter);
    functor.exitFor();
  }

}

#endif

