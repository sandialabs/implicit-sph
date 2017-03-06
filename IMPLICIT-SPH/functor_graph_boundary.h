#pragma once
#ifndef __FUNCTOR_GRAPH_BOUNDARY_H__
#define __FUNCTOR_GRAPH_BOUNDARY_H__

#include <iostream>
#include <string>
#include "utils.h"
#include "functor.h"

namespace LAMMPS_NS {

  // create a graph with a static profile and set it up
  template<class PairIsph>
  class FunctorOuterGraphBoundary : public FunctorOuter<PairIsph> {
  public:
    FunctorOuterGraphBoundary(PairIsph *isph, Epetra_Map *map, 
                              double **normal,
                              Epetra_CrsGraph **graph)
      : FunctorOuter<PairIsph>(isph),
        _map(map),
        _normal(normal),
        _graph(graph) { }
    void enterFor();
    void operator()(const int ii);
    void exitFor();

  protected:
    FUNCTOR_REMOVE_THIS_POINTER;

    // input 
    Epetra_Map *_map;

    double **_normal;
    
    // outputs that will be created
    Epetra_CrsGraph **_graph;

    // internal utility objects; stores global tags
    Epetra_IntSerialDenseVector _row, _idx;
  };
  template<class PairIsph> inline void 
  FunctorOuterGraphBoundary<PairIsph>::enterFor() {

    // at this moment, graph construction does not allow filter which means 
    // each row members in the graph represents all effective neighbors regardless
    // of their types
    if (_filter != NULL)
      _pair->error->all(FLERR, "FunctorGraphBoundary:: Filter is not allowed"); 

    // count graph structure 
    _row.Size(_inum);
    for (int ii=0;ii < _inum;++ii) {
      bool has_normal = util.dotVectors(_dim, _normal[ii], _normal[ii]) > 0.5;
      if (has_normal)
        _row[ii] = _numneigh[_ilist[ii]] + 1;
      else 
        _row[ii] = 0;
    }
    _idx.Size(_row.InfNorm());

    // create a static graph
    (*_graph) = new Epetra_CrsGraph(Copy, *_map, _row.Values(), true);
  }

  template<class PairIsph> inline void 
  FunctorOuterGraphBoundary<PairIsph>::operator()(const int ii) {
    int i = _ilist[ii];
    int itype = _type[i];

    // neighbor around i
    int *jlist = _firstneigh[i];
    int jnum = _numneigh[i];

    // # of columns counter
    size_t cnt = 0;
    int *tag_in_boundary = _idx.Values();

    // boundary interface depends on the existence of normal vectors
    bool has_normal = util.dotVectors(_dim, _normal[i], _normal[i]) > 0.5;
    if (has_normal) {
      for (int jj=0;jj<jnum;++jj) {
        // j th particle; always put neighmask two bits are reserved for LAMMPS
        int j = (jlist[jj] & NEIGHMASK);
        int jtype = _type[j];
        
        // compute distance between partile i and j
        double rsq = 0.0;
        for (int k=0;k<_dim;++k) {
          double r = _x[i][k] - _x[j][k];
          rsq += (r*r);
        }
        
        // integrate
        if (rsq < LOOKUP(_pair->cutsq, itype, jtype)) 
          tag_in_boundary[cnt++] = _tag[j];
      }
      tag_in_boundary[cnt++] = _tag[i];

      // add the connectivity for ith particle; this allows duplication
      (*_graph)->InsertGlobalIndices(_tag[i], (int)cnt, tag_in_boundary);
    }
  }

  template<class PairIsph> inline void 
  FunctorOuterGraphBoundary<PairIsph>::exitFor() {
    // finalize connectivity construction; here redundant entrees are removed
    (*_graph)->FillComplete();
    (*_graph)->OptimizeStorage();
  }
}

#endif

