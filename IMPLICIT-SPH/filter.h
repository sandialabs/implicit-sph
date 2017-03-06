#pragma once
#ifndef __FILTER_H__
#define __FILTER_H__

namespace LAMMPS_NS {

  // universal filter
  typedef class Filter FilterAny;

  class Filter {
  public:
  
    Filter() { }

    // eval itype only
    virtual bool yes(const int itype) const {
      return true;
    }
    virtual bool yes(const int itype, const int jtype) const {
      return true;
    }
    virtual bool same(const int itype, const int jtype) const { 
      return (itype == jtype);
    }
    virtual bool different(const int itype, const int jtype) const { 
      return (itype != jtype);
    }

    virtual ~Filter(){};
  };

  class FilterBinary : public FilterAny {
  protected:
    int _pair_yes[2];
  public:
    FilterBinary() : FilterAny() { 
      _pair_yes[0] =  127;
      _pair_yes[1] =  127;
    }

    void setPairYes(const int itype) {
      _pair_yes[0] = itype;
    }
    void setPairYes(const int itype, const int jtype) {
      _pair_yes[0] = itype;
      _pair_yes[1] = jtype;
    }

    virtual bool yes(const int itype) const {
      return (itype & _pair_yes[0]);
    }
    virtual bool yes(const int itype, const int jtype) const {
      return ((itype & _pair_yes[0]) && 
              (jtype & _pair_yes[1]));
    }
  };

  class FilterMatch : public FilterAny {
  protected:
    int _pair_yes[2];
  public:
    FilterMatch() : FilterAny() { 
      _pair_yes[0] =  127;
      _pair_yes[1] =  127;
    }

    void setPairYes(const int itype) {
      _pair_yes[0] = itype;
    }
    void setPairYes(const int itype, const int jtype) {
      _pair_yes[0] = itype;
      _pair_yes[1] = jtype;
    }

    virtual bool yes(const int itype) const {
      return (itype == _pair_yes[0]);
    }
    virtual bool yes(const int itype, const int jtype) const {
      return ((itype == _pair_yes[0]) && 
              (jtype == _pair_yes[1]));
    }
  };

  class FilterMatchBinary : public FilterAny {
  protected:
    int _pair_yes[2];
  public:
    FilterMatchBinary() : FilterAny() { 
      _pair_yes[0] =  127;
      _pair_yes[1] =  127;
    }

    void setPairYes(const int itype) {
      _pair_yes[0] = itype;
    }
    void setPairYes(const int itype, const int jtype) {
      _pair_yes[0] = itype;
      _pair_yes[1] = jtype;
    }

    virtual bool yes(const int itype) const {
      return (itype == _pair_yes[0]);
    }
    virtual bool yes(const int itype, const int jtype) const {
      return ((itype == _pair_yes[0]) && 
              (jtype &  _pair_yes[1]));
    }
  };

  class FilterExclusiveOr : public FilterBinary {
  public:
  
    FilterExclusiveOr() : FilterBinary() { }

    virtual bool yes(const int itype, const int jtype) const {
      return different(itype, jtype);
    }
  };

}

#endif
