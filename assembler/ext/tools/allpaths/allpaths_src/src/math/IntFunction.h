///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
// author: Filipe Ribeiro      03/2011
// 
//  A function defined in the integer interval [x_min, x_max]. 
//
//   It's your job to make sure you are within bounds.
//
//          


#ifndef _MATH__INT_FUNCTION_H
#define _MATH__INT_FUNCTION_H

#include <deque>

#include "MainTools.h"


#define __INT_FUNCTION_BINARY_VERSION__ 2


template<class T>
class IntFunction
{
  bool _neg_inf;
  bool _pos_inf;
  int _x0;
  deque<T> _v;

public:
  IntFunction(const int x0 = 0, 
              const int x1 = 0,
              const T v = 0)
    : _neg_inf(false), 
      _pos_inf(false),
      _x0(x0),
      _v((x1 >= x0) ? x1 - x0 + 1 : 0, v)
  {
    ForceAssertGe(x1, x0);
  }

  ~IntFunction()
  {
    // cout << "~int_func(): " << x_min() << ", " << x_max()
    //     << " _neg.sz(): " << _negative.size() 
    //     << " _pos_sz(): " << _positive.size() << endl;
  }
  
  int size() const { return _v.size(); }

  int x_min() const { return _x0; }
  int x_max() const { return _x0 + size() - 1; }
  T f_x_min() const { return _v.front(); }
  T f_x_max() const { return _v.back(); }

  int x_f_min() const;
  int x_f_max() const;
  T f_min() const { return (*this)[x_f_min()]; }
  T f_max() const { return (*this)[x_f_max()]; }

  T sum_in(const int x0, const int x1) const 
  {
    T sum = 0;
    for (int x = x0; x < x1; x++) sum += (*this)[x];
    return sum;
  }
  T sum_below(const int x1) const { return sum_in(x_min(), x1); }
  T sum_above(const int x0) const { return sum_in(x0, x_max()); }
  T sum() const { return sum_in(x_min(), x_max()); }
  

  void expand_neg_infinity() { _neg_inf = true; }
  void expand_pos_infinity() { _pos_inf = true; }
  void expand_infinity() { _neg_inf = _pos_inf = true; }

  T operator [] (const int x) const
  {
    if (x < x_min()) return (_neg_inf) ? f_x_min() : 0;
    if (x > x_max()) return (_pos_inf) ? f_x_max() : 0;
    return _v[x - _x0];
  }

  T & operator [] (const int x)
  {
    if (x > x_max()) {
      const T f1 = (_pos_inf) ? f_x_max() : 0;
      _v.insert(_v.end(), x - x_max(), f1);
    }
    if (x < x_min()) { 
      const T f0 = (_neg_inf) ? f_x_min() : 0;
      const int n = x_min() - x;
      for (int i = 0; i != n; i++) _v.push_front(f0);
      //_v.insert(_v.rend(), x_min() - x, f0);   
      _x0 = x;
    }
    return _v[x - _x0];
  }
  
  IntFunction & operator += (const IntFunction & f)
  {
    for (int x = f.x_min(); x <= f.x_max(); x++) (*this)[x] += f[x];
    return *this;
  }

  virtual void to_text_file(const String & fn) const;


public:
  // ---- SELF_SERIALIZABLE method
  size_t writeBinary(BinaryWriter& writer) const;

  // ---- SELF_SERIALIZABLE method
  void readBinary(BinaryReader& reader);

  // ---- SELF_SERIALIZABLE method
  static size_t externalSizeof() { return 0; }
};


SELF_SERIALIZABLE(IntFunction<int64_t>);
SELF_SERIALIZABLE(IntFunction<int32_t>);
SELF_SERIALIZABLE(IntFunction<int16_t>);
SELF_SERIALIZABLE(IntFunction<int8_t>);
SELF_SERIALIZABLE(IntFunction<uint64_t>);
SELF_SERIALIZABLE(IntFunction<uint32_t>);
SELF_SERIALIZABLE(IntFunction<uint16_t>);
SELF_SERIALIZABLE(IntFunction<uint8_t>);
SELF_SERIALIZABLE(IntFunction<float>);
SELF_SERIALIZABLE(IntFunction<double>);




template<class T>
int IntFunction<T>::x_f_max() const
{
  int x_f_max = x_min();
  T f_max = (*this)[x_f_max];

  const int x1 = x_max();
  for (int x = x_f_max + 1; x <= x1; x++) {
    const T f = (*this)[x];
    if (f > f_max) {
      x_f_max = x;
      f_max = f;
    }
  }
  return x_f_max;
}


template<class T>
int IntFunction<T>::x_f_min() const
{
  int x_f_min = x_min();
  T f_min = (*this)[x_f_min];

  const int x1 = x_max();
  for (int x = x_f_min + 1; x <= x1; x++) {
    const T f = (*this)[x];
    if (f < f_min) {
      x_f_min = x;
      f_min = f;
    }
  }
  return x_f_min;
}



template<class T>
void IntFunction<T>::to_text_file(const String & fn) const
{
  const int x0 = x_min();
  const int x1 = x_max();

  ofstream os;
  os.open(fn.c_str());
  os << "# x_min = " << x0 << endl;
  os << "# x_max = " << x1 << endl;
  os << "# 1:x  2:f_x" << endl;
  os << fixed;

  for (int x = x0; x <= x1; x++) {
    os << setw(10) << x << " "
       << setw(16) << (*this)[x]
       << endl;
  }
  os.close();
}




// ---- SELF_SERIALIZABLE method
template<class T>
size_t IntFunction<T>::writeBinary(BinaryWriter& writer) const
{
  size_t len = 0;

  const int version = __INT_FUNCTION_BINARY_VERSION__;
  len += writer.write(version);

  len += writer.write(_neg_inf);
  len += writer.write(_pos_inf);
  len += writer.write(_x0);
  const size_t n = _v.size();
  len += writer.write(n);
  len += writer.writeItr(_v.begin(), _v.end());

  return len;
}


// ---- SELF_SERIALIZABLE method
template<class T>
void IntFunction<T>::readBinary(BinaryReader& reader)
{ 
  int version;
  reader.read(&version);
  if (version != __INT_FUNCTION_BINARY_VERSION__) {
    cout << Date() << "**** binary data on disk is from a different code version." << endl;
    ForceAssertEq(version, __INT_FUNCTION_BINARY_VERSION__);
  }

  reader.read(&_neg_inf);
  reader.read(&_pos_inf);
  reader.read(&_x0);
  size_t n; 
  reader.read(&n); 
  _v.resize(n, 0);
  reader.readItr(_v.begin(), _v.end());
}









template<class T>
class IntFunctionPrimitive
{
  const IntFunction<T> & _f;
  IntFunction<T>         _f_sum;

public:
  IntFunctionPrimitive(const IntFunction<T> & f)
    : _f(f),
      _f_sum(f.x_min(), f.x_max())
  {
    const int x0 = _f.x_min();
    const int x1 = _f.x_max();

    _f_sum[x0] = _f[x0];
    for (int x = x0 + 1; x <= x1; x++)
      _f_sum[x] = _f_sum[x - 1] + _f[x];
  }

  T f_sum(const int a, const int b) const;

};




template<class T>
T IntFunctionPrimitive<T>::f_sum(const int a, const int b) const
{
  if (a == b) return _f[a];
  
  const int x0 = _f.x_min();
  const int x1 = _f.x_max();
  const T f0 = _f.f_x_min();
  const T f1 = _f.f_x_max();

    
  const int aa = (a < b) ? a - 1 : b;
  const int bb = (a < b) ? b     : a - 1;

  const int n1a = (aa > x1) ? aa - x1 : 0;
  const int n1b = (bb > x1) ? bb - x1 : 0;

  const int n0a = (aa < x0) ? aa - x0 : 0;
  const int n0b = (bb < x0) ? bb - x0 : 0;

  const T sum_a = (aa >= x1) ? _f_sum.f_x_max() : (aa <= x0) ? _f_sum.f_x_min() : _f_sum[aa];
  const T sum_b = (bb >= x1) ? _f_sum.f_x_max() : (bb <= x0) ? _f_sum.f_x_min() : _f_sum[bb];


  return (sum_b - sum_a) + f1 * (n1b - n1a) + f0 * (n0b - n0a);
}









#endif
