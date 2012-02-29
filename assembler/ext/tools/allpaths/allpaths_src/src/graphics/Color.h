// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research

// The purpose of this file is to define an rgb color class, and to define
// a bunch of colors in terms of it.

#ifndef COLOR_H
#define COLOR_H

#include <iostream>

class color {
  
 public:
  
  color( float r, float g, float b ) : r_(r), g_(g), b_(b) { }
  
  color( ) :
    r_(0.0),
    g_(0.0),
    b_(0.0)
    { }
  
  float R( ) const { return r_; }
  float G( ) const { return g_; }
  float B( ) const { return b_; }
  
  friend bool operator==( const color& lhs, const color& rhs )
    { return ( lhs.r_ == rhs.r_ && lhs.g_ == rhs.g_ && lhs.b_ == rhs.b_); }
  
  friend bool operator!=( const color& lhs, const color& rhs )
    { return !( lhs == rhs ); }
  
  friend std::ostream& operator<<( std::ostream& o, const color& c);
  
 private:
  
  float r_, g_, b_;
};

const color white  (1, 1, 1), black(0, 0, 0);
const color red    (1, 0, 0), darkred  (0.5, 0, 0);
const color green  (0, 1, 0), darkgreen(0, 0.5, 0);
const color blue   (0, 0, 1), darkblue (0, 0, 0.5);

const color magenta(1, 0, 1), yellow(1, 1, 0), cyan(0, 1, 1);
const color orange(1, 0.5, 0);
const color brown(0.7, 0.5, 0);
const color purple(0.5, 0, 0.5);

const color pink(1, 0.5, 0.5);
const color lighterpink(1, 0.8, 0.8);

const color gray(0.5, 0.5, 0.5);
const color darkgray(0.4, 0.4, 0.4);
const color darkergray(0.3, 0.3, 0.3);
const color lightgray(0.6, 0.6, 0.6);
const color lightergray(0.8, 0.8, 0.8);

const color azure(0, 0.5, 1);
const color orangered(1, 0.27, 0);
const color seashell(1, 0.96, 0.93);
const color limegreen(0.2, 0.8, 0.2);
const color wheat(0.96, 0.87, 0.7);



const color & MakeUpColor(int num);


#endif
