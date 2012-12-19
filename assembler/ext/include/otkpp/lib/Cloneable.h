
#ifndef CLONEABLE_H

/// Implements the base class for all "cloneable" classes.
template< typename D, typename B > class Cloneable : public B
{
 public:
  virtual B *clone() const
  {
    return new D(*static_cast< const D * >(this));
  }
};

#define CLONEABLE_H

#endif
