//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef COMPARATOR_H_
#define COMPARATOR_H_

class Compare {
   public:
      bool operator() (std::string * lhs, std::string * rhs) const {
          return *lhs < *rhs;
      }
};

#endif /* COMPARATOR_H_ */
