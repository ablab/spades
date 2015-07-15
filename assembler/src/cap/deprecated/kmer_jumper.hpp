//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

namespace cap {

template <class TargetType>
class KmerJumper {
 public:
  virtual void SetTransition(char symbol, TargetType *link) const = 0;
  virtual const TargetType &GetTransition(char symbol) const = 0;
  virtual TargetType *GetTransitionLink(char symbol) const = 0;
  virtual bool HasTransition(char symbol) const = 0;
  virtual char Arity() const = 0;
};

template <class TargetType>
class MultiKmerJumper : public KmerJumper<TargetType> {
  TargetType *transitions[4];

 public:
  MultiKmerJumper() : transitions({NULL, NULL, NULL, NULL}) {
  }
  
  inline void SetTransition(char symbol, TargetType *link) {
    transitions[symbol] = link;
  }
  inline const TargetType &GetTransition(char symbol) const {
    return *(transitions[symbol]);
  }
  inline TargetType GetTransitionLink(char symbol) const {
    return transitions[symbol];
  }
  inline bool HasTransition(char symbol) const {
    return transitions[symbol] != NULL;
  }
  inline char Arity() const {
    return 4;
  }
};

template <class TargetType>
class SingleKmerJumper : public KmerJumper<TargetType> {
  TargetType *transition;

 public:
  SingleKmerJumper() : transition(NULL) {
  }

  inline void SetTransition(char symbol, TargetType *link) const {
    transition = link;
  }
  inline const TargetType &GetTransition(char symbol) const {
    return *transition;
  }
  inline TargetType *GetTransitionLink(char symbol) const {
    return transition;
  }
  inline bool HasTransition(char symbol) const {
    return transition != NULL;
  }
  inline char Arity() const {
    return 1;
  }

};

}
