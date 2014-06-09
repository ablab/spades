//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include "standard_base.hpp"
#include "boost/function.hpp"

namespace func {

//to use with boost::function-s
template<class T>
void Composition(T t, boost::function<void(T)> f1,
		boost::function<void(T)> f2) {
	if (f1)
		f1(t);
	if (f2)
		f2(t);
}

template<class A, class B>
class Func {
public:
	typedef function<B(A)> function_t;

	virtual B Apply(A a) const = 0;

	virtual ~Func() {
	}
};

template<class T>
class AndOperator;

template<class T>
class OrOperator;

template<class T>
class NotOperator;

template<class T>
class Predicate: public Func<T, bool> {
public:
	bool Apply(T t) const {
		return Check(t);
	}

	virtual bool Check(T t) const = 0;

	virtual ~Predicate() {
	}
};

template<class T>
class AdaptorPredicate: public Predicate<T> {
    typedef boost::function<bool(T)> pred_func_t;
    pred_func_t pred_f_;
public:
    AdaptorPredicate(pred_func_t pred_f) :
        pred_f_(pred_f) {
    }

    bool Check(T t) const {
        return pred_f_(t);
    }
};

//template<class T>
//const shared_ptr<Predicate<T>> operator &&(const shared_ptr<Predicate<T>>& a, const shared_ptr<Predicate<T>>& b) {
//	return AndOperator<T>(a, b);
//}
//
//template<class T>
//const shared_ptr<Predicate<T>> operator ||(const shared_ptr<Predicate<T>>& a, const shared_ptr<Predicate<T>>& b) {
//	return OrOperator<T>(a, b);
//}
//
//template<class T>
//const shared_ptr<Predicate<T>> operator !(const shared_ptr<Predicate<T>>& a) {
//	return NotOperator<T>(a);
//}

template<class T>
const shared_ptr<Predicate<T>> And(const shared_ptr<Predicate<T>>& a,
		const shared_ptr<Predicate<T>>& b) {
	return make_shared<AndOperator<T>>(a, b);
}

template<class T>
const shared_ptr<Predicate<T>> Or(const shared_ptr<Predicate<T>>& a,
		const shared_ptr<Predicate<T>>& b) {
	return make_shared<OrOperator<T>>(a, b);
}

template<class T>
const shared_ptr<Predicate<T>> Not(const shared_ptr<Predicate<T>>& a) {
	return make_shared<NotOperator<T>>(a);
}

template<class T>
class AlwaysTrue: public Predicate<T> {
public:

	bool Check(T /*t*/) const {
		return true;
	}

};

template<class T>
class AlwaysFalse: public Predicate<T> {
public:

	bool Check(T /*t*/) const {
		return false;
	}

};

template<class T>
class NotOperator: public Predicate<T> {
	shared_ptr<Predicate<T>> a_;

public:
	NotOperator(const shared_ptr<Predicate<T>>& a) :
			a_(a) {
	}

	bool Check(T t) const {
		return !a_->Check(t);
	}
};

template<class T>
class AndOperator: public Predicate<T> {
	shared_ptr<Predicate<T>> a_;
	shared_ptr<Predicate<T>> b_;

public:
	AndOperator(const shared_ptr<Predicate<T>>& a,
			const shared_ptr<Predicate<T>>& b) :
			a_(a), b_(b) {
	}

	bool Check(T t) const {
		return a_->Check(t) && b_->Check(t);
	}
};

template<class T>
class OrOperator: public Predicate<T> {
	shared_ptr<Predicate<T>> a_;
	shared_ptr<Predicate<T>> b_;

public:
	OrOperator(const shared_ptr<Predicate<T>>& a,
			const shared_ptr<Predicate<T>>& b) :
			a_(a), b_(b) {
	}

	bool Check(T t) const {
		return a_->Check(t) || b_->Check(t);
	}
};

}
