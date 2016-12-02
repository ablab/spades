#pragma once

#include "function_traits.hpp"

#include <memory>
#include <functional>

namespace func {

template<class T>
class AbstractPredicate {
public:
    typedef T checked_type;

    virtual bool Check(T t) const = 0;

    bool operator()(T t) const { return Check(t); }

    virtual ~AbstractPredicate() {}
};

template<typename T>
class TypedPredicate {
    struct TypedPredicateConcept {
        virtual ~TypedPredicateConcept() { };

        virtual bool operator()(T x) const = 0;
    };

    template<class P>
    struct TypedPredicateModel : TypedPredicateConcept {
        TypedPredicateModel(P p)
                : data_(std::move(p)) { }

        virtual bool operator()(T x) const override {
            return data_(x);
        }

        P data_;
    };

    std::shared_ptr<const TypedPredicateConcept> self_;

public:
    typedef T checked_type;

    template<typename P>
    TypedPredicate(P p)
            : self_(std::make_shared<TypedPredicateModel<P>>(std::move(p))) { }

    bool operator()(T x) const {
        return self_->operator()(x);
    }
};

template<typename T>
class AlwaysTrueOperator {
public:
    typedef T checked_type;

    bool operator()(T) const {
        return true;
    }
};

template<typename T>
class AlwaysFalseOperator {
    typedef T checked_type;

public:
    bool operator()(T) const {
        return false;
    }
};

template<typename T>
class AndOperator {
public:
    typedef T checked_type;

    AndOperator(TypedPredicate<T> lhs, TypedPredicate<T> rhs)
            : lhs_(std::move(lhs)),
              rhs_(std::move(rhs)) { }

    bool operator()(T x) const {
        return lhs_(x) && rhs_(x);
    }

private:
    const TypedPredicate<T> lhs_, rhs_;
};

template<typename T>
class OrOperator {
public:
    typedef T checked_type;

    OrOperator(TypedPredicate<T> lhs, TypedPredicate<T> rhs)
            : lhs_(std::move(lhs)), rhs_(std::move(rhs)) { }

    bool operator()(T x) const {
        return lhs_(x) || rhs_(x);
    }

private:
    const TypedPredicate<T> lhs_, rhs_;
};

template<typename T>
class NotOperator {
public:
    typedef T checked_type;

    NotOperator(const TypedPredicate<T> p)
            : p_(std::move(p)) { }

    bool operator()(T x) const {
        return !p_(x);
    }

private:
    const TypedPredicate<T> p_;
};

template<class P,
        bool = function_traits<P>::arity == 1 &&
               std::is_same<typename function_traits<P>::return_type, bool>::value>
struct is_predicate : public std::true_type {
};

template<class P>
struct is_predicate<P, false> : public std::false_type {
};

template<class TP1, class TP2,
        typename _T1 = typename function_traits<TP1>::template arg<0>::type,
        typename _T2 = typename function_traits<TP2>::template arg<0>::type,
        typename =
        typename std::enable_if<std::is_same<_T1, _T2>::value &&
                                is_predicate<TP1>::value && is_predicate<TP2>::value
        >::type>
TypedPredicate<_T1> And(TP1 lhs, TP2 rhs) {
    return AndOperator<_T1>(lhs, rhs);
}

template<class TP1, class TP2,
        typename _T1 = typename function_traits<TP1>::template arg<0>::type,
        typename _T2 = typename function_traits<TP2>::template arg<0>::type,
        typename =
        typename std::enable_if<std::is_same<_T1, _T2>::value &&
                                is_predicate<TP1>::value && is_predicate<TP2>::value
        >::type>
TypedPredicate<_T1> Or(TP1 lhs, TP2 rhs) {
    return OrOperator<_T1>(lhs, rhs);
}

template<class TP,
        typename _T = typename function_traits<TP>::template arg<0>::type,
        typename =
        typename std::enable_if<is_predicate<TP>::value>::type>
TypedPredicate<_T> Not(TP p) {
    return NotOperator<_T>(p);
}

template<class T>
TypedPredicate<T> AlwaysTrue() {
    return AlwaysTrueOperator<T>();
}

template<class T>
TypedPredicate<T> AlwaysFalse() {
    return AlwaysFalseOperator<T>();
}

} // namespace func
