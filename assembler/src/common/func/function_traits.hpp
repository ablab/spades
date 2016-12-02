#pragma once

#include <functional>

namespace func {

template<class F>
struct function_traits;

// function pointer
template<class R, class... Args>
struct function_traits<R(*)(Args...)> : public function_traits<R(Args...)> {
};

// member function pointer
template<class C, class R, class... Args>
struct function_traits<R(C::*)(Args...)> : public function_traits<R(C &, Args...)> {
};

// const member function pointer
template<class C, class R, class... Args>
struct function_traits<R(C::*)(Args...) const> : public function_traits<R(C &, Args...)> {
};

// member object pointer
template<class C, class R>
struct function_traits<R(C::*)> : public function_traits<R(C &)> {
};

template<class R, class... Args>
struct function_traits<R(Args...)> {
    using return_type = R;

    static constexpr std::size_t arity = sizeof...(Args);

    template<std::size_t N>
    struct arg {
        static_assert(N < arity, "invalid argument index");
        using type = typename std::tuple_element<N, std::tuple<Args...>>::type;
    };
};

template<class F>
struct function_traits<F &> : public function_traits<F> {
};

template<class F>
struct function_traits<F &&> : public function_traits<F> {
};

// functors & default implementation
template<class F>
struct function_traits {
private:
    using call_type = function_traits<decltype(&F::operator())>;

public:
    using return_type = typename call_type::return_type;

    // Remeber to get rid of this argument
    static constexpr std::size_t arity = call_type::arity - 1;

    template<std::size_t N>
    struct arg {
        static_assert(N < arity, "invalid argument index");
        // Remeber to get rid of this argument
        using type = typename call_type::template arg<N + 1>::type;
    };
};

} // namespace func
