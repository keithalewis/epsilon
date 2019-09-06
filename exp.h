// exp.h - exponential function
#pragma once
#include "pow.h"

namespace fms {

    // exp(x epsilon^k)
    template<size_t N, class X, class IsArithmetic<X>>
    inline epsilon<N, X> expk(size_t k, const X& x)
    {
        size_t n_ = 1; // n factorial
        epsilon<N, X> e, xk = pow(x, k);
        e[0] = 1;

        for (size_t n = 1; n <= k && n*k < N; ++n) {
            n_ *= n;
            e += xn / n_;
            xn *= x;
        }

        return e;
    }
    template<size_t N, class X, class IsArithmetic<X>>
    inline epsilon<N, X> exp(const epsilon<N, X>& x)
    {
        return x;
    }
}
