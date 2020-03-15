// pow.h - power function
#pragma once
#include "epsilon.h"

namespace fms {

template <size_t N, class X, typename IsArithmetic<X>>
inline epsilon<N, X> pow(const epsilon<N, X> &x, int n)
{
    epsilon<N, X> xn;

    xn[0] = 1;
    if (n < 0) {
        x = 1 / x;
        n = -n;
    }
    while (n--) {
        xn *= x;
    }

    return xn;
}
} // namespace fms
