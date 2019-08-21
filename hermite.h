// hermite.h - Hermite polynomials
#pragma once

namespace fms {

    template<class X>
    inline constexpr X Hermite(size_t n, const X& x)
    {
        return n == 0 ? X(1)
             : n == 1 ? x
                      : x * Hermite(n - 1, x) - (n - 1) * Hermite(n - 2, x);
    }

}
