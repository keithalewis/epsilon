// hermite.h - Hermite polynomials
#pragma once

namespace fms {

    template<class X>
    inline constexpr X Hermite(size_t n, const X& x)
    {
        if (n == 0)
            return 1;
        else if (n == 1)
            return x;
        else
            return x*Hermite(n - 1, x) - (n - 1)*Hermite(n - 2, x); 
    }

}
