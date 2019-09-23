// exp.h - exponential function
#pragma once
#include "epsilon.h"

namespace fms {

    template<class X>
    inline X exp(const X& x)
    {
        //!!! normalize x
        // Use frexp and ldrexp.
        // Requires log_e 2.
        X ex = X(1);
        X xn_(x); // x^n/n!

        size_t n = 1;
        while (xn_ + 1 != X(1)) {
            ex += xn_;
            xn_ *= x / ++n;
        }

        return ex;
    }

}
