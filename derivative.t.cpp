#include <cassert>
#include "epsilon.h"

using namespace fms;

template<class X>
X p(const X& x)
{
    return 1 + 2 * x + 3 * x * x;
}
template<class X>
X dp(const X& x)
{
    return 2 + 6 * x;
}
template<class X>
X ddp([[maybe_unused]] const X& x)
{
    return 6;
}

int test_derivative()
{
    {
        epsilon<5, double> e0;
        auto e = 1.0 + e0 + e0*e0 + e0*e0*e0 + e0*e0*e0*e0;
        e *= e;
        assert(e[0] == 1);
        assert(e[1] == 2);
        assert(e[2] == 6);
        assert(e[3] == 24);
        assert(e[4] == 120);

    }
    {
        double x = 2;
        auto x2 = epsilon<4>(x) * epsilon<4>(x);
        assert(x2[0] == x * x);
        assert(x2[1] == 2 * x);
        assert(x2[2] == 2);
        assert(x2[3] == 0);
    }
 
    return 0;
}

static int derivative = test_derivative();