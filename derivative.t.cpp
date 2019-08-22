#include <cassert>
#include "epsilon.h"

using namespace fms;

template<class X>
X p(const X& x)
{
    return 1. + 2. * x + 3. * x * x;
}
template<class X>
X dp(const X& x)
{
    return 2. + 6. * x;
}
template<class X>
X ddp([[maybe_unused]] const X& x)
{
    return 6.;
}

int test_derivative()
{
    {
        epsilon<5, double> e{ 1, 1, 1, 1, 1 };
        e *= e;
        assert(e[0] == 1);
        assert(e[1] == 2);
        assert(e[2] == 4);
        assert(e[3] == 8);
        assert(e[4] == 16);

    }
    {
        epsilon<5, double> e(0.);
        auto e1 = e;
        e1 *= e; assert(e1 == epsilon<5>({ 0, 0, 2, 0, 0 }));
        e1 *= e; assert(e1 == epsilon<5>({ 0, 0, 0, 6, 0 }));
        e1 *= e; assert(e1 == epsilon<5>({ 0, 0, 0, 0, 24 }));
        e1 *= e; assert(e1 == epsilon<5>({ 0, 0, 0, 0, 0 }));
    }
    {
        double x = 1;
        epsilon<3, double> e(x);
        auto pe = p(e);
        assert(pe[0] == p(x));
        assert(pe[1] == dp(x));
        assert(pe[2] == ddp(x));
    }

    return 0;
}

static int derivative = test_derivative();