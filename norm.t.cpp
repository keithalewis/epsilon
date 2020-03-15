// norm.t.cpp - Test epsilon norm
#include "epsilon.h"
#include <cassert>

using namespace fms;

template <size_t N, class X> int test_norm()
{
    {
        epsilon<N, X> e;
        assert(0 == e.norm());
    }
    {
        X x = 1;
        epsilon<N, X> e(x, 1);
        assert(x + X(1) * (N > 1) == e.norm());
        assert(sqrt(x * x + X(1) * (N > 1)) == e.norm(2));
    }

    return 0;
}
int test_norm_double1 = test_norm<1, double>();
int test_norm_double2 = test_norm<2, double>();
int test_norm_double3 = test_norm<3, double>();
