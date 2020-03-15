#include "epsilon.h"
#include <cassert>

using namespace fms;

template <size_t N, class X = double> int test_op_minus()
{
    epsilon<N, X> x(1, 2);
    auto y = -x;
    auto x1 = epsilon<N, X>(1, 2);
    assert(x == x1);
    assert(y[0] == -x[0]);
    assert(y[1] == -x[1]);

    return 0;
}
int test_op_minus_ = test_op_minus<3>();

template <size_t N, class X = double> int test_sub()
{
    epsilon<N, X> e0(0, 1);
    e0 -= e0;
    for (size_t i = 0; i < N; ++i) {
        assert(e0[i] == X(0));
    }

    return 0;
}

static int t1 = test_sub<1, double>();
static int t2 = test_sub<2, double>();
static int t3 = test_sub<3, double>();
