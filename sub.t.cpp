#include <cassert>
#include "epsilon.h"

using namespace fms;

template<size_t N, class X = double>
int test_sub()
{
    epsilon<N, X> e0(X(0));
    e0 -= e0;
    for (size_t i = 0; i < N; ++i) {
        assert(e0[i] == X(0));
    }

    return 0;
}

static int t1 = test_sub<1, double>();
static int t2 = test_sub<2, double>();
static int t3 = test_sub<3, double>();
