#include <cassert>
#include "epsilon.h"

using namespace fms;

template<size_t N, class X = double>
int test_mul()
{
    epsilon<N, X> e0(X(0));
    auto e1 = e0;
    e1 *= e0;
    for (size_t i = 0; i < N; ++i) {
        assert(e1[i] == X(i == 2 ? 2 : 0));
    }
    e1 *= e0;
    for (size_t i = 0; i < N; ++i) {
        assert(e1[i] == X(i == 3 ? 6 : 0));
    }

    return 0;
}

static int t1 = test_mul<1, double>();
static int t2 = test_mul<2, double>();
static int t3 = test_mul<3, double>();

