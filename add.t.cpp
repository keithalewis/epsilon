#include <cassert>
#include "epsilon.h"

using namespace fms;


template<size_t N, class X = double>
int test_add()
{
    epsilon<N, X> e0;
    e0 += e0;
    for (size_t i = 0; i < N; ++i) {
        assert(e0[i] == X(i == 1 ? 2 : 0));
    }

    return 0;
}

static int t1 = test_add<1, double>();
static int t2 = test_add<2, double>();
static int t3 = test_add<3, double>();
