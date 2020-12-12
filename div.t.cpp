#include <cassert>
#include "epsilon.h"

using namespace fms;

template<size_t N, class X>
int test_div()
{
    epsilon<N,X> e1(1);
    auto e2 = e1 / e1;
    assert(e2 == X(1));

    return 0;
}

static int t1 = test_div<1, double>();
static int t2 = test_div<2, double>();
static int t3 = test_div<3, double>();
