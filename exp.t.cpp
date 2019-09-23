#include <cassert>
#include "exp.h"

using namespace fms;

template<class X = double>
int test_exp()
{
    X x = 1;
    X ex = fms::exp<X>(x);

    X dx = ex - ::exp(x);
    assert(fabs(dx) <= 2*std::numeric_limits<X>::epsilon());

    return 0;
}

static int t1 = test_exp<double>();
