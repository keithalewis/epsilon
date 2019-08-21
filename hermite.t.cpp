#include <cassert>
#include "hermite.h"

using namespace fms;

template<class X>
int test_hermite()
{
    X x;
    x = 1;
    assert(Hermite(2, x) == x * x - 1);
    x = 2;
    assert(Hermite(2, x) == x * x - 1);

    return 0;
}

static int hermite = test_hermite<double>();
