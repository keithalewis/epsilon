// test expressions
#include <assert.h>
#include "epsilon.h"

using namespace fms;

template<class X, class Y, class Z>
auto test_expression()
{
    X x = 0;
    Y y = 1;
    auto z = Z(3) + epsilon<2, Z>{};

    return x + y * z - z/z;
}

#pragma warning(push)
#pragma warning(disable: 4244 4267)
static int t = [] {
    auto t1 = test_expression<double, double, double>();

    return 0;
}();
#pragma warning(pop)