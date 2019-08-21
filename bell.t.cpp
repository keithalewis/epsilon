#include <cassert>
#include "bell.h"

using namespace fms;

template<class X>
int test_bell()
{
    X a[] = { 1, 1, 1, 1, 1 };
    {
        Bell<X*, X> B(a);
        assert(B[5] == 52);
        assert(B[4] == 15);
        assert(B[3] == 5);
        assert(B[2] == 2);
        assert(B[1] == 1);
        assert(B[0] == 1);
    }
    {
        Bell<X*, X> B(a);
        assert(B[0] == 1);
        assert(B[1] == 1);
        assert(B[2] == 2);
        assert(B[3] == 5);
        assert(B[4] == 15);
        assert(B[5] == 52);
    }

    return 0;
}

int bell = test_bell<double>();

