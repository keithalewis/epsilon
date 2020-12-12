#include <cassert>
#include "epsilon.h"

using namespace fms;

template<class X>
int test_constructor()
{
    {
        epsilon<1, X> e;
        auto e2(e);
        e = e2;
        assert(e == e2);
        assert(e.size() == 1);
        assert(e[0] == 0);
    }
    {
        epsilon<2, X> e;
        auto e2(e);
        e = e2;
        assert(e == e2);
        assert(e.size() == 2);
        assert(e[0] == 0);
        assert(e[1] == 1);
    }
    {
        X x = 2;
        epsilon<3, X> e(x);
        auto e2(e);
        e = e2;
        assert(e == e2);
        assert(e.size() == 3);
        assert(e[0] == x);
        assert(e[1] == 1);
        assert(e[2] == 0);
    }
    /*
    {
        epsilon<N, X> e0;
        for (int i = 0; i < N; ++i) {
            assert(e0[i] == X(i == 1 ? 1 : 0));
        }
        auto e1 = e0;
        assert(e1 == e0);
        assert(e0 == e1);
        assert(!(e0 != e1));
        e0 = e1;
        assert(e0 == e1);
    }
    */
 
    return 0;
}

static int t1 = test_constructor<double>();
