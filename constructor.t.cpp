#include "epsilon.h"
#include <cassert>

using namespace fms;

template <size_t N, class X> int test_constructor()
{
    {
        epsilon<N, X> e0(0, 1);
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
    {
        X e[N];
        for (int i = 0; i < N; ++i) {
            e[i] = X(i);
        }
        epsilon<N, X> e0(e);
        for (int i = 0; i < N; ++i) {
            assert(e0[i] == X(i));
        }
        auto e1{e0};
        assert(e1 == e0);
        e0 = e1;
        assert(e0 == e1);
    }
    {
        epsilon<3, X> e0{0, 1, 2};
        assert(e0[0] == 0);
        assert(e0[1] == 1);
        assert(e0[2] == 2);
    }

    return 0;
}

static int t1 = test_constructor<1, double>();
static int t2 = test_constructor<2, double>();
static int t3 = test_constructor<3, double>();
static int t4 = test_constructor<4, double>();
