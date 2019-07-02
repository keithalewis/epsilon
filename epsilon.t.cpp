// epsilon.t.cpp - test epsilon
#include <cassert>
#include "epsilon.h"
using namespace fms;

template<size_t N, class X>
void test_constructor()
{
    epsilon<N,X> e0;
    for (int i = 0; i < N; ++i) {
        assert(e0[i] == X(i == 1 ? 1 : 0));
    }
    auto e1 = e0;
    assert (e1 == e0);
    assert (e0 == e1);
}

int main()
{
    test_constructor<1, double>();
    test_constructor<2, double>();
    test_constructor<3, double>();
    test_constructor<4, double>();

    return 0;
}