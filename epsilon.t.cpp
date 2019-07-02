// epsilon.t.cpp - test epsilon
#include <cassert>
#include "epsilon.h"
using namespace fms;

template<size_t N, class X>
void test_constructor()
{
    {
        epsilon<N,X> e0;
        for (int i = 0; i < N; ++i) {
            assert(e0[i] == X(i == 1 ? 1 : 0));
        }
        auto e1 = e0;
        assert (e1 == e0);
        assert (e0 == e1);
        assert (!(e0 != e1));
        e0 = e1;
        assert(e0 == e1);
    }
    {
        X e[N];
        for (int i = 0; i < N; ++i) {
            e[i] = X(i);
        }
        epsilon<N,X> e0(e);
        for (int i = 0; i < N; ++i) {
            assert(e0[i] == X(i));
        }
        auto e1{e0};
        assert (e1 == e0);
        e0 = e1;
        assert (e0 == e1);
    }
    {
        epsilon<3,X> e0{0,1,2};
        assert (e0[0] == 0);
        assert (e0[1] == 1);
        assert (e0[2] == 2);
    }
}

template<size_t N, class X = double>
void test_add()
{
    epsilon<N,X> e0(X(0));
    e0 += e0;
    for (size_t i = 0; i < N; ++i) {
        assert (e0[i] == X(i == 1 ? 2 : 0));
    }
}
template<size_t N, class X = double>
void test_sub()
{
    epsilon<N, X> e0(X(0));
    e0 -= e0;
    for (size_t i = 0; i < N; ++i) {
        assert(e0[i] == X(0));
    }
}

template<size_t N, class X = double>
void test_mul()
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
}

template<class X>
X poly(const X& x)
{
    return 1. + 2.*x + x*x;
}

void test_derivative()
{
    epsilon<3,double> e(1);
    auto pe = 1. + 2.*e + 3.*e*e;
    /*
    assert (pe[0] == 1 + 2*1 + 1*1);
    assert (pe[1] == 2 + 2*1);
    assert (pe[2] == 2);
    */
}

int main()
{
    test_constructor<1, double>();
    test_constructor<2, double>();
    test_constructor<3, double>();
    test_constructor<4, double>();
    test_add<1, double>();
    test_add<2, double>();
    test_add<3, double>();
    test_add<4, double>();
    test_sub<1, double>();
    test_sub<2, double>();
    test_sub<3, double>();
    test_sub<4, double>();
    test_mul<1, double>();
    test_mul<2, double>();
    test_mul<3, double>();
    test_mul<4, double>();

    test_derivative();

    return 0;
}