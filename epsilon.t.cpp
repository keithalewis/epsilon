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
X p(const X& x)
{
    return 1. + 2.*x + 3.*x*x;
}
template<class X>
X dp(const X& x)
{
    return 2. + 6.*x;
}
template<class X>
X ddp(const X& x)
{
    return 6.;
}

void test_derivative()
{
    {
        epsilon<5,double> e{1,1,1,1,1};
        e *= e;
        assert(e[0] == 1);
        assert(e[1] == 2);
        assert(e[2] == 4);
        assert(e[3] == 8);
        assert(e[4] == 16);

    }
    {
        epsilon<5,double> e(1.);
        auto e1 = e;
        e1 *= e;
        e1 *= e;
        e1 *= e;
    }
    {
        double x = 1;
        epsilon<3,double> e(x);
        auto pe = p(e);
        assert(pe[0] == p(x));
        assert(pe[1] == dp(x));
        assert(pe[2] == ddp(x));
    }
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