#include "multi_epsilon.h"
#include <cassert>
using namespace fms;
static int test_multinomial()
{
    {
        auto a_ = multi_epsilon::add_epsilon({0, 0}, 2);
        auto a = 2 * a_[0] + 3 * a_[1];
        auto a2 = a * a;
        assert(a2[multi_epsilon::rep({2, 0}, 2)] == 4);
        assert(a2[multi_epsilon::rep({0, 2}, 2)] == 9);
        assert(a2[multi_epsilon::rep({1, 1}, 2)] == 12);
        auto a0 = a2 - a2;
        for (int i = 0; i < pow(2 + 1, 2); i++)
            assert(a0[i] == 0);

        /*multinomial a({
                { i{ 1, 0 }, 2 }, // 2 x_0
                { i{ 0, 1 }, 3 }, // + 3 x_1
                });
        auto a2 = a * a; // 4 x_0^2 + 9 x_1^2 + 12 x_0 x_1
        assert(a2.size() == 3);
        auto a20 = a2[i{ 2, 0 }];
        assert(a20 == 4);
        auto a02 = a2[i{ 0, 2 }];
        assert(a02 == 9);
        auto a11 = a2[i{ 1, 1 }];
        assert(a11 == 12);

        auto a0 = a2 - a2;
        for (const auto& [k, v] : a0) {
                assert(v == 0);
        }*/
    }
    {
        auto a_ = multi_epsilon::add_epsilon({0, 0}, 1);
        auto a = 2 + 3 * a_[0];
        auto b_ = multi_epsilon::add_epsilon({0, 0}, 1);
        auto b = 4 * b_[0] + 5 * b_[1];
        auto ab = a + b;
        assert(ab[multi_epsilon::rep({0, 0}, 1)] == 2);
        assert(ab[multi_epsilon::rep({1, 0}, 1)] == 7);
        assert(ab[multi_epsilon::rep({0, 1}, 1)] == 5);

        /*multinomial a({
                { { i{ 0, 0 } }, 2 }, // 2
                { { i{ 1, 0 } }, 3 }  // + 3 x_0
                });
        multinomial b({
                { i{ 1, 0 }, 4 }, // 4 x_0
                { i{ 0, 1 }, 5 }, // + 5 x_1
                });
        auto ab = a + b;
        assert(ab.size() == 3);
        auto a00 = ab[i{ 0, 0 }];
        assert(2 == a00);
        auto a10 = ab[i{ 1, 0 }];
        assert(7 == a10);
        auto a01 = ab[i{ 0, 1 }];
        assert(5 == a01);*/
    }

    return 0;
}
static int test_multinomial_ = test_multinomial();

template <class X> X f1(const X &x)
{
    return x * x * x;
}
// f1_x = 3 x^2
// f1_xx = 6 x
// f1_xxx = 6
// f1_xxxx = 0
static int test_multinomial_derivatve1()
{
    double x = 2;
    auto X = multi_epsilon::add_epsilon({x}, 3)[0];
    auto F1 = f1(X);
    auto f1_0 = F1[multi_epsilon::rep({0}, 3)];
    assert(f1_0 == x * x * x);
    auto f1_1 = F1[multi_epsilon::rep({1}, 3)];
    assert(f1_1 == 3 * x * x);
    auto f1_2 = F1[multi_epsilon::rep({2}, 3)];
    assert(f1_2 == 6 * x / 2);
    auto f1_3 = F1[multi_epsilon::rep({3}, 3)];
    assert(f1_3 == 6 / (2 * 3));

    return 0;
    /*using i = multi_index;
    double x = 2;
    auto X = constant<1>(x) + epsilon<1>(0);
    auto F1 = f1(X);
    auto f1_0 = F1[i{ 0 }];
    assert(f1_0 == x*x*x);
    auto f1_1 = F1[i{ 1 }];
    assert(f1_1 == 3 * x * x);
    auto f1_2 = F1[i{ 2 }];
    assert(f1_2 == 6 * x/2);
    auto f1_3 = F1[i{ 3 }];
    assert(f1_3 == 6/(2*3));

    assert(F1.find(i{ 4 }) == F1.end());

    return 0;*/
}
static int test_multinomial_derivatve1_ = test_multinomial_derivatve1();

static int test_epsilon()
{
    auto X = multi_epsilon::add_epsilon({0, 0}, 2);
    assert(X[0] != X[1]);
    auto c = 1 + X[0] - X[0];
    auto e0 = X[0];
    auto e1 = X[1];
    auto e01 = e0 * e1;
    auto e10 = e1 * e0;
    auto ce = (c + e0) * (c + e1);

    assert(ce[fms::multi_epsilon::rep({0, 0}, 2)] == 1);
    assert(ce[fms::multi_epsilon::rep({0, 1}, 2)] == 1);
    assert(ce[fms::multi_epsilon::rep({1, 0}, 2)] == 1);
    assert(ce[fms::multi_epsilon::rep({1, 1}, 2)] == 1);

    auto ce2 = ce * ce;
    assert(ce2[fms::multi_epsilon::rep({0, 0}, 2)] == 1);
    assert(ce2[fms::multi_epsilon::rep({1, 0}, 2)] == 2);
    assert(ce2[fms::multi_epsilon::rep({0, 1}, 2)] == 2);
    assert(ce2[fms::multi_epsilon::rep({2, 0}, 2)] == 1);
    assert(ce2[fms::multi_epsilon::rep({1, 1}, 2)] == 4);
    assert(ce2[fms::multi_epsilon::rep({0, 2}, 2)] == 1);
    assert(ce2[fms::multi_epsilon::rep({2, 1}, 2)] == 2);
    assert(ce2[fms::multi_epsilon::rep({1, 2}, 2)] == 2);
    assert(ce2[fms::multi_epsilon::rep({2, 2}, 2)] == 1);

    return 0;
    /*using i = fms::multi_index;

    auto c = constant<2>(1);
    auto e0 = epsilon<2>(0);
    auto e1 = epsilon<2>(1);
    assert(e0 != e1);
    auto e01 = e0 * e1;
    auto e10 = e1 * e0;
    assert(e01 == e10);
    auto ce = (c + e0) * (c + e1); // c*c + c*e0 + c*e1 + e0*e1

    assert(ce.size() == 4);

    assert(ce.contains({ 0, 0 }));
    assert(ce[i({ 0, 0 })] == 1);

    assert(ce.contains({ 0, 1 }));
    assert(ce[i({ 0, 1 })] == 1);

    assert(ce.contains({ 1, 0 }));
    assert(ce[i({ 1, 0 })] == 1);

    assert(ce.contains({ 1, 1 }));
    assert(ce[i({ 1, 1 })] == 1);

    auto ce2 = ce * ce; // (1 + e0 + e1 + e0 e1) +
                                            // (e0 + e0^2 + e0 e1 + e0^2 e1) +
                                            // (e1 + e0 e1 + e1^2 + e0 e1^2) +
                                            // (e0 e1 + e0^2 e1 + e0 e1^2 + e0^2
    e1^2)
                                            // = 1 + 2e0 + 2e1 + e0^2 + 4 e0 e1
    + e1^2 + 2 e0^2 e1 + 2 e0 e1^2 + e0^2 e1^2 assert(ce2[i({ 0, 0 })] == 1);
    assert(ce2[i({ 1, 0 })] == 2);
    assert(ce2[i({ 0, 1 })] == 2);
    assert(ce2[i({ 2, 0 })] == 1);
    assert(ce2[i({ 1, 1 })] == 4);
    assert(ce2[i({ 0, 2 })] == 1);
    assert(ce2[i({ 2, 1 })] == 2);
    assert(ce2[i({ 1, 2 })] == 2);
    assert(ce2[i({ 2, 2 })] == 1);


    return 0;*/
}
static int test_epsilon_ = test_epsilon();

template <class X> X f2(const X &x, const X &y)
{
    return x * x * x + x * y;
}
// f2_x = 3 x^2 + y
// f2_xx = 6 x
// f2_xxx = 6
// f2_xxxx = 0
// f2_y = x
// f2_yy = 0
// f2_xy = 1
static int test_multinomial_derivatve2()
{
    double x = 2;
    double y = 3;
    auto X = multi_epsilon::add_epsilon({x, y}, 3)[0];
    auto Y = multi_epsilon::add_epsilon({x, y}, 3)[1];
    auto F2 = f2(X, Y);
    auto f2_0 = F2[multi_epsilon::rep({0, 0}, 3)];
    assert(f2_0 == x * x * x + x * y);
    auto f2_x = F2[multi_epsilon::rep({1, 0}, 3)];
    assert(f2_x == 3 * x * x + y);
    auto f2_xx = F2[multi_epsilon::rep({2, 0}, 3)];
    assert(f2_xx == 6 * x / 2);
    auto f2_xxx = F2[multi_epsilon::rep({3, 0}, 3)];
    assert(f2_xxx == 6 / (2 * 3));
    auto f2_y = F2[multi_epsilon::rep({0, 1}, 3)];
    assert(f2_y == x);
    auto f2_xy = F2[multi_epsilon::rep({1, 1}, 3)];
    assert(f2_xy == 1);
    return 0;
}
static int test_multinomial_derivatve2_ = test_multinomial_derivatve2();