#include <cassert>
#include "exp.h"
#include <iostream>
using namespace fms;

template<class X>
int test_exp()
{
    X x = 1;
    X ex = fms::exp<X>(x);

    X dx = ex - ::exp(x);
    assert(fabs(dx) <= 2*std::numeric_limits<X>::epsilon());

    return 0;
}
static int test_exp_double = test_exp<double>();

template<size_t N, class X>
int test_expn()
{
	X x = 1;
	epsilon<N, X> x_(x, 1);
	auto ex = fms::exp(x_);

	X expx = ::exp(x);	
	assert(fabs(ex[0] - expx) <= 2 * std::numeric_limits<X>::epsilon());
	if constexpr (N > 1) {
		assert(fabs(ex[1] - expx) <= 2 * std::numeric_limits<X>::epsilon());
	}

	return 0;
}
static int test_expn_double = test_expn<2,double>();

template<size_t N, class X>
int test_expn_exp2n_at_small()
{
	X x = 0.1;
	epsilon<N, X> x_(x, 1);
	auto ex = fms::exp(x_);

	X expx = ::exp(x);
	assert(fabs(ex[0] - expx) <= 2 * std::numeric_limits<X>::epsilon());
	if constexpr (N > 1) {
		assert(fabs(ex[1] - expx) <= 2 * std::numeric_limits<X>::epsilon());
	}

	ex = fms::exp2(x_);
	
	assert(fabs(ex[0] - expx) <= 0 * std::numeric_limits<X>::epsilon());
	if constexpr (N > 1) {
		assert(fabs(ex[1] - expx) <= 0 * std::numeric_limits<X>::epsilon());
	}

	return 0;
}
static int test_expn_exp2n_at_small_double = test_expn_exp2n_at_small<2, double>();

template<size_t N, class X>
int test_expn_exp2n_at_large()
{
	X x = 10;
	epsilon<N, X> x_(x, 1);
	auto ex = fms::exp(x_);

	X expx = ::exp(x);
	
	assert(fabs(ex[0] - expx) <= 16384 * std::numeric_limits<X>::epsilon());
	if constexpr (N > 1) {
		assert(fabs(ex[1] - expx) <= 16384 * std::numeric_limits<X>::epsilon());
	}

	ex = fms::exp2(x_);	
	assert(fabs(ex[0] - expx) <= 0 * std::numeric_limits<X>::epsilon());
	if constexpr (N > 1) {
		assert(fabs(ex[1] - expx) <= 0 * std::numeric_limits<X>::epsilon());
	}

	return 0;
}
static int test_expn_exp2n_at_large_double = test_expn_exp2n_at_large<2, double>();

template<size_t N, class X>
int test_expn_higher_order()
{
	epsilon<N, X> x_({ 10,2,0.1 });
	auto ex = fms::exp2(x_);

	std::vector<X> expx({ ::exp(10),::exp(10) * 2, ::exp(10) * 4.1 });

	assert(fabs(ex[0] - expx[0]) <= 0 * std::numeric_limits<X>::epsilon());
	if constexpr (N > 1) {
		assert(fabs(ex[1] - expx[1]) <= 0 * std::numeric_limits<X>::epsilon());		
		assert(fabs(ex[2] - expx[2]) <= 65536 * std::numeric_limits<X>::epsilon());
	}

	ex = fms::exp3(x_);
	assert(fabs(ex[0] - expx[0]) <= 0 * std::numeric_limits<X>::epsilon());
	if constexpr (N > 1) {
		assert(fabs(ex[1] - expx[1]) <= 0 * std::numeric_limits<X>::epsilon());
		assert(fabs(ex[2] - expx[2]) <= 0 * std::numeric_limits<X>::epsilon());
	}
	return 0;
}
static int test_expn_higher_order_double = test_expn_higher_order<3, double>();


template<class X>
X p(const X& x)
{
	return 1 + 2 * x + 3 * x * x;
}
template<class X>
X dp(const X& x)
{
	return 2 + 6 * x;
}
template<class X>
X ddp([[maybe_unused]] const X& x)
{
	return 6;
}

/*
template<class X>
X q(const X& x)
{
	return fms::exp(-x * x);
}
template<class X>
X dq(const X& x)
{
	return -2*x* fms::exp(-x * x);
}
template<class X>
X ddq([[maybe_unused]] const X& x)
{
	return (4*x*x-2) * fms::exp(-x * x);
}

int test_exp_derivative()
{
	
	double x = 1;
	epsilon<3, double> e(x, 1);
	auto qe = q(e);
	assert(fabs(qe[0] - q(x)) <= 1 * std::numeric_limits<double>::epsilon());
	assert(fabs(qe[1] - dq(x)) <= 1 * std::numeric_limits<double>::epsilon());
	assert(fabs(qe[2] - ddq(x)) <= 1 * std::numeric_limits<double>::epsilon());
	return 0;
}

static int exp_derivative = test_exp_derivative();*/