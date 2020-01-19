#include <cassert>
#include "normal.h"
#include <iostream>
using namespace fms;

template<class X>
int test_erf()
{
	X x = 1;
	X ex = fms::erf<X>(x);

	X dx = ex - ::erf(x);
	assert(fabs(dx) <= 1 * std::numeric_limits<X>::epsilon());

	return 0;
}
static int test_erf_double = test_erf<double>();


template<size_t N, class X>
int test_erfn()
{
	X x = 1;
	epsilon<N, X> x_(x, 1);
	auto ex = fms::erf(x_);

	X erfx = ::erf(x);	
	assert(fabs(ex[0] - erfx) <= 1 * std::numeric_limits<X>::epsilon());
	if constexpr (N > 1) {
		assert(fabs(ex[1] - M_2_SQRTPI * ::exp(-1)) <= 1 * std::numeric_limits<X>::epsilon());
		assert(fabs(ex[2] + 2*M_2_SQRTPI * ::exp(-1)) <= 1 * std::numeric_limits<X>::epsilon());
	}

	x_ = epsilon<N, X>({ 1,1,1 });
	ex = fms::erf(x_);

	erfx = ::erf(x);
	assert(fabs(ex[0] - erfx) <= 1 * std::numeric_limits<X>::epsilon());
	if constexpr (N > 1) {
		assert(fabs(ex[1] - M_2_SQRTPI * ::exp(-1)) <= 1 * std::numeric_limits<X>::epsilon());
		assert(fabs(ex[2] + M_2_SQRTPI * ::exp(-1)) <= 1 * std::numeric_limits<X>::epsilon());
	}
	return 0;
}
static int test_erfn_double = test_erfn<3, double>();

template<size_t N, class X>
int test_erf2n()
{
	X x = 1;
	epsilon<N, X> x_(x, 1);
	auto ex = fms::erf2(x_);

	X erfx = ::erf(x);
	assert(fabs(ex[0] - erfx) <= 1 * std::numeric_limits<X>::epsilon());
	if constexpr (N > 1) {
		assert(fabs(ex[1] - M_2_SQRTPI * ::exp(-1)) <= 1 * std::numeric_limits<X>::epsilon());
		assert(fabs(ex[2] + 2 * M_2_SQRTPI * ::exp(-1)) <= 1 * std::numeric_limits<X>::epsilon());
	}

	x_ = epsilon<N, X>({ 1,1,1 });
	ex = fms::erf2(x_);

	erfx = ::erf(x);
	assert(fabs(ex[0] - erfx) <= 1 * std::numeric_limits<X>::epsilon());
	if constexpr (N > 1) {
		assert(fabs(ex[1] - M_2_SQRTPI * ::exp(-1)) <= 1 * std::numeric_limits<X>::epsilon());
		assert(fabs(ex[2] + M_2_SQRTPI * ::exp(-1)) <= 1 * std::numeric_limits<X>::epsilon());
	}
	return 0;
}
static int test_erf2n_double = test_erf2n<3, double>();