#include <cassert>
#include "exp.h"

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
