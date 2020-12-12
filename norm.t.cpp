// norm.t.cpp - Test epsilon norm
#include <cassert>
#include "epsilon.h"

using namespace fms;

template<size_t N, class X>
int test_norm()
{
	{
		epsilon<N,X> e;
		if constexpr (N > 1) {
			assert(1 == e.norm());
		}
	}
	{
		X x = 1;
		auto e = x + epsilon<N, X>{};
		assert(x + X(1)*(N > 1) == e.norm());
		assert(sqrt(x*x + X(1) * (N > 1)) == e.norm(2));
	}

	return 0;
}
int test_norm_double1 = test_norm<1, double>();
int test_norm_double2 = test_norm<2, double>();
int test_norm_double3 = test_norm<3, double>();
