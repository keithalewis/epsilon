// exp.h - exponential function
#pragma once
#include <cstring>
#include <cmath>
#include "epsilon.h"
#include "bell.h"

namespace fms {

	// exp(sum_{n >= 0} a_n e^n/n!) = exp(a_0) sum_{n > 0} B_n(a_1,...,a_n) e^n/n!
	template<size_t N, class X>
	inline epsilon<N,X> exp(epsilon<N,X> x) 
	{
		static epsilon<N, X> B;

		double ex0 = ::exp(x[0]);
		Bell(N, &x[0], &B[0]);
		std::swap(x, B);
		x *= ex0;

		return x;
	}
}
