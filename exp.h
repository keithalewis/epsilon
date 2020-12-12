// exp.h - exponential function
#pragma once
#include <vector>
#include "bell.h"

namespace fms {

	// exp(sum_{n >= 0} a_n x^n/n!) = sum_{n>0} B_n(a_1,...,a_n) x^n/n!
	template<class X>
	inline X exp(const X& x) 
	{
		size_t n = x.size();
		std::vector<X> B(n);
		Bell(n, x, &B[0]);
		X ex = X(1); // exp(x)
		X xn = X(1); // n^n
		size_t n_ = 1; // n!
		for (size_t k = 0; k < n; ++k) {
			ex += B[k] * xn / X(n_);
			xn *= x;
			if (k) {
				n_ *= k;
			}
		}

		return ex;
	}
}
