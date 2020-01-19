// exp.h - exponential function
#pragma once
#include "epsilon.h"
#include <vector>
#include <cmath>
namespace fms {

    template<class X>
    inline X exp(const X& x)
    {
        X ex = 1;
        X xn_(x); // x^n/n!
		
        int n = 1;
        while (fabs(xn_) + X(1) != X(1)) {
            ex += xn_;
            xn_ *= x / ++n;			
        }

        return ex;
    }
	template<>
	inline double exp(const double& x){
		return ::exp(x);
	}

	// Return exp(x J^I)
	//!!this definition isn't compatible with later exp2()
	//!!template<size_t I, size_t N, class X = double>
	//!!inline epsilon<N, X> expIN(const X& x)
	template<size_t N, class X=double>
	inline epsilon<N,X> expIN(const X& x, const size_t& I)
	{
		auto factorial = [](const size_t& i) {return std::tgamma(i + 1); };
		static_assert (0 != N);
		if (0 == I){
			return epsilon<N, X>(exp(x));
		}
		//to standard Toeplitz matrix representation
		auto x_ = x / factorial(I);
		std::vector<X> xn(N);
		X xn_ = 1; // xn_ = x^n/n!
		xn[0] = 1;
		for (size_t i = 1; i * I < N; ++i) {
			xn_ *= x_ / i;
			xn[i * I] = xn_;
			//back to our representation
			xn[i * I] *= factorial(i * I);
		}

		return epsilon<N, X>(xn.data());
	}

	template<size_t N, class X = double>
	epsilon<N, X> exp2(const epsilon<N, X>& x) {
		epsilon<N, X> res(expIN<N, X>(x[0], 0));
		for (size_t i = 1; i < N; i++)
			res *= expIN<N, X>(x[i], i);
		return res;
	}

	template<size_t N, class X = double>
	epsilon<N, X> exp3(const epsilon<N, X>& x) {
		epsilon<N, X> res(expIN<N, X>(x[0], 0));
		res *= exp(x - epsilon<N, X>(x[0]));
		return res;
	}

}
