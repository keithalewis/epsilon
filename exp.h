// exp.h - exponential function
#pragma once
#include "epsilon.h"
#include <vector>
namespace fms {

    template<class X>
    inline X exp(const X& x)
    {
        //!!! normalize x
        // Use frexp and ldrexp.
		//frexp(x): x= arg * 2^E
		//arg is in the range (-1;-0.5], [0.5; 1)
		//e^(arg*2^E) = (e^arg) ^ (2^E)
        // Requires log_e 2.
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
		static_assert (0 != N);
		//!!if constexpr (0 == I) {
		if (0 == I){
			return epsilon<N, X>(exp(x));
		}

		std::vector<X> xn(N);
		X xn_ = 1; // xn_ = x^n/n!
		xn[0] = 1;
		for (size_t i = 1; i * I < N; ++i) {
			xn_ *= x / i;
			xn[i * I] = xn_;
		}

		return epsilon<N, X>(xn.data());
	}

	template<size_t N, class X = double>
	epsilon<N, X> exp2(const epsilon<N, X>& x) {
		//!!epsilon<N, X> res(expIN<0, N, X>(x[0]));
		epsilon<N, X> res(expIN<N, X>(x[0], 0));
		for (size_t i = 1; i < N; i++)
			//!!res *= expIN<i, N, X>(x[i]);
			res *= expIN<N, X>(x[i], i);
		return res;
	}

}
