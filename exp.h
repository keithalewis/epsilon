// exp.h - exponential function
#pragma once
#include "epsilon.h"

namespace fms {

    template<class X>
    inline X exp(const X& x)
    {
        //!!! normalize x
        // Use frexp and ldrexp.
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
    
    // x = x0 I + x1 J + x2 J^2 + ...
    // exp(x) = exp(x0 I) exp(x1 J) exp(x2 J^2) ...
    //        = expN<0>(x0) * expN<1>(x1) * ...
    template<size_t N, class X = double>
    inline epsilon<N,X> expN(size_t i, const X& x)
    {
	if (0 == i) {
	    return exp(x)*epsilon<N,X>(1); // epsilon(1) is the identity matrix
	}
	if (1 == i) {
	    std::vector<X> xn(N);
	    X xn_ = 1; // xn_ = x^n/n!
	    xn[0] = xn_;
	    for (size_t i = 1; i < N; ++i) {
		xn_ *= x/i;
		xn[i] = xn_;
	    }
		
	    return epsilon<N,X>(xn.data());
	}
	// ... this does not work!!!
	// We don't know N so we don't know how many cases to implement
	    
	return epsilon<N,X>();
    }

    // Return exp(x J^I)
    template<size_t I, size_t N, class X = double>
    inline epsilon<N,X> expIN(const X& x)
    {
	static_assert (0 != N);
	if constexpr (0 == I) {
            return exp(x)*epsilon<N,X>(1);
	}
	
	 std::vector<X> xn(N);
	 X xn_ = 1; // xn_ = x^n/n!
	 xn[0] = 1;
	 for (i = 1; i*I < N; ++i) {
	     xn_ *= x/i;
	     xn[i*I] = xn_;
	 }
	 
	 return epsilon<N,X>(xn.data());
    }
	
    // Implement exp(x) using expIN
}
