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
	
	template<size_t N, class Y = double, class = IsArithmetic<Y>>
	inline epsilon<N, Y> exp2(const epsilon<N, Y>& x) {
		using X = epsilon<N, Y>;
		X res;
		res += ::exp(x.a[0]);
		//decompose x by each of its element
		// a b c   a 0 0   0 b 0   0 0 c
		// 0 a b = 0 a 0 + 0 0 b + 0 0 0
		// 0 0 a   0 0 a   0 0 0   0 0 0
		for (size_t i = 1; i < N; i++) {
			auto num = x[i];
			//write num=A * 2^E
			//with 0.5<=abs(A)<1
			int E;
			auto A = std::frexp(num, &E);
			if (E > 0) {
				//epsilon_A=[0,...,A,...,0]
				X ex = 1;

				// (epsilon_A)^n/n!
				X xn_;
				xn_.a[i] = A;
				X xn__(xn_);
				int n = 1;
				while (xn_ + X(1) != X(1)) {
					ex += xn_;
					xn_ *= xn__ / ++n;
				}
				auto power = std::ldexp(1.0, E);
				while (power) {
					res *= ex;
					power--;
				}
			}
			else {
				//epsilon_num=[0,...,num,...,0]
				X ex = 1;
				// (epsilon_num)^n/n!
				X xn_;
				xn_.a[i] = num;

				X xn__(xn_);

				int n = 1;
				while (fabs(xn_) + X(1) != X(1)) {
					ex += xn_;
					xn_ *= xn__ / ++n;
				}
				res *= ex;
			}
			
		}
		return res;
	}
}
