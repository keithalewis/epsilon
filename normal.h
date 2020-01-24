#pragma once
#include "epsilon.h"
#include "exp.h"
#include <vector>
#include <cmath>
#include "corecrt_math_defines.h"
namespace fms {
	//erf(x)=2/sqrt(pi)*(x-x^3/3+x^5/(5*2!)-x^7/(7*3!)+...
	template<class X>
	inline X erf(const X& x)
	{
		X ex = 0.0;
		X xn_(x); // x^(2n+1)/(n!*(2n+1))

		int n = 0;
		while (xn_ / (2 * n + 1) + ex != ex) {
			if (n % 2 == 0)
				ex += xn_ / (2 * n + 1);
			else
				ex -= xn_ / (2 * n + 1);
			xn_ *= x * x / ++n;
		}
		return M_2_SQRTPI * ex;
	}

	template<>
	inline double erf(const double& x) {
		return std::erf(x);
	}

	//erf(x)=2/sqrt(pi) * Int_0^x exp(-t*t) dt
	//		=2/sqrt(pi) * Int_0^(aI+B) exp(-t*t) dt
	//		=2/sqrt(pi) [Int_0^aI exp(-t*t) dt + exp(-a*a) B + exp(-a*a)'/2 B^2 + ...]
	//		=erf(a) I + 2/sqrt(pi) [exp(-a*a) B + exp(-a*a)'/2 B^2 + ...]
	template<size_t N, class X = double>
	epsilon<N, X> erf2(const epsilon<N, X>& x) {
		epsilon<N, X> res(erf(x[0]));
		epsilon<N, X> B = x - epsilon<N, X>(x[0]);

		epsilon<N, X> Bn_(B); // B^n/n!

		epsilon<N, double> e(x[0], 1);
		auto fe = f(e);

		int n = 1;
		while (n < N && fabs(Bn_) + X(1) != X(1)) {
			res += Bn_ * M_2_SQRTPI * fms::Dexp(n - 1, x[0]);
			//res += Bn_ * M_2_SQRTPI * fe[n-1];
			Bn_ *= B / ++n;
		}
		return res;
	}

	template<class X>
	X f(X i) { return exp(-i * i); };

	//d^n/(dx^n)[exp(-x^2)] = 2^n e^(-x^2) (-x)^n n! sum_(k=0)^floor(n/2) ((-4)^(-k) x^(-2 k))/(k! (-2 k + n)!) for (n element Z and n>=0)
	template<class X = double>
	X Dexp(size_t N, const X & x) {

		auto factorial = [](const size_t& i) {return std::tgamma(i + 1); };
		X res = 0;
		for (size_t k = 0; k <= std::floor(N / 2.0); k++) {
			res += pow(-4, -k) * pow(x, -2 * k) / factorial(k) / factorial(-2 * k + N);
		}
		return res * std::ldexp(1.0, N) * ::exp(-x * x) * pow(-x, N);
	}

	//PHI(x)=1/2+erf(x/sqrt2)/2
	template<class X>
	inline X cdf(const X& x)
	{
		return X(0.5) + erf(x * M_SQRT1_2) * 0.5;
	}

	template<class X>
	inline X cdf2(const X& x)
	{
		return X(0.5) + erf2(x * M_SQRT1_2) * 0.5;
	}

}