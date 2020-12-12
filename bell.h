// bell.h - Complete Bell polynomials
// B_{n+1}(a_1,...,a_{n+1}) = sum_{k = 0}^n C(n,k) B_{n-k}(a_1,...,a_{n-k}) a_{k+1}.
// They satisfy exp(sum_{n>0} a_n x^n/n!) = sum_{n>=0} B_n(a_1,...,a_n) x^n/n!
// Taking a derivative with respect to x and equating equal powers results in this formula.
#pragma once
#include <iterator>
#include <vector>

namespace fms {

    // In place computation of Bell polynomials
    // B_{n+1}(x_1,...,x_n) = sum_{k=0}^n C(n,k) B_{n-k}(x_1,...,x_{n-k}) x_{k+1}, B_0 = 1.
    // B_n(x_1,...,x_{n-1}) = sum_{k=0}^{n-1} C(n-1,k) B_{n-1-k}(x_1,...,x_{n-1-k}) x_{k+1}, B_0 = 1.
    // Note (x_0,...,x_{n-1}) results in (B_0,...,B_n)
    template<class X>
    inline void Bell(size_t N, const X* x, X* B)
    {
        if (N > 0) {
            B[0] = X(1);
        }

        for (size_t n = 1; n <= N; ++n) {
            // *++B = sum(take(n, choose(n)*reverse(take(n, B))*x));
            B[n] = 0;
            X Cnk = 1; // really C(n-1,k)
            for (size_t k = 0; k < n; ++k) {
                B[n] += Cnk * B[n - 1 - k] * x[k];
                Cnk *= n - 1 - k;
                Cnk /= k + 1;
            }
        }
    }

}