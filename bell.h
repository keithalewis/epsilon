// bell.h - Complete Bell polynomials
// B_{n+1}(a_1,...,x_{n+1}) = sum_{k = 0}^n C(n,k) B_{n-k}(a_1,...,a_{n-k}) a_{k+1}.
// They satisfy exp(sum_{n>0} a_n x^n/n!) = sum_{n>=0} B_n(a_1,...,a_n) x^n/n!
// Taking a derivative with respect to x and equating equal powers results in this formula.
#pragma once
#include <iterator>
#include <vector>

namespace fms {

    template<class S, class X>
    class Bell {
        S a; // sequence of coefficients
        std::vector<X> B; // memoized Bell values
    public:
        Bell(S a)
            : a(a)
        {
            B.push_back(1);
        }
        X operator[](size_t n) {

            if (n >= B.size()) {
                for (size_t k = B.size(); k <= n; ++k) {
                    X Bk = 0;
                    X Ckj = 1; // really C(k-1,j)
                    S ak = a;
                    for (size_t j = 0; j < k && ak; ++j) {
                        Bk += Ckj * B[k - 1 - j] * (*ak);
                        Ckj *= k - 1 - j;
                        Ckj /= j + 1;
                        ++ak;
                    }
                    B.push_back(Bk);
                }
            }

            return B[n];
        }
    };

}