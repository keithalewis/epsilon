// bell.h - Complete Bell polynomials
// B_{n+1}(a_1,...,x_{n+1}) = sum_{k = 0}^n C(n,k) B_{n-k}(a_1,...,a_{n-k}) a_{k+1}
// eetermined by exp(sum_{n>0} a_n x^n/n!) = sum_{n>=0} B_n(a_1,...,a_n) x^n/n!
#pragma once
#include <iterator>
#include <vector>

namespace fms {

    template<class S, class X = std::iterator_traits<S>::value_type>
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
            if (n < B.size())
                return B[n];
            for (size_t k = B.size(); k <= n; ++k) {
                X Bk = 0;
                X Ckj = 1; // really C(k-1,j)
                S ak = a;
                for (size_t j = 0; j < k /*&& ak*/; ++j) {
                    Bk += Ckj * operator[](k - 1 - j) * (*ak);
                    Ckj *= k - 1 - j;
                    Ckj /= j + 1;
                    ++ak;
                }
                B.push_back(Bk);
            }

            return B[n];
        }
    };

}