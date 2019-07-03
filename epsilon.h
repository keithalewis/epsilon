// epsilon.h - Machine precision derivative 
#pragma once
#include <valarray>

namespace fms {

    template<size_t N, class X = double>
    class epsilon {
        // a[0], ..., a[N-1] <-> sum_{k < N} a_k epsilon^k/k!
        std::valarray<X> a;
    public:
        // xI + epsilon
        epsilon(X x = 0)
            : a(X(0), N)
        {
            a[0] = x;
            if (N > 1)
                a[1] = 1;
        }
        epsilon(const X* px)
            : a(px, N)
        { }
        epsilon(std::initializer_list<X> il)
            : a{il}
        { }

        bool operator==(const epsilon& b) const
        {
            return (a == b.a).min() == true;
        }
        bool operator!=(const epsilon& b) const
        {
            return !operator==(b);
        }
        // bool operator< ...

        // underlying raw valarray
        const std::valarray<X>& array() const
        {
            return a;
        }
        // n-th derivative
        const X& operator[](size_t n) const
        {
            return a[n];
        }

        epsilon& operator-()
        {
            operator*=(X(-1));

            return *this;
        }

        epsilon& operator+=(const epsilon& b)
        {
            a += b.a;

            return *this;
        }
        epsilon& operator+=(const X& b)
        {
            a[0] += b;

            return *this;
        }
        epsilon& operator-=(const epsilon& b)
        {
            a -= b.a;

            return *this;
        }
        epsilon& operator-=(const X& b)
        {
            a[0] -= b;

            return *this;
        }
        // sum a_j e^j/j! sum b_k e^k/k!
        // = sum_n sum_{j + k = n} a_j b_{n-j} e^n/j!k!
        // = sum_n [sum_{j + k = n} C(n,j) a_j b_{n-j}] e^n/n!
        epsilon& operator*=(const epsilon& b)
        {
            std::valarray<X> c(X(0), N);

            for (int n = 0; n < N; ++n) {
                X Cnk = 1;
                for (int k = 0; k <= n; ++k) {
                    c[n] += Cnk*a[k]*b.a[n-k];
                    Cnk *= n - k;
                    Cnk /= k + 1;
                }
            }
            std::swap(c,a);

            return *this;
        }
        epsilon& operator*=(const X& b)
        {
            a[0] *= b;

            return *this;
        }

        // a/b = c iff a = b*c = sum_n [sum C(n,j) b_{n-j} c_j] e^n/n!
        // a_0 = b_0 c_0,
        //   c_0 = a_0/b_0.
        // a_1 = b_1 c_0 + b_0 c_1,
        //   c_1 = (a_1 - b_1 c_0)/b_0.
        // ...
        // a_n = C(n,0) b_n c_0 + C(n,1) b_{n-1} c_1 + .... + C(n,n) b_0 c_n,
        //   c_n = (a_n - C(n,0) b_n c_0 - ... - C(n, n-1) b_1 c_{n-1})/b_0.
        epsilon& operator/=(const epsilon& b)
        {
            std::valarray<X> c(X(0), N);
            X b0 = b[0];
            for (int n = 0; n < N; ++n) {
                X Cnk = 1;
                c[n] = a[n];
                for (int k = 0; k < n; ++k) {
                    c[n] -= Cnk*b.a[k]*c[n-k];
                    Cnk *= n - k;
                    Cnk /= k + 1;
                }
                c[n] /= b0;
            }
            std::swap(c, a);

            return *this;
        }
        epsilon& operator/=(const X& b)
        {
            a[0] /= b;

            return *this;
        }

    };

} // namespace fms

template<size_t N, class X = double>
inline fms::epsilon<N,X> operator+(fms::epsilon<N,X> a, const fms::epsilon<N,X>& b)
{
    return a += b;
}
template<size_t N, class X = double>
inline fms::epsilon<N, X> operator+(fms::epsilon<N, X> a, const X& b)
{
    return a += b;
}
template<size_t N, class X = double>
inline fms::epsilon<N, X> operator+(const X& a, fms::epsilon<N, X> b)
{
    return b += a;
}

template<size_t N, class X = double>
inline fms::epsilon<N, X> operator-(fms::epsilon<N, X> a, const fms::epsilon<N, X>& b)
{
    return a -= b;
}
template<size_t N, class X = double>
inline fms::epsilon<N, X> operator-(fms::epsilon<N, X> a, const X& b)
{
    return a -= b;
}
template<size_t N, class X = double>
inline fms::epsilon<N, X> operator-(const X& a, fms::epsilon<N, X> b)
{
    return -(b -= a);
}

template<size_t N, class X = double>
inline fms::epsilon<N, X> operator*(fms::epsilon<N, X> a, const fms::epsilon<N, X>& b)
{
    return a *= b;
}
template<size_t N, class X = double>
inline fms::epsilon<N, X> operator*(fms::epsilon<N, X> a, const X& b)
{
    return a *= b;
}
template<size_t N, class X = double>
inline fms::epsilon<N, X> operator*(const X& a, fms::epsilon<N, X> b)
{
    return b *= a;
}
