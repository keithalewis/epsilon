// epsilon.h - machine precision derivatives 
// f(x + epsilon) = sum_{n >= 0} f^{(n)}(x) epsilon^n/n!

#pragma once
#include <algorithm>
#include <cmath>
#include <compare>
#include <concepts>
#include <type_traits>
#include <valarray>

template<std::floating_point X>
inline auto operator<=>(const std::valarray<X>& a, const std::valarray<X>& b)
{
    if (a.size() != b.size()) {
        return X(a.size()) <=> X(b.size());
    }

    return a.min() <=> b.max();
}

namespace fms {

    template<size_t N, std::floating_point X = double>
    class epsilon {
        // a[0], ..., a[N-1] represents sum_{n < N} a_n epsilon^n/n!
        std::valarray<X> a;
    public:
        // x + epsilon
        epsilon(const X& x = X(0))
            : a(X(0), N)
        {
            if constexpr (N > 0) {
                a[0] = x;
            }
            if constexpr (N > 1) {
                a[1] = X(1);
            }
        }
        // from array
        epsilon(size_t n, const X* a_)
            : a(N)
        {
            for (size_t i = 0; i < n and i < N; ++i) {
                a[i] = a_[i];
            }
        }
        epsilon(const epsilon&) = default;
        epsilon& operator=(const epsilon&) = default;
        epsilon(epsilon&&) = default;
        epsilon& operator=(epsilon&&) = default;
        ~epsilon()
        { }

        // not all 0
        explicit operator bool() const
        {
            return (a == X(0)).min() == false;
        }
        /*
        auto operator<=>(const epsilon&) const = default;
        */
        bool operator==(const epsilon& b) const 
        {
            return (a == b.a).min() == true;
        }
        bool operator<(const epsilon& b) const
        {
            return std::lexicographical_compare(a.begin(), a.end(), b.a.begin(), b.a.end(), std::less<X>{});
        }
        bool operator==(const X& x) const
        {
            if (a[0] != x) {
                return false;
            }
            for (size_t i = 1; i < size(); ++i) {
                if (a[i] != 0) {
                    return false;
                }
            }

            return true;
        }
        size_t size() const noexcept
        {
            return N;
        }
        // derivatives
        X operator[](size_t n) const
        {
            return n < N ? a[n] : 0;
        }
        X& operator[](size_t n)
        {
            return a[n];
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
        // (sum a_j e^j/j!) (sum b_k e^k/k!)
        // = sum_n sum_{j + k = n} a_j b_{n-j} e^n/j!k!
        // = sum_n [sum_{j + k = n} C(n,j) a_j b_{n-j}] e^n/n!
        epsilon& operator*=(const epsilon& b)
        {
            std::valarray<X> c(X(0), N);

            for (size_t n = 0; n < N; ++n) {
                X Cnk = 1;
                for (size_t k = 0; k <= n; ++k) {
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
            a *= b;

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
            for (size_t n = 0; n < N; ++n) {
                X Cnk = 1;
                c[n] = a[n];
                for (size_t k = 0; k < n; ++k) {
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
            a /= b;

            return *this;
        }

        // L-p norm
		X norm(X p = 1) const
		{
            if (p == std::numeric_limits<X>::infinity()) {
                return abs(a).max();
            }
            if (p == 1) {
                return abs(a).sum();
            }

			return std::pow(pow(abs(a), p).sum(), X(1)/p);
		}

    };

} // namespace fms

//
// Global operators
//

// change size
template<size_t M, size_t N, std::floating_point X>
inline fms::epsilon<M> resize(const fms::epsilon<N,X>& a)
{
    return fms::epsilon<M>(N, &a[0]);
}

template<size_t N, std::floating_point X>
inline X fabs(const fms::epsilon<N,X>& a)
{
	return a.norm(X(1));
}
template<size_t N, std::floating_point X>
inline X norm(const fms::epsilon<N, X>& a, const X& p = 1)
{
	return a.norm(p);
}

template<size_t N, std::floating_point X>
inline fms::epsilon<N,X> operator+(fms::epsilon<N,X> a, const fms::epsilon<N,X>& b)
{
    return a += b;
}
template<size_t N, std::floating_point X, std::floating_point Y>
inline fms::epsilon<N, X> operator+(fms::epsilon<N, X> a, const Y& b)
{
    return a += b;
}
template<size_t N, std::floating_point X, std::floating_point Y>
inline fms::epsilon<N, X> operator+(const Y& a, fms::epsilon<N, X> b)
{
    return b += a;
}

template<size_t N, std::floating_point X>
inline fms::epsilon<N, X> operator-(fms::epsilon<N, X> a, const fms::epsilon<N, X>& b)
{
    return a -= b;
}
template<size_t N, std::floating_point X, std::floating_point Y>
inline fms::epsilon<N, X> operator-(fms::epsilon<N, X> a, const Y& b)
{
    return a -= b;
}
template<size_t N, std::floating_point X, std::floating_point Y>
inline fms::epsilon<N, X> operator-(const Y& a, fms::epsilon<N, X> b)
{
    return -(b -= a);
}

template<size_t N, std::floating_point X>
inline fms::epsilon<N, X> operator*(fms::epsilon<N, X> a, const fms::epsilon<N, X>& b)
{
    return a *= b;
}
template<size_t N, std::floating_point X, std::floating_point Y>
inline fms::epsilon<N, X> operator*(fms::epsilon<N, X> a, const Y& b)
{
    return a *= b;
}
template<size_t N, std::floating_point X, std::floating_point Y>
inline fms::epsilon<N, X> operator*(const Y& a, fms::epsilon<N, X> b)
{
    return b *= a;
}

template<size_t N, std::floating_point X>
inline fms::epsilon<N, X> operator/(fms::epsilon<N, X> a, const fms::epsilon<N, X>& b)
{
    return a /= b;
}
template<size_t N, std::floating_point X, std::floating_point Y>
inline fms::epsilon<N, X> operator/(fms::epsilon<N, X> a, const Y& b)
{
    return a /= b;
}
template<size_t N, std::floating_point X, std::floating_point Y>
inline fms::epsilon<N, X> operator/(const Y& a, fms::epsilon<N, X> b)
{
    return (a + fms::epsilon<N,X>()) /= b;
}

template<size_t N, std::floating_point X>
inline fms::epsilon<N, X> operator-(fms::epsilon<N, X> x)
{
    return X(-1) * x;
}


