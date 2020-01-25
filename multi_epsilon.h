#pragma once
#include <cstring>
#include <cassert>
#include <iostream>
#include <valarray>
#include <vector>
#include <numeric>
#include <algorithm>
using std::cout;
using std::endl;
using std::vector;
namespace fms {
	class multi_epsilon {
	public:
		size_t m;//how many variables
		size_t Size;//size of m_lpBuf
		size_t N;//N=2 if only consider 1st order derivative
		
		const double zero = 0;

		//m variable,  derivatives up to order n
		multi_epsilon(size_t m = 0, size_t n=0)
			: Size(pow(n + 1, m)), m_lpBuf(0.0, pow(n+1,m)), N(n+1), m(m) {};

		multi_epsilon(const multi_epsilon& rhs)
			:Size(rhs.Size), m_lpBuf(rhs.m_lpBuf), N(rhs.N), m(rhs.m)
		{};

		
		multi_epsilon(const std::initializer_list<double>& rhs, size_t m, size_t n)
			:m_lpBuf(rhs), Size(rhs.size()), N(n+1), m(m) { };

		multi_epsilon(multi_epsilon&& rhs) noexcept {
			this->Size = rhs.Size;
			this->N = rhs.N;
			this->m = rhs.m;
			this->m_lpBuf = std::move(rhs.m_lpBuf);
		};

		~multi_epsilon() {
			Size = 0;
		};


		multi_epsilon& operator=(const multi_epsilon& rhs) {
			this->Size = rhs.Size;
			this->m_lpBuf = rhs.m_lpBuf;
			this->N = rhs.N;
			this->m = rhs.m;
			return *this;
		};


		multi_epsilon& operator=(multi_epsilon&& rhs) noexcept {
			this->Size = rhs.Size;
			this->m_lpBuf = std::move(rhs.m_lpBuf);
			this->N = rhs.N;
			this->m = rhs.m;
			return *this;
		};



		bool operator == (const multi_epsilon& rhs) const {
			if (this->Size != rhs.Size) return false;
			if (this->m != rhs.m) return false;
			if (this->N != rhs.N) return false;
			return (m_lpBuf == rhs.m_lpBuf).min() == true;
		}

		bool operator !=(const multi_epsilon& rhs) const {
			return !operator==(rhs);
		}

		multi_epsilon& operator += (const multi_epsilon& rhs)
		{
			assert(rhs.Size == this->Size);
			assert(rhs.m == this->m);
			assert(rhs.N == this->N);
			m_lpBuf += rhs.m_lpBuf;
			return *this;
		}

		//this+rhs*I
		multi_epsilon& operator += (const double& rhs)
		{
			m_lpBuf[0] += rhs;
			return *this;
		}

		multi_epsilon& operator -= (const multi_epsilon& rhs)
		{
			assert(rhs.Size == this->Size);
			assert(rhs.m == this->m);
			assert(rhs.N == this->N);
			m_lpBuf -= rhs.m_lpBuf;
			return *this;
		}

		//this-rhs*I
		multi_epsilon& operator -= (const double& rhs)
		{
			m_lpBuf[0] -= rhs;
			return *this;
		}

		//matrix multiplication
		//only need to calculate the first row
		multi_epsilon& operator *= (const multi_epsilon& rhs)
		{
			//if (!rhs.m_lpBuf) return *this;
			assert(rhs.Size == this->Size);
			assert(rhs.m == this->m);
			assert(rhs.N == this->N);
			std::valarray<double> temp((double)0, this->Size);
			size_t i = 0;
			{
				temp = std::valarray<double>((double)0, Size);
				for (size_t j = 0; j < Size; j++) {
					double t = 0;
					for (size_t k = 0; k <= j; k++)
						t += operator()(i, k) * rhs(k, j);
					temp[j] = t;
				}
				//for (size_t j = i; j < this->Size; j++)
				//	m_lpBuf[j] = temp[j];
				std::swap(m_lpBuf, temp);
			}
			return *this;
		}

		//this*rhs
		multi_epsilon& operator *= (const double& rhs)
		{
			m_lpBuf *= rhs;
			return *this;
		}

		//matrix division
		multi_epsilon& operator /= (const multi_epsilon& rhs)
		{
			//if (!rhs.m_lpBuf) return *this;
			assert(rhs.Size == this->Size);
			assert(rhs.m == this->m);
			assert(rhs.N == this->N);
			multi_epsilon rhs_inv = rhs.inverse();
			operator*=(rhs_inv);
			return *this;
		}

		//this/rhs
		multi_epsilon& operator /= (const double& rhs)
		{
			m_lpBuf /= rhs;
			return *this;
		}

		//-1*this
		multi_epsilon operator -()
		{
			multi_epsilon res(*this);
			res*=(-1.0);
			return res;
		}

		const double& operator ()(size_t i, size_t j) const
		{
			assert(i < Size);
			assert(j < Size);
			assert(j >= i);
			size_t k = N;
			while (k <= Size) {
				if ((i % k) > (j % k))
					return zero;
				k *= N;
			}
			return m_lpBuf[j - i];
		}

		double& operator ()(size_t i, size_t j)
		{
			assert(i < Size);
			assert(j < Size);
			assert(j >= i);
			//if (i > j) return 0;//visiting lower triangle element
			size_t k = N;
			while (k <= Size) {
				assert((i % k) <= (j % k));
				k *= N;
			}
			return m_lpBuf[j - i];
		}

		const double& operator [](size_t i) const {
			return m_lpBuf[i];
		}

		double& operator [](size_t i) {
			return m_lpBuf[i];
		}

		//inverse(this)
		multi_epsilon inverse() const {
			multi_epsilon res(m, N-1);
			res[0] = 1.0 / m_lpBuf[0];
			for (size_t i = 1; i < Size; i++) {
				double cum_prod = 0;
				for (size_t j = 0; j < i; j++)
					cum_prod += res[j] * this->operator()(j, i);
				res[i] = -cum_prod / m_lpBuf[0];
			}
			return res;
		}

		

		//static method identity
		//n th order derivative
		//dimension of matrix=n+1
		static multi_epsilon identity(size_t m, size_t n) {
			multi_epsilon result(m,n);
			result += 1;
			return result;
		}


		//epsilon \otimes epsilon \otimes I_2
		//return 6
		static size_t rep(const std::vector<size_t>& order, const size_t n) {
			size_t res=order[0]+1;
			for (size_t i = 1; i < order.size(); i++) {
				res = (n+1)*(res-1)+order[i]+1;
			}
			return res - 1;
		}

		static vector<multi_epsilon> add_epsilon(const vector<double>& x, const size_t n) {
			vector<multi_epsilon> result;
			vector<size_t> order(x.size(), 0);
			for (size_t i = 0; i < x.size(); i++) {
				order[i] = 1;
				multi_epsilon temp = identity(x.size(), n);
				temp *= x[i];
				temp[rep(order, n)] = 1.0;
				result.emplace_back(std::move(temp));
				order[i] = 0;
			}
			return result;
		}
		
		//print current matrix
		void print() const {
			for (size_t i = 0; i < Size; i++) {
				for (size_t j = 0; j < Size; j++)
					if (j < i) std::cout << 0 << ' ';
					else std::cout << operator()(i, j) << ' ';
				std::cout << std::endl;
			}
		}
	private:
		std::valarray<double> m_lpBuf; // data container
	};
}
inline fms::multi_epsilon operator + (fms::multi_epsilon A, const fms::multi_epsilon& B)
{
	return A += B;
}
inline fms::multi_epsilon operator + (fms::multi_epsilon A, const double& B)
{
	return A += B;
}
inline fms::multi_epsilon operator + (const double& B, fms::multi_epsilon A)
{
	return A += B;
}
inline fms::multi_epsilon operator - (fms::multi_epsilon A, const fms::multi_epsilon& B)
{
	return A -= B;
}
inline fms::multi_epsilon operator - (fms::multi_epsilon A, const double& B)
{
	return A -= B;
}
inline fms::multi_epsilon operator - (const double& B, fms::multi_epsilon A)
{
	return A -= B;
}
inline fms::multi_epsilon operator * (fms::multi_epsilon A, const fms::multi_epsilon& B)
{
	A *= B;
	return A;
}
inline fms::multi_epsilon operator * (fms::multi_epsilon A, const double& B)
{
	A *= B;
	return A;
}
inline fms::multi_epsilon operator * (const double& B, fms::multi_epsilon A)
{
	return A *= B;
}
inline fms::multi_epsilon operator / (fms::multi_epsilon A, const fms::multi_epsilon& B)
{
	return A /= B;
}
inline fms::multi_epsilon operator / (fms::multi_epsilon A, const double& B)
{
	return A /= B;
}
inline fms::multi_epsilon operator / (const double& B, fms::multi_epsilon A)
{
	return A /= B;
}
