#pragma once
#include <cstring>
#include <cassert>
#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>
using std::cout;
using std::endl;
namespace fms {
	class TriangularMatrix {
	public:
		size_t  Size;
		//a_0 a_1 a_2 a_3
		// 0  a_4 a_5 a_6
		// 0   0  a_7 a_8
		// 0   0   0  a_9
		//matrix element in memory, 1D array
		TriangularMatrix(size_t m = 0) {			
			Size = m;
			this->m_lpBuf = new double[(m + 1) * m / 2];
			std::memset(m_lpBuf, 0, (m + 1) * m / 2 * sizeof(double));
		};
		TriangularMatrix(const TriangularMatrix& rhs) {			
			this->Size = rhs.Size;
			this->m_lpBuf = new double[(Size + 1) * Size / 2];
			for (int i = 0; i < (Size * (Size + 1) / 2); ++i)
				this->m_lpBuf[i] = rhs.m_lpBuf[i];
		};
		
		TriangularMatrix(TriangularMatrix&& rhs) noexcept {			
			this->Size = rhs.Size;
			this->m_lpBuf = rhs.m_lpBuf;
			rhs.m_lpBuf = nullptr;
		};

		~TriangularMatrix() {
			Size = 0;
			delete[] this->m_lpBuf;
			this->m_lpBuf = nullptr;
		};


		TriangularMatrix& operator=(const TriangularMatrix& rhs) {
			this->Size = rhs.Size;
			delete[] this->m_lpBuf;
			this->m_lpBuf = new double[(Size + 1) * Size / 2];
			for (int i = 0; i < (Size * (Size + 1) / 2); ++i)
				this->m_lpBuf[i] = rhs.m_lpBuf[i];
			return *this;
		};

		
		TriangularMatrix& operator=(TriangularMatrix&& rhs) noexcept {
			this->Size = rhs.Size;
			delete [] this->m_lpBuf;
			this->m_lpBuf = rhs.m_lpBuf;
			rhs.m_lpBuf = nullptr;
			return *this;
		};
		
		double& operator ()(size_t i, size_t j) const
		{
			assert(i < Size);
			assert(j < Size);
			assert(j >= i);
			//if (i > j) return 0;//visiting lower triangle element
			return *(m_lpBuf + Size * i - i * (i - 1) / 2 + j - i);
		}

		bool operator == (const TriangularMatrix& rhs) const{
			if (this->Size != rhs.Size) return false;
			for (size_t i = 0; i < Size; i++)
				if (this->m_lpBuf[i] != rhs.m_lpBuf[i])
					return false;
			return true;			
		}

		bool operator !=(const TriangularMatrix& rhs) const {
			return !operator==(rhs);
		}

		TriangularMatrix& operator += (const TriangularMatrix& rhs)
		{
			if (!rhs.m_lpBuf) return *this;
			assert(rhs.Size == this->Size);
			for (int i = 0; i < (Size + 1) * Size / 2; i++)
				this->m_lpBuf[i] += rhs.m_lpBuf[i];

			return *this;
		}

		//this+rhs*I
		TriangularMatrix& operator += (const double& rhs)
		{
			for (int i = 0; i < Size; i++)
				operator()(i, i) += rhs;
			return *this;
		}

		TriangularMatrix& operator -= (const TriangularMatrix& rhs)
		{
			if (!rhs.m_lpBuf) return *this;
			assert(rhs.Size == this->Size);
			for (int i = 0; i < (Size + 1) * Size / 2; i++)
				this->m_lpBuf[i] -= rhs.m_lpBuf[i];

			return *this;
		}

		//this-rhs*I
		TriangularMatrix& operator -= (const double& rhs)
		{
			for (int i = 0; i < Size; i++)
				operator()(i, i) -= rhs;
			return *this;
		}

		//matrix multiplication
		TriangularMatrix& operator *= (const TriangularMatrix& rhs)
		{
			//if (!rhs.m_lpBuf) return *this;
			assert(rhs.Size == this->Size);
			double* temp = new double[this->Size];
			for (size_t i = 0; i < Size; i++) {
				std::memset(temp, 0, sizeof(double) * this->Size);
				for (size_t j = i; j < Size; j++) {
					double t = 0;
					for (size_t k = i; k <= j; k++)
						t += operator()(i, k) * rhs(k, j);
					temp[j] = t;
				}
				for (size_t j = i; j < this->Size; j++)
					operator()(i, j) = temp[j];
			}
			delete[] temp;
			return *this;
		}

		//this*rhs
		TriangularMatrix& operator *= (const double& rhs)
		{
			for (int i = 0; i < (Size + 1) * Size / 2; i++)
				m_lpBuf[i] *= rhs;
			return *this;
		}

		//matrix division
		TriangularMatrix& operator /= (const TriangularMatrix& rhs)
		{
			//if (!rhs.m_lpBuf) return *this;
			assert(rhs.Size == this->Size);
			TriangularMatrix rhs_inv = rhs.inverse();
			double* temp = new double[this->Size];
			for (size_t i = 0; i < Size; i++) {
				std::memset(temp, 0, sizeof(double) * this->Size);
				for (size_t j = i; j < Size; j++) {
					double t = 0;
					for (size_t k = i; k <= j; k++)
						t += operator()(i, k) * rhs_inv(k, j);
					temp[j] = t;
				}
				for (size_t j = i; j < this->Size; j++)
					operator()(i, j) = temp[j];
			}
			delete[] temp;
			return *this;
		}

		//this/rhs
		TriangularMatrix& operator /= (const double& rhs)
		{
			for (int i = 0; i < (Size + 1) * Size / 2; i++)
				m_lpBuf[i] /= rhs;
			return *this;
		}

		//-1*this
		TriangularMatrix& operator -()
		{
			operator*=(-1);
			return *this;
		}

		//inverse(this)
		//Gauss reduction method
		TriangularMatrix inverse() const {
			TriangularMatrix identity(Size);
			if (!Size) return *this;
			identity += 1;
			for (size_t i = Size - 1; i >= 0 && i <= Size - 1; i--) {//consider the ith row
				for (size_t j = i; j < Size; j++)
					identity(i, j) /= operator()(i, i);
				operator()(i, i) = 1;
				if (i == 0) continue;
				for (size_t j = i - 1; j >= 0 && j <= i - 1; j--) {
					for (size_t k = i; k < Size; k++)
						identity(j, k) -= operator()(j, i) * identity(i, k);
					operator()(j, i) = 0;
				}
			}
			//for (size_t i = 0; i < (Size + 1) * Size / 2; i++)
			//    this->m_lpBuf[i] = identity.m_lpBuf[i];
			return identity;
		}

		//direct product
		TriangularMatrix outer(const TriangularMatrix & rhs) const{
			TriangularMatrix result(this->Size * rhs.Size);
			for (size_t i = 0; i < Size; i++) {
				for (size_t j = i; j < Size; j++) {
					for (size_t r_i=0; r_i < rhs.Size; r_i++) {
						for (size_t r_j=r_i; r_j < rhs.Size; r_j++) {
							if (i * rhs.Size + r_i <= j * rhs.Size + r_j) {
								result.operator()(i * rhs.Size + r_i, j * rhs.Size + r_j) = operator()(i, j) * rhs(r_i, r_j);
							}
						}
					}
				}
			}
			return result;
		}

		//static method simple epsilon
		//for order (n-1) th derivative
		//dimension of matrix=n
		// 0 1 0
		// 0 0 1
		// 0 0 0 e.g. n=3, pow=1
		// 0 0 1
		// 0 0 0
		// 0 0 0 e.g. n=3, pow=2
		static TriangularMatrix simple_epsilon(size_t n, size_t pow=1) {
			TriangularMatrix result(n);
			for (size_t i = 0; i+pow < n; i++)
				result(i, i + pow) = 1;
			return result;
		}

		//static method identity
		//(n-1) th order derivative
		//dimension of matrix=n
		static TriangularMatrix identity(size_t n) {
			TriangularMatrix result(n);
			result += 1;
			return result;
		}

		//static method multi_epsilon
		//M variable
		//N th order derivative
		static std::vector<TriangularMatrix> multi_epsilon(size_t M, size_t N) {
			std::vector<TriangularMatrix> result;
			for (size_t i = 0; i < M; i++) {
				TriangularMatrix t = identity(i*(N+1));
				if (i != 0)
					t = t.outer(simple_epsilon(N + 1));
				else
					t = simple_epsilon(N + 1);
				if (i!=M-1)
					result.emplace_back(t.outer(identity((M - i - 1)*(N+1))));
				else
					result.emplace_back(t);
			}
			return result;
		}


		//get specific epsilon representation
		//e.g. order=[1,1]
		//N=2 total order
		//epsilon^1.outer(epsilon^1)= 0 0 0 1
		//							  0 0 0 0
		//							  0 0 0 0
		//							  0 0 0 0
		//return 3
		static size_t rep(const std::vector<size_t>& order, size_t N){
			assert(*std::max_element(order.begin(), order.end()) <= N);
			TriangularMatrix t(1);
			t += 1;
			for (auto& it : order) {
				if (it != 0) t = t.outer(simple_epsilon(N + 1, it));
				else t = t.outer(identity(N + 1));
			}
			int pos = 0;
			while (t.m_lpBuf[pos] != 1 && pos<t.Size) pos++;
			return pos;
		}



		//print current matrix
		void print() {
			for (size_t i = 0; i < Size; i++) {
				for (size_t j = 0; j < Size; j++)
					if (j < i) std::cout << 0 << ' ';
					else std::cout << operator()(i, j) << ' ';
				std::cout << std::endl;
			}
		}
	private:
		double* m_lpBuf; // pointer to data
		
	};
}
inline fms::TriangularMatrix operator + (fms::TriangularMatrix A, const fms::TriangularMatrix& B)
{
	return A += B;
}
inline fms::TriangularMatrix operator + (fms::TriangularMatrix A, const double& B)
{
	return A += B;
}
inline fms::TriangularMatrix operator + (const double& B, fms::TriangularMatrix A)
{
	return A += B;
}
inline fms::TriangularMatrix operator - (fms::TriangularMatrix A, const fms::TriangularMatrix& B)
{
	return A -= B;
}
inline fms::TriangularMatrix operator - (fms::TriangularMatrix A, const double& B)
{
	return A -= B;
}
inline fms::TriangularMatrix operator - (const double& B, fms::TriangularMatrix A)
{
	return A -= B;
}
inline fms::TriangularMatrix operator * (fms::TriangularMatrix A, const fms::TriangularMatrix& B)
{	
	A *= B;
	return A;
}
inline fms::TriangularMatrix operator * (fms::TriangularMatrix A, const double& B)
{
	A *= B;
	return A;
}
inline fms::TriangularMatrix operator * (const double& B, fms::TriangularMatrix A)
{
	return A *= B;
}
inline fms::TriangularMatrix operator / (fms::TriangularMatrix A, const fms::TriangularMatrix& B)
{
	return A /= B;
}
inline fms::TriangularMatrix operator / (fms::TriangularMatrix A, const double& B)
{
	return A /= B;
}
inline fms::TriangularMatrix operator / (const double& B, fms::TriangularMatrix A)
{
	return A /= B;
}
inline std::vector<fms::TriangularMatrix> multi_epsilon(const std::vector<double>& x, size_t N) {
	std::vector<fms::TriangularMatrix> result;
	std::vector<fms::TriangularMatrix> epsilon = fms::TriangularMatrix::multi_epsilon(x.size(), N);
	for (size_t i = 0; i < x.size();i++) {
		result.emplace_back(x[i]+epsilon[i]);
	}
	return result;
}