#pragma once
#include <cstring>
#include <cassert>
#include <iostream>
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
			rhs.m_lpBuf = NULL;
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
			delete[] this->m_lpBuf;
			this->m_lpBuf = rhs.m_lpBuf;
			rhs.m_lpBuf = NULL;
			return *this;
		};
		~TriangularMatrix() {
			Size = 0;
			delete[] this->m_lpBuf;
			this->m_lpBuf = NULL;
		};
		double& operator ()(size_t i, size_t j) const
		{
			assert(i < Size);
			assert(j >= i);
			//if (i > j) return 0;//visiting lower triangle element
			return *(m_lpBuf + Size * i - i * (i - 1) / 2 + j - i);
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
