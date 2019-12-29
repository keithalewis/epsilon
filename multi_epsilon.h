#pragma once
#include <cstring>
#include <cassert>
#include <iostream>
#include <functional>
#include "TriangularMatrix.h"
namespace fms {
	class multi_epsilon{
	public:
		multi_epsilon(const size_t N) : N(N) {};
		~multi_epsilon() {};

		void calc(const std::function<fms::TriangularMatrix(...)> &f, std::vector<double>& x) {
			std::vector<fms::TriangularMatrix> epsilon = fms::TriangularMatrix::multi_epsilon(x.size(), N);
			for (size_t i = 0; i < x.size(); i++) {
				mat.emplace_back(x[i] + epsilon[i]);
			}
			
		}
		
		double get_derivative(std::vector<size_t> order) {

		};
	private:
		size_t N;
		std::vector<double> x;
		std::vector<fms::TriangularMatrix> mat;
		fms::TriangularMatrix result;
	};
}