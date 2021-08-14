// exp.t.cpp - text exp
#include <cassert>
#include "exp.h"

using namespace fms;

int fms_exp_test = []() {
	{
		double x = 1;
		epsilon<3, double> x_(x);
		auto y_ = exp(x_);
		assert(y_[0] == exp(1));
		assert(y_[1] == exp(1));
		assert(y_[2] == exp(1));
	}

	return 0;
}();