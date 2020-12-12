#include <cassert>
#include "bell.h"

using namespace fms;

template<class X>
int test_bell()
{
	X B[6];
	{
		X x[] = { 1, 1, 1, 1, 1, 1 };
		Bell(5, x, B);
		assert(B[0] == 1);
		assert(B[1] == 1);
		assert(B[2] == 2);
		assert(B[3] == 5);
		assert(B[4] == 15);
		assert(B[5] == 52);
	}
	{
		X x[] = { 2, 3, 4 };
		Bell(3, x, B);
		assert(B[0] == 1);
		assert(B[1] == x[0]);
		assert(B[2] == x[0] * x[0] + x[1]);
		assert(B[3] == x[0] * x[0] * x[0] + 3 * x[0] * x[1] + x[2]);
	}

	return 0;
}

int bell = test_bell<double>();

