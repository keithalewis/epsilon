#include <cassert>
#include "multi_index.h"

using namespace fms;

int test_multi_index()
{
	{
		multi_index i{ 1, 2, 3 };
		multi_index j{ 4, 5 };
		auto ij = i * j;
		assert(ij.size() == 3);
		assert(ij[0] = i[0] + j[0]);
		assert(ij[1] = i[1] + j[1]);
		assert(ij[2] = i[2]);
	}
	{
		multi_index i{ 1, 2, 3 };
		multi_index j{ 4, 5 };
		auto ij = i / j;
		assert(ij.size() == 3);
		assert(ij[0] = i[0] - j[0]);
		assert(ij[1] = i[1] - j[1]);
		assert(ij[2] = i[2]);

		auto ji = j / i;
		assert(ji.size() == 3);
		assert(ji[0] = j[0] - i[0]);
		assert(ji[1] = j[1] - i[1]);
		assert(ji[2] =      - i[2]);
	}

	return 0;
}
int test_multi_index_ = test_multi_index();

int test_multinomial()
{
	using i = multi_index;
	{
		multinomial a({
			{ i{ 1, 0 }, 2 }, // 2 x_0
			{ i{ 0, 1 }, 3 }, // + 3 x_1
			});
		auto a2 = a * a; // 4 x_0^2 + 9 x_1^2 + 12 x_0 x_1
		assert(a2.size() == 3);
		auto a20 = a2[i{ 2, 0 }];
		assert(a20 == 4);
		auto a02 = a2[i{ 0, 2 }];
		assert(a02 == 9);
		auto a11 = a2[i{ 1, 1 }];
		assert(a11 == 12);

		auto a0 = a2 - a2;
		for (const auto& [k, v] : a0) {
			assert(v == 0);
		}
	}
	{
		multinomial a({
			{ { i{ 0, 0 } }, 2 }, // 2
			{ { i{ 1, 0 } }, 3 }  // + 2 x_0
			});
		multinomial b({
			{ i{ 1, 0 }, 4 }, // 4 x_0
			{ i{ 0, 1 }, 5 }, // + 5 x_1
			});
		auto ab = a + b;
		assert(ab.size() == 3);
		auto a00 = ab[i{ 0, 0 }];
		assert(2 == a00);
		auto a10 = ab[i{ 1, 0 }];
		assert(7 == a10);
		auto a01 = ab[i{ 0, 1 }];
		assert(5 == a01);
	}

	return 0;
}
int test_multinomial_ = test_multinomial();