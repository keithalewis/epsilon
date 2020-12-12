#include "multi_epsilon.h"
#include <cassert>
#include <iostream>
using fms::multi_epsilon;
using std::cout;
using std::endl;

int mul_test() {
	//test + number
	multi_epsilon a(1,1);
	a += 1;
	assert(a(0, 0) == 1);
	assert(a(0, 1) == 0);
	assert(a(1, 1) == 1);

	//test + matrix
	a += a;
	assert(a(0, 0) == 2);
	assert(a(0, 1) == 0);
	assert(a(1, 1) == 2);

	//test / number
	a /= 2;
	assert(a(0, 0) == 1);
	assert(a(0, 1) == 0);
	assert(a(1, 1) == 1);

	//test - matrix
	a = 2 * a - a;
	assert(a(0, 0) == 1);
	assert(a(0, 1) == 0);
	assert(a(1, 1) == 1);

	//test - number
	a = 2 * a - 1;
	assert(a(0, 0) == 1);
	assert(a(0, 1) == 0);
	assert(a(1, 1) == 1);

	//test inverse
	multi_epsilon b({ 1,2,3,4 }, 1, 3);
	int t = 0;
	b = b.inverse();
	assert(b[0] == 1.0);
	assert(b[1] == -2.0);
	assert(b[2] == 1.0);
	assert(b[3] == 0);
	b = b.inverse();

	//test = and * matrix
	a = b * b;
	assert(a[0] == 1.0);
	assert(a[1] == 4.0);
	assert(a[2] == 10.0);
	assert(a[3] == 20.0);

	//test / matrix
	multi_epsilon c = a / b;
	t = 0;
	assert(c == b);

	//test ==
	assert(c == c);

	
	return 0;
}

template<class X>
X pp(X x, X y) {
	return x / y;
}

template<class X>
X dpp_dx([[maybe_unused]] X x, X y) {
	return (X)1.0 / y;
}

template<class X>
X dpp_dy(X x, X y) {
	return -x / (y * y);
}

template<class X>
X d2pp_dx_dy([[maybe_unused]] X x, X y) {
	return -(X)1.0 / (y * y);
}

int mul_test2() {
	double x = 1.0;
	double y = 2.0;
	auto epsilon = multi_epsilon::add_epsilon({ x,y }, 1);
	auto e_x = epsilon[0];
	auto e_y = epsilon[1];
	auto dpp = pp(e_x, e_y);
	//test pp|(x=1.0,y=2.0)
	auto index = multi_epsilon::rep({ 0,0 }, 1);
	assert(dpp(0, index) == pp(x, y));
	//test dpp/dx|(x=1.0,y=2.0)
	index = multi_epsilon::rep({ 1,0 }, 1);
	assert(dpp(0, index) == dpp_dx(x, y));
	//test dpp/dy|(x=1.0,y=2.0)
	index = multi_epsilon::rep({ 0,1 }, 1);
	assert(dpp(0, index) == dpp_dy(x, y));
	//test d^2pp/dxdy|(x=1.0,y=2.0)
	index = multi_epsilon::rep({ 1,1 }, 1);
	assert(dpp(0, index) == d2pp_dx_dy(x, y));
	return 0;
}

template<class X>
X q(const X& x, const X& y, const X& z) {
	return (x * x + y * y) / z;
}

template<class X>
X dq_dx([[maybe_unused]] const X& x, [[maybe_unused]] const X& y, const X& z) {
	return ((X)2 * x) / z;
}

template<class X>
X dq_dy([[maybe_unused]] const X& x, const X& y, const X& z) {
	return ((X)2 * y) / z;
}

template<class X>
X dq_dz(const X& x, const X& y, const X& z) {
	return -(x * x + y * y) / (z * z);
}

template<class X>
X d2q_dx_dy([[maybe_unused]] const X& x, [[maybe_unused]] const X& y, [[maybe_unused]] const X& z) {
	return (X)0;
}

template<class X>
X d2q_dx_dz(const X& x, [[maybe_unused]] const X& y, const X& z) {
	return -((X)2 * x) / (z * z);
}

template<class X>
X d3q_dx_dy_dz([[maybe_unused]] const X& x, [[maybe_unused]] const X& y, [[maybe_unused]] const X& z) {
	return (X)0;
}

int mul_test3() {
	size_t index;
	double x = 1.0;
	double y = 2.0;
	double z = 3.0;
	auto epsilon = multi_epsilon::add_epsilon({ x,y,z }, 1);
	auto e_x = epsilon[0];
	auto e_y = epsilon[1];
	auto e_z = epsilon[2];
	auto dq = q(e_x, e_y, e_z);
	//test q|(x=1.0,y=2.0,z=3.0)
	index = multi_epsilon::rep({ 0,0,0 }, 1);
	assert(fabs(dq(0, index) - q(x, y, z)) <= .001);
	//test dq/dx|(x=1.0,y=2.0,z=3.0)
	index = multi_epsilon::rep({ 1,0,0 }, 1);
	assert(dq(0, index) == dq_dx(x, y, z));
	//test dq/dy|(x=1.0,y=2.0,z=3.0)
	index = multi_epsilon::rep({ 0,1,0 }, 1);
	assert(dq(0, index) == dq_dy(x, y, z));
	//test dq/dz|(x=1.0,y=2.0,z=3.0)
	index = multi_epsilon::rep({ 0,0,1 }, 1);
	assert(dq(0, index) == dq_dz(x, y, z));
	//test d^2q/dxdy|(x=1.0,y=2.0,z=3.0)
	index = multi_epsilon::rep({ 1,1,0 }, 1);
	assert(dq(0, index) == d2q_dx_dy(x, y, z));
	//test d^2q/dxdz|(x=1.0,y=2.0,z=3.0)
	index = multi_epsilon::rep({ 1,0,1 }, 1);
	assert(dq(0, index) == d2q_dx_dz(x, y, z));
	//test d^3q/dxdydz|(x=1.0,y=2.0,z=3.0)
	index = multi_epsilon::rep({ 1,1,1 }, 1);
	assert(dq(0, index) == d3q_dx_dy_dz(x, y, z));
	cout << "Test suceessfully ends." << endl;
	return 0;
}
template<class X>
X p(const X& x)
{
	return 1 + 2 * x + 3 * x * x;
}
template<class X>
X dp_dx(const X& x)
{
	return (X)2 + (X)6 * x;
}
template<class X>
X d2p_dx2([[maybe_unused]] const X& x)
{
	return (X)6;
}
int mul_test4() {
	size_t index;
	double x = 1.0;
	auto epsilon = multi_epsilon::add_epsilon({ x }, 3);
	auto e_x = epsilon[0];
	auto dp = p(e_x);
	index = multi_epsilon::rep({ 0 }, 3);
	assert(dp(0, index) == p(x));
	index = multi_epsilon::rep({ 1 }, 3);
	assert(dp(0, index) == dp_dx(x));
	index = multi_epsilon::rep({ 2 }, 3);
	assert(dp(0, index) == d2p_dx2(x) / 2.0);
	index = multi_epsilon::rep({ 3 }, 3);
	assert(dp(0, index) == 0);
	return 0;
}

static int a = mul_test();
static int b = mul_test2();
static int c = mul_test3();
static int d = mul_test4();