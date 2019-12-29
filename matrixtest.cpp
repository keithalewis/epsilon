#include "TriangularMatrix.h"
#include <cassert>
#include <iostream>
using fms::TriangularMatrix;
using std::cout;
using std::endl;
int test() {
	//test + number
	TriangularMatrix a(2);
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
	TriangularMatrix b(4);
	int t = 0;
	for (int i = 0; i < 4; i++)
		for (int j = i; j < 4; j++)
		{
			b(i, j) = ++t;
		}
	b = b.inverse();
	assert(b(0, 0) == 1.0);
	assert(fabs(b(0, 1)+.4)<.0001);
	assert(fabs(b(0, 2)+0.075)<0.0001);
	assert(fabs(b(0, 3)+0.0525)<.0001);
	assert(fabs(b(1, 1)-.2)<0.0001);
	assert(fabs(b(1, 2)+.15)<0.0001);
	assert(fabs(b(1, 3)+.005)<0.0001);
	assert(fabs(b(2, 2)-.125)<.0001);
	assert(fabs(b(2, 3)+.1125)<.0001);
	assert(fabs(b(3, 3)-.1)<.0001);
	b = b.inverse();

	//test = and * matrix
	a = b * b;
	assert(a(0, 0) == 1.0);
	assert(a(0, 1) == 12);
	assert(a(0, 2) == 39);
	assert(a(0, 3) == 85);
	assert(a(1, 1) == 25);
	assert(fabs(a(1, 2) - 78.0) < 0.001);
	assert(fabs(a(1, 3) - 159)<0.001);
	assert(a(2, 2) == 64);
	assert(a(2, 3) == 162);
	assert(a(3, 3) == 100);

	//test / matrix
	TriangularMatrix c = a / b;
	t = 0;
	for (int i = 0; i < 4; i++) {
		for (int j = i; j < 4; j++)
			assert(fabs(c(i, j) - ++t)<0.001);
	}

	//test !=
	assert(c != b);

	//test ==
	assert(c == c);

	//test outer
	//  0 1   1 0
	//a=0 0 b=0 1
	a = TriangularMatrix(2);
	a(0, 1) = 1;
	b = TriangularMatrix(2);
	b += 1;
	//  	  0 0 1 0
	//		  0 0 0 1
	//		  0 0 0 0
	//c=a ¡Á b=0 0 0 0 
	c = a.outer(b);
	assert(c.Size == 4);
	assert(c(0, 2) == 1);
	assert(c(1 , 3) == 1);
	//  	  0 1 0 0
	//		  0 0 0 0
	//		  0 0 0 1
	//c=b ¡Á a=0 0 0 0 
	c = b.outer(a);
	assert(c(0, 1) == 1);
	assert(c(2, 3) == 1);

	//test multi_epsilon
	std::vector<TriangularMatrix>epsilon_list = TriangularMatrix::multi_epsilon(3, 1);
	assert(epsilon_list[0] == TriangularMatrix({0,0,0,0,1,0,0,0,
												0,0,0,0,1,0,0,
												0,0,0,0,1,0,
												0,0,0,0,1,
												0,0,0,0,
												0,0,0,
												0,0,
												0}));
	assert(epsilon_list[1] == TriangularMatrix({0,0,1,0,0,0,0,0,
												0,0,1,0,0,0,0,
												0,0,0,0,0,0,
												0,0,0,0,0,
												0,0,1,0,
												0,0,1,
												0,0,
												0}));
	assert(epsilon_list[2] == TriangularMatrix({0,1,0,0,0,0,0,0,
												0,0,0,0,0,0,0,
												0,1,0,0,0,0,
												0,0,0,0,0,
												0,1,0,0,
												0,0,0,
												0,1,
												0}));
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
	return -x / (y*y);
}

template<class X>
X d2pp_dx_dy([[maybe_unused]] X x, X y) {
	return -(X)1.0 / (y * y);
}

int test2() {
	size_t index = TriangularMatrix::rep(std::vector<size_t>({ 1,1 }),1);
	assert(index == 3);
	double x = 1.0;
	double y = 2.0;
	auto epsilon = TriangularMatrix::add_epsilon({ x,y }, 1);
	auto e_x = epsilon[0];
	auto e_y = epsilon[1];
	auto dpp=pp(e_x, e_y);
	//test pp|(x=1.0,y=2.0)
	index =TriangularMatrix::rep({ 0,0 }, 1);
	assert(dpp(0, index) == pp(x,y));
	//test dpp/dx|(x=1.0,y=2.0)
	index = TriangularMatrix::rep({ 1,0 }, 1);
	assert(dpp(0, index) == dpp_dx(x,y));
	//test dpp/dy|(x=1.0,y=2.0)
	index = TriangularMatrix::rep({ 0,1 }, 1);
	assert(dpp(0, index) == dpp_dy(x,y));
	//test d^2pp/dxdy|(x=1.0,y=2.0)
	index = TriangularMatrix::rep({ 1,1 }, 1);
	assert(dpp(0, index) == d2pp_dx_dy(x,y));
	return 0;
}

template<class X>
X q(const X& x, const X& y, const X& z) {
	return (x*x+y*y)/z;
}

template<class X>
X dq_dx([[maybe_unused]] const X& x, [[maybe_unused]] const X& y, const X& z) {
	return ((X)2*x)/z;
}

template<class X>
X dq_dy([[maybe_unused]] const X& x, const X& y, const X& z) {
	return ((X)2*y)/z;
}

template<class X>
X dq_dz(const X& x, const X& y, const X& z) {
	return -(x * x + y * y) / (z*z);
}

template<class X>
X d2q_dx_dy([[maybe_unused]] const X& x, [[maybe_unused]] const X& y, [[maybe_unused]] const X& z) {
	return (X)0;
}

template<class X>
X d2q_dx_dz(const X& x, [[maybe_unused]] const X& y, const X& z) {
	return -((X)2 * x) / (z*z);
}

template<class X>
X d3q_dx_dy_dz([[maybe_unused]] const X& x, [[maybe_unused]] const X& y, [[maybe_unused]] const X& z) {
	return (X)0;
}

int test3() {
	size_t index;
	double x = 1.0;
	double y = 2.0;
	double z = 3.0;
	auto epsilon = TriangularMatrix::add_epsilon({ x,y,z }, 1);
	auto e_x = epsilon[0];
	auto e_y = epsilon[1];
	auto e_z = epsilon[2];
	auto dq = q(e_x, e_y,e_z);
	//test q|(x=1.0,y=2.0,z=3.0)
	index = TriangularMatrix::rep({ 0,0,0 }, 1);
	assert(fabs(dq(0, index) - q(x, y,z))<.001);
	//test dq/dx|(x=1.0,y=2.0,z=3.0)
	index = TriangularMatrix::rep({ 1,0,0 }, 1);
	assert(dq(0, index) == dq_dx(x, y,z));
	//test dq/dy|(x=1.0,y=2.0,z=3.0)
	index = TriangularMatrix::rep({ 0,1,0 }, 1);
	assert(dq(0, index) == dq_dy(x, y,z));
	//test dq/dz|(x=1.0,y=2.0,z=3.0)
	index = TriangularMatrix::rep({ 0,0,1 }, 1);
	assert(dq(0, index) == dq_dz(x, y, z));
	//test d^2q/dxdy|(x=1.0,y=2.0,z=3.0)
	index = TriangularMatrix::rep({ 1,1,0 }, 1);
	assert(dq(0, index) == d2q_dx_dy(x, y,z));
	//test d^2q/dxdz|(x=1.0,y=2.0,z=3.0)
	index = TriangularMatrix::rep({ 1,0,1 }, 1);
	assert(dq(0, index) == d2q_dx_dz(x, y, z));
	//test d^3q/dxdydz|(x=1.0,y=2.0,z=3.0)
	index = TriangularMatrix::rep({ 1,1,1 }, 1);
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
int test4() {
	size_t index;
	double x = 1.0;
	auto epsilon = TriangularMatrix::add_epsilon({ x }, 3);
	auto e_x = epsilon[0];
	auto dp = p(e_x);
	index= TriangularMatrix::rep({ 0 }, 3);
	assert(dp(0, index) == p(x));
	index = TriangularMatrix::rep({ 1 }, 3);
	assert(dp(0, index) == dp_dx(x));
	index = TriangularMatrix::rep({ 2 }, 3);
	assert(dp(0, index) == d2p_dx2(x)/2.0);
	index = TriangularMatrix::rep({ 3 }, 3);
	assert(dp(0, index) == 0);
	return 0;
}

static int a = test();
static int b = test2();
static int c = test3();
static int d = test4();