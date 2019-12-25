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
	cout << "Test suceessfully ends."<<endl;
	return 0;
}

int a = test();