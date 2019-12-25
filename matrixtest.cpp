#include "TriangularMatrix.h"
#include <cassert>
#include <iostream>
using fms::TriangularMatrix;
using std::cout;
using std::endl;
int test() {
	TriangularMatrix a(2);
	a += 1;
	assert(a(0, 0) == 1);
	assert(a(0, 1) == 0);
	assert(a(1, 1) == 1);
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
	for (int i = 0; i < 4; i++) {
		for (int j = i; j < 4; j++)
			cout << a(i, j) << ' ';
		cout << endl;
	}
	cout << "Test suceessfully ends.";
	return 0;
}

int a = test();