#include "TriangularMatrix.h"
#include <cassert>
#include <iostream>
using fms::TriangularMatrix;
using std::cout;
using std::endl;
int test() {
	TriangularMatrix a(1);
	assert(a(0,0) == 0);
	a = TriangularMatrix(4);
	assert(a(0, 0) == 0);
	assert(a(0, 1) == 0);
	assert(a(1, 1) == 0);
	TriangularMatrix b(4);
	int t = 0;
	for (int i = 0; i < 4; i++)
		for (int j = i; j < 4; j++)
		{
			b(i, j) = ++t;
			a(i, j) = t;
		}
	a *=b;
	b.inverse();
	for (int i = 0; i < 4; i++) {
		for (int j = i; j < 4; j++)
			cout << b(i, j) << ' ';
		cout << endl;
	}
	cout << "Test suceessfully ends.";
	return 0;
}

int a = test();