# epsilon

This library can be used to compute derivatives of arbitrary order to machine precision.
The class `epsilon` implements Toeplitz matrices together with standard numberical operations on them.
```
double x = 1.23;
epsilon<3> x_(x);
auto y = f(x_);
```
Now y[0] == f(x), y[1] == f'(x), and y[2] == f''(x).
