# epsilon

This library can be used to compute derivatives of arbitrary order to machine precision.
The class `epsilon` implements Toeplitz matrices together with standard numerical operations on them.

You must define functions that take generic arguments. For example:
```
template<class X>
X poly(X x)
{
    return x*x + 2.*x + 3.;
}
```
Now you can use `epsilon` to compute derivatives.
```
double x = 1.23;
epsilon<3> x_(x);
auto y = f(x_);
assert (y[0] == f(x));
assert (y[1] == 2*x + 2);
assert (y[2] == 2);
```

