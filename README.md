# epsilon

This library can be used to compute derivatives of arbitrary order to machine precision.
The class `epsilon` implements Toeplitz matrices together with standard numerical operations on them.

If ε is not 0 and ε<sup>2</sup> = 0 then for any twice differentiable function f,
f(x + ε) = f(x) + f'(x) ε + 1/2 f''(x) ε<sup>2</sup>.

You must define functions that take generic arguments. For example:
```
template<class X>
X f(X x)
{
    return x*x + 2*x + 3;
}
```
Now you can use `epsilon` to compute derivatives.
```
double x = 1.23;
auto y = f(epsilon<3>(x));
assert (y[0] == f(x));
assert (y[1] == 2*x + 2);
assert (y[2] == 2);
```

