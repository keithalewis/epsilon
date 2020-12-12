![C/C++ CI](https://github.com/keithalewis/epsilon/workflows/C/C++%20CI/badge.svg)

# epsilon

This library can be used to compute derivatives of any order to machine precision.

If ε is not 0 and ε<sup>2</sup> = 0 then for any sufficiently differentiable function _f_,
_f_(_x_ + ε) = _f_(_x_) + _f_'(_x_) ε + _f_''(_x_) ε<sup>2</sup>/2 + &middot; &middot; &middot; =  _f_(_x_) + _f_'(_x_) ε.
The coefficient of ε is the derivative of _f_.

For example, if _f_(_x_) = _x_<sup>2</sup> 
then (_x_ + ε)<sup>2</sup> = _x_<sup>2</sup> + 2 _x_ ε + ε<sup>2</sup> = _x_<sup>2</sup> + 2 _x_ ε,
so the derivative of _x_<sup>2</sup> is 2 _x_. No need to take limits of difference
quotients.

Of course ε cannot be a real number, but the 2 x 2 matrix ε = [0 1; 0 0] satisfies these conditions.
This is a special case of a Toeplitz matrix. The class `epsilon` implements the
arithmetic datatype required to compute derivatives to any order.

In order to use this library you must define functions that take generic arguments. For example:
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
auto y = f(epsilon<3>(x, 1));
assert (y[0] == f(x));
assert (y[1] == 2*x + 2); // f'(x)
assert (y[2] == 2);       // f''(x)
```

