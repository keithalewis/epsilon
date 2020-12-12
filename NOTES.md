# NOTES

Actually use epsilon.

```
template<class X>
inline X f(X x) { ... return some function of x; }

double x = 1.23;
auto y = f(x + epsilon<3>{});
double y2 = y[2]; // f''(x)
```