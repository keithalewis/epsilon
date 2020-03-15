#include <cassert>
class AA {
  public:
    int a;
    AA(int a) : a(a){};
    AA &operator-()
    {
        a *= -1;
        return *this;
    }
};
int test_a_trivial_thing()
{
    AA aa(3);
    auto b = -aa;
    auto c = -aa;
    assert(b.a != c.a);
    return 1;
}
static int test_result = test_a_trivial_thing();