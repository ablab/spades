// Require initializer lists, initializer list assignment (fails on
// Clang 3.0 even though initializer lists are otherwise supported),
// variadic templates and r-value references.
#include <initializer_list>

int  x { 42 };

template <typename... Args>
void
foo (int&&, Args&&... args)
{
  x = { 0 };
}
