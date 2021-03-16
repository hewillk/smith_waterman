## Biomodern.SWAligner

`Biomodern.SWAligner` is modern C++ Smith-Waterman wrapper which using [SSW][SSW] Library as a underlying code base.

## Compilers
- GCC 10.2

## Example

```cpp
#include <iostream>
#include "istring.hpp"

int main() {
  using namespace biomodern;
  using namespace biomodern::utility;
  const auto read = "1102111000031323330032232203332323"_is;
  const auto ref = "2131100321101000010313231313001322323213"_is;
  const auto prof = SWAligner::get_profile(read, mat);
  const auto [score, ref_beg, ref_end, read_beg, read_end, cigar] = SWAligner::align(prof, ref);
  std::cout << "score: " << score << "\n";
  std::cout << "ref_beg: " << ref_beg << "\n";
  std::cout << "ref_end: " << ref_end << "\n";
  std::cout << "read_beg: " << read_beg << "\n";
  std::cout << "read_end: " << read_end << "\n";
  std::cout << "cigar: " << cigar << "\n";
}
```

[SSW]: https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library
