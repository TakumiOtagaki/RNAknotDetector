#pragma once

#include <algorithm>
#include <utility>

#include "entanglement.h"

namespace rna {

inline std::pair<int, int> SortedPair(const BasePair &pair) {
  return {std::min(pair.i, pair.j), std::max(pair.i, pair.j)};
}

}  // namespace rna
