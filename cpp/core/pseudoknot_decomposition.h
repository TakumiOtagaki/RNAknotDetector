#pragma once

#include <vector>

#include "entanglement.h"

namespace rna {

std::vector<BasePair> ExtractMainLayerFromBasePairs(
    const std::vector<BasePair> &base_pairs);

}  // namespace rna
