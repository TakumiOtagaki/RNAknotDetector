#pragma once

#include <vector>

#include "entanglement.h"

namespace rna {

std::vector<Surface> BuildSurfaces(const std::vector<ResidueCoord> &coords,
                                   const std::vector<Loop> &loops,
                                   const SurfaceBuildOptions &options);

}  // namespace rna
