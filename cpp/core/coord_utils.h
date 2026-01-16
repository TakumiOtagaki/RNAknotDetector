#pragma once

#include <algorithm>
#include <vector>

#include "entanglement.h"

namespace rna {

struct CoordMap {
  int n_res = 0;
  std::vector<Vec3> coords;
  std::vector<char> has_coord;
};

struct Segment {
  int id = 0;
  Vec3 a;
  Vec3 b;
};

inline CoordMap BuildCoordMap(const std::vector<ResidueCoord> &coords, int atom_index) {
  int max_index = 0;
  for (const auto &res : coords) {
    max_index = std::max(max_index, res.res_index);
  }
  CoordMap map;
  map.n_res = max_index;
  map.coords.resize(max_index + 1);
  map.has_coord.assign(max_index + 1, 0);
  for (const auto &res : coords) {
    if (res.res_index <= 0 || res.res_index > max_index) {
      continue;
    }
    if (atom_index < 0 || atom_index >= static_cast<int>(res.atoms.size())) {
      continue;
    }
    map.coords[res.res_index] = res.atoms[atom_index];
    map.has_coord[res.res_index] = 1;
  }
  return map;
}

}  // namespace rna
