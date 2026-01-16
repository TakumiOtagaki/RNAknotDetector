#include "entanglement.h"
#include "pseudoknot_decomposition.h"

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <unordered_set>

namespace rna {
namespace {

std::vector<int> BuildPairMap(const std::vector<BasePair> &base_pairs, int n_res) {
  std::vector<int> pair_map(n_res + 1, 0);
  for (const auto &bp : base_pairs) {
    if (bp.i <= 0 || bp.j <= 0 || bp.i > n_res || bp.j > n_res) {
      throw std::invalid_argument("Base pair index out of range");
    }
    if (bp.i == bp.j) {
      throw std::invalid_argument("Base pair cannot be self-paired");
    }
    int i = std::min(bp.i, bp.j);
    int j = std::max(bp.i, bp.j);
    if (pair_map[i] != 0 || pair_map[j] != 0) {
      throw std::invalid_argument("Residue paired multiple times");
    }
    pair_map[i] = j;
    pair_map[j] = i;
  }
  return pair_map;
}

bool IsPaired(const std::vector<int> &pair_map, int idx) {
  return pair_map[idx] != 0;
}


std::vector<int> CollectUnpaired(const std::vector<int> &pair_map, int start, int end) {
  std::vector<int> residues;
  for (int k = start; k <= end; ++k) {
    if (!IsPaired(pair_map, k)) {
      residues.push_back(k);
    }
  }
  return residues;
}

// Find immediate child base pairs inside (i, j).
// A child pair is the first paired region encountered when scanning the
// interval; nested pairs inside that region are ignored.
// Used to count how many stems close the loop bounded by (i, j).
std::vector<BasePair> FindChildPairs(const std::vector<int> &pair_map, int i, int j) {
  std::vector<BasePair> child_pairs;
  int depth = 0;
  for (int idx = i + 1; idx <= j - 1; ++idx) {
    if (!IsPaired(pair_map, idx)) {
      continue;
    }
    int partner = pair_map[idx];
    if (idx < partner) {
      if (depth == 0) {
        child_pairs.push_back(BasePair{idx, partner, BasePair::Type::kUnclassified});
      }
      depth++;
    } else if (idx > partner) {
      depth--;
    }
  }
  return child_pairs;
}

// Classify loop by counting immediate child pairs within (i, j).
// closing_pairs includes the outer pair (i, j) and each immediate child pair.
// boundary collects unpaired residues that lie on the loop boundary.
// Rules: 0 child -> hairpin, 1 child -> internal/bulge/stacking, 2+ -> multi.
LoopKind ClassifyLoop(const std::vector<int> &pair_map,
                      int i,
                      int j,
                      std::vector<int> *boundary,
                      std::vector<BasePair> *closing_pairs) {
  closing_pairs->clear();
  closing_pairs->push_back(BasePair{i, j, BasePair::Type::kUnclassified});

  std::vector<BasePair> child_pairs = FindChildPairs(pair_map, i, j);
  closing_pairs->insert(closing_pairs->end(), child_pairs.begin(), child_pairs.end());

  if (child_pairs.empty()) {
    *boundary = CollectUnpaired(pair_map, i + 1, j - 1);
    return LoopKind::kHairpin;
  }
  if (child_pairs.size() == 1) {
    int k = std::min(child_pairs[0].i, child_pairs[0].j);
    int l = std::max(child_pairs[0].i, child_pairs[0].j);
    auto left = CollectUnpaired(pair_map, i + 1, k - 1);
    auto right = CollectUnpaired(pair_map, l + 1, j - 1);
    boundary->clear();
    boundary->reserve(left.size() + right.size());
    boundary->insert(boundary->end(), left.begin(), left.end());
    boundary->insert(boundary->end(), right.begin(), right.end());
    return LoopKind::kInternal;
  }

  *boundary = CollectUnpaired(pair_map, i + 1, j - 1);
  return LoopKind::kMulti;
}

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

CoordMap BuildCoordMap(const std::vector<ResidueCoord> &coords, int atom_index) {
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

std::vector<int> BuildSkipResidues(const Loop &loop) {
  std::vector<int> skip;
  if (loop.closing_pairs.empty()) {
    return skip;
  }
  if (loop.kind == LoopKind::kHairpin) {
    int i = std::min(loop.closing_pairs[0].i, loop.closing_pairs[0].j);
    int j = std::max(loop.closing_pairs[0].i, loop.closing_pairs[0].j);
    for (int k = i; k <= j; ++k) {
      skip.push_back(k);
    }
    return skip;
  }
  if (loop.kind == LoopKind::kInternal) {
    if (loop.closing_pairs.size() < 2) {
      int i = std::min(loop.closing_pairs[0].i, loop.closing_pairs[0].j);
      int j = std::max(loop.closing_pairs[0].i, loop.closing_pairs[0].j);
      for (int k = i; k <= j; ++k) {
        skip.push_back(k);
      }
      return skip;
    }
    int i = std::min(loop.closing_pairs[0].i, loop.closing_pairs[0].j);
    int j = std::max(loop.closing_pairs[0].i, loop.closing_pairs[0].j);
    int k = std::min(loop.closing_pairs[1].i, loop.closing_pairs[1].j);
    int l = std::max(loop.closing_pairs[1].i, loop.closing_pairs[1].j);
    for (int idx = i; idx <= k; ++idx) {
      skip.push_back(idx);
    }
    for (int idx = l; idx <= j; ++idx) {
      skip.push_back(idx);
    }
    return skip;
  }
  return skip;
}

int64_t HitKey(int loop_id, int segment_id) {
  return (static_cast<int64_t>(loop_id) << 32) | static_cast<uint32_t>(segment_id);
}

}  // namespace


// Build closed elements (loops) from base pairs.
// Each loop corresponds to an outer closing pair (i, j) and immediate child pairs.
// boundary_residues holds unpaired residues on the loop boundary (minimal set).
// Pseudoknots are not supported.
std::vector<Loop> BuildLoops(const std::vector<BasePair> &base_pairs,
                             int n_res,
                             const LoopBuildOptions &options) {
  if (n_res <= 0) {
    throw std::invalid_argument("n_res must be positive");
  }
  // Placeholder: assumes pseudoknot-free input; no validation for crossing pairs.
  std::vector<BasePair> filtered_pairs = base_pairs;
  if (options.main_layer_only) {
    filtered_pairs = ExtractMainLayer(base_pairs);
  }
  std::vector<int> pair_map = BuildPairMap(filtered_pairs, n_res);

  std::vector<Loop> loops;
  int loop_id = 1;
  for (int i = 1; i <= n_res; ++i) {
    int j = pair_map[i];
    if (j == 0 || i > j) {
      continue;
    }
    std::vector<int> boundary;
    std::vector<BasePair> closing_pairs;
    LoopKind kind = ClassifyLoop(pair_map, i, j, &boundary, &closing_pairs);
    if (kind == LoopKind::kMulti && !options.include_multi) {
      continue;
    }
    Loop loop;
    loop.id = loop_id++;
    loop.kind = kind;
    loop.closing_pairs = std::move(closing_pairs);
    loop.boundary_residues = std::move(boundary);
    loops.push_back(std::move(loop));
  }
  return loops;
}

std::vector<BasePair> ExtractMainLayer(const std::vector<BasePair> &base_pairs) {
  if (base_pairs.empty()) {
    return {};
  }
  return ExtractMainLayerFromBasePairs(base_pairs);
}

std::vector<Surface> BuildSurfaces(const std::vector<ResidueCoord> &coords,
                                   const std::vector<Loop> &loops,
                                   const SurfaceBuildOptions &options) {
  CoordMap map = BuildCoordMap(coords, options.atom_index);
  std::vector<Surface> surfaces;
  surfaces.reserve(loops.size());
  std::cerr << "[debug] BuildSurfaces: loops=" << loops.size()
            << " n_res=" << map.n_res << "\n";
  for (const auto &loop : loops) {
    Surface surface;
    surface.loop_id = loop.id;
    surface.kind = loop.kind;
    surface.skip_residues = BuildSkipResidues(loop);

    if (loop.closing_pairs.empty()) {
      std::cerr << "[debug] loop=" << loop.id << " kind=" << static_cast<int>(loop.kind)
                << " closing_pairs=0\n";
    } else {
      std::cerr << "[debug] loop=" << loop.id << " kind=" << static_cast<int>(loop.kind)
                << " closing_pairs=" << loop.closing_pairs.size() << " first_pair=("
                << loop.closing_pairs[0].i << "," << loop.closing_pairs[0].j << ")\n";
    }

    std::vector<int> boundary_indices;
    boundary_indices.reserve(loop.boundary_residues.size() +
                             loop.closing_pairs.size() * 2);
    std::vector<char> seen(map.n_res + 1, 0);
    auto add_index = [&](int res_index) {
      if (res_index <= 0 || res_index > map.n_res) {
        return;
      }
      if (seen[res_index]) {
        return;
      }
      seen[res_index] = 1;
      boundary_indices.push_back(res_index);
    };
    if (loop.kind == LoopKind::kMulti) {
      std::vector<std::pair<int, int>> pairs;
      pairs.reserve(loop.closing_pairs.size());
      for (const auto &pair : loop.closing_pairs) {
        int i = std::min(pair.i, pair.j);
        int j = std::max(pair.i, pair.j);
        pairs.emplace_back(i, j);
      }
      std::sort(pairs.begin(), pairs.end());
      if (!pairs.empty()) {
        int l = pairs.front().first;
        int i_branch = 0;
        int j_branch = 0;
        for (const auto &pair : pairs) {
          if (pair.first > l) {
            i_branch = pair.first;
            j_branch = pair.second;
            break;
          }
        }
        if (i_branch > 0) {
          for (int idx = l; idx <= i_branch - 1; ++idx) {
            add_index(idx);
          }
          add_index(i_branch);
          add_index(j_branch);
        } else {
          add_index(l);
          add_index(pairs.front().second);
        }
      }
    } else {
      auto add_range = [&](int start, int end) {
        if (start > end) {
          return;
        }
        for (int idx = start; idx <= end; ++idx) {
          add_index(idx);
        }
      };
      if (loop.closing_pairs.empty()) {
        for (int res_index : loop.boundary_residues) {
          add_index(res_index);
        }
      } else if (loop.kind == LoopKind::kHairpin) {
        int i = std::min(loop.closing_pairs[0].i, loop.closing_pairs[0].j);
        int j = std::max(loop.closing_pairs[0].i, loop.closing_pairs[0].j);
        add_range(i, j);
      } else if (loop.kind == LoopKind::kInternal) {
        int i = std::min(loop.closing_pairs[0].i, loop.closing_pairs[0].j);
        int j = std::max(loop.closing_pairs[0].i, loop.closing_pairs[0].j);
        if (loop.closing_pairs.size() >= 2) {
          int h = std::min(loop.closing_pairs[1].i, loop.closing_pairs[1].j);
          int l = std::max(loop.closing_pairs[1].i, loop.closing_pairs[1].j);
          add_range(i, h - 1);
          add_index(h);
          add_index(l);
          add_range(l + 1, j - 1);
          add_index(i);
          add_index(j);
        } else {
          add_range(i, j);
        }
      } else {
        for (int res_index : loop.boundary_residues) {
          add_index(res_index);
        }
        for (const auto &pair : loop.closing_pairs) {
          add_index(pair.i);
          add_index(pair.j);
        }
      }
    }

    std::vector<Vec3> boundary_points;
    boundary_points.reserve(boundary_indices.size());
    for (int res_index : boundary_indices) {
      if (!map.has_coord[res_index]) {
        continue;
      }
      boundary_points.push_back(map.coords[res_index]);
    }
    std::cerr << "[debug] loop=" << loop.id
              << " boundary_indices=" << boundary_indices.size()
              << " boundary_points=" << boundary_points.size()
              << " skip_residues=" << surface.skip_residues.size() << "\n";
    surface.plane = FitPlane(boundary_points, options.eps_collinear);
    surface.polygon = ProjectPolygon(boundary_points, surface.plane);
    std::cerr << "[debug] loop=" << loop.id
              << " plane_valid=" << surface.plane.valid
              << " polygon_valid=" << surface.polygon.valid << "\n";
    surfaces.push_back(std::move(surface));
  }
  return surfaces;
}

Result EvaluateEntanglement(const std::vector<ResidueCoord> &coords,
                            const std::vector<Surface> &surfaces,
                            const EvaluateOptions &options) {
  Result result;
  CoordMap map = BuildCoordMap(coords, options.atom_index);
  if (map.n_res <= 1) {
    return result;
  }
  std::vector<Segment> segments;
  segments.reserve(map.n_res);
  for (int i = 1; i < map.n_res; ++i) {
    if (!map.has_coord[i] || !map.has_coord[i + 1]) {
      continue;
    }
    segments.push_back(Segment{i, map.coords[i], map.coords[i + 1]});
  }
  std::cerr << "[debug] EvaluateEntanglement: surfaces=" << surfaces.size()
            << " segments=" << segments.size() << "\n";
  const std::unordered_set<int> debug_segments = {46, 89, 143};

  std::unordered_set<int64_t> hit_keys;
  for (const auto &surface : surfaces) {
    if (!surface.plane.valid || !surface.polygon.valid) {
      std::cerr << "[debug] loop=" << surface.loop_id
                << " skipped: invalid surface\n";
      continue;
    }
    std::vector<char> skip_mask(map.n_res + 1, 0);
    for (int idx : surface.skip_residues) {
      if (idx > 0 && idx <= map.n_res) {
        skip_mask[idx] = 1;
      }
    }
    int candidate_segments = 0;
    int plane_hits = 0;
    for (const auto &segment : segments) {
      int i = segment.id;
      bool watch_segment = debug_segments.count(i) > 0;
      if (skip_mask[i] || skip_mask[i + 1]) {
        if (watch_segment) {
          std::cerr << "[debug] loop=" << surface.loop_id
                    << " segment=" << i << " skipped_by_mask\n";
        }
        continue;
      }
      candidate_segments++;
      Vec3 intersection;
      if (!SegmentPlaneIntersection(segment.a, segment.b, surface.plane, options.eps_plane,
                                    &intersection)) {
        if (watch_segment) {
          std::cerr << "[debug] loop=" << surface.loop_id
                    << " segment=" << i << " plane_miss\n";
        }
        continue;
      }
      plane_hits++;
      Vec3 d = Sub(intersection, surface.plane.c);
      Vec2 q{Dot(d, surface.plane.e1), Dot(d, surface.plane.e2)};
      bool in_poly = PointInPolygon2D(q, surface.polygon, options.eps_polygon);
      if (watch_segment) {
        double min_x = 0.0;
        double min_y = 0.0;
        double max_x = 0.0;
        double max_y = 0.0;
        bool bbox_init = false;
        for (const auto &v : surface.polygon.vertices) {
          if (!bbox_init) {
            min_x = max_x = v.x;
            min_y = max_y = v.y;
            bbox_init = true;
            continue;
          }
          min_x = std::min(min_x, v.x);
          min_y = std::min(min_y, v.y);
          max_x = std::max(max_x, v.x);
          max_y = std::max(max_y, v.y);
        }
        std::cerr << "[debug] loop=" << surface.loop_id
                  << " segment=" << i << " plane_hit"
                  << " in_polygon=" << in_poly
                  << " q=(" << q.x << "," << q.y << ")"
                  << " poly_n=" << surface.polygon.vertices.size();
        if (bbox_init) {
          std::cerr << " poly_bbox=[(" << min_x << "," << min_y << "),("
                    << max_x << "," << max_y << ")]";
        }
        std::cerr << "\n";
      }
      if (!in_poly) {
        continue;
      }
      int64_t key = HitKey(surface.loop_id, segment.id);
      if (hit_keys.insert(key).second) {
        result.hits.push_back(HitInfo{surface.loop_id, segment.id, intersection});
      }
    }
    std::cerr << "[debug] loop=" << surface.loop_id
              << " candidates=" << candidate_segments
              << " plane_hits=" << plane_hits
              << " hits=" << result.hits.size() << "\n";
  }
  result.K = static_cast<int>(result.hits.size());
  return result;
}

}  // namespace rna
