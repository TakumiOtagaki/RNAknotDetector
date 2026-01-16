#include "entanglement.h"
#include "coord_utils.h"
#include "loop_utils.h"
#include "pair_utils.h"
#include "pseudoknot_decomposition.h"

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <unordered_set>

namespace rna {
namespace {

// 例えば hairpin closing (i, j) の polygon に関して、線分 (i, j) などはスキップするべき。そういう操作。
std::vector<int> BuildSkipResidues(const Loop &loop) {
  std::vector<int> skip;
  if (loop.closing_pairs.empty()) {
    return skip;
  }
  if (loop.kind == LoopKind::kHairpin) {
    auto [i, j] = SortedPair(loop.closing_pairs[0]);
    for (int k = i; k <= j; ++k) {
      skip.push_back(k);
    }
    return skip;
  }
  if (loop.kind == LoopKind::kInternal) {
    if (loop.closing_pairs.size() < 2) {
      auto [i, j] = SortedPair(loop.closing_pairs[0]);
      for (int k = i; k <= j; ++k) {
        skip.push_back(k);
      }
      return skip;
    }
    auto [i, j] = SortedPair(loop.closing_pairs[0]);
    auto [k, l] = SortedPair(loop.closing_pairs[1]);
    for (int idx = i; idx <= k; ++idx) {
      skip.push_back(idx);
    }
    for (int idx = l; idx <= j; ++idx) {
      skip.push_back(idx);
    }
    return skip;
  }
  if (loop.kind == LoopKind::kMulti) {
    int min_res = std::numeric_limits<int>::max();
    int max_res = std::numeric_limits<int>::min();
    for (const auto &pair : loop.closing_pairs) {
      auto [i, j] = SortedPair(pair);
      min_res = std::min(min_res, i);
      max_res = std::max(max_res, j);
      skip.push_back(i);
      skip.push_back(j);
    }
    if (min_res <= max_res) {
      for (int idx = min_res; idx <= max_res; ++idx) {
        skip.push_back(idx);
      }
    }
    return skip;
  }
  return skip;
}


// HitKey は (loop_id, segment_id) の組を一意に識別するキーを作るための関数です。
// EvaluateEntanglement 内で unordered_set<int64_t> hit_keys に入れて、同じループ×線分の重複カウントを防ぐ用途に使っています。
// 具体的にはこう使っています（entanglement.cpp）:
// - 交差が見つかったら HitKey(loop_id, segment_id) を作成
// - そのキーが未登録なら hits に追加、既にあればスキップ
// - つまり “ユニークなヒット数 = K” を保証する仕組みです
//
// もし loop_id や segment_id が 32bit を超えることがあるなら衝突するので、その場合は別のキー形式に変える必要があります。
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

  for (const auto &loop : loops) {
    Surface surface;
    surface.loop_id = loop.id;
    surface.kind = loop.kind;
    surface.closing_pairs = loop.closing_pairs;
    surface.skip_residues = BuildSkipResidues(loop);

    if (loop.closing_pairs.empty()) {
      std::cerr << "[debug] loop=" << loop.id << " kind=" << static_cast<int>(loop.kind)
                << " closing_pairs=0\n";
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
        pairs.push_back(SortedPair(pair));
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
      // add_range を最小スコープで定義。閉区間 [start, end] の残基を追加。
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
        auto [i, j] = SortedPair(loop.closing_pairs[0]);
        add_range(i, j);
      } else if (loop.kind == LoopKind::kInternal) {
        auto [i, j] = SortedPair(loop.closing_pairs[0]);
        if (loop.closing_pairs.size() >= 2) {
          auto [h, l] = SortedPair(loop.closing_pairs[1]);
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
        std::cerr << "[warning] loop=" << loop.id
                  << " unknown loop kind for boundary construction\n";
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
      auto log_with_loop = [&](const char *status) {
        std::cerr << "[debug] loop=" << surface.loop_id
                  << " type=" << static_cast<int>(surface.kind)
                  << " pairs=";
        for (size_t idx = 0; idx < surface.closing_pairs.size(); ++idx) {
          const auto &bp = surface.closing_pairs[idx];
          std::cerr << "(" << bp.i << "," << bp.j << ")";
          if (idx + 1 < surface.closing_pairs.size()) {
            std::cerr << ",";
          }
        }
        std::cerr << " segment=(" << i << "," << i + 1 << ") " << status;
      };

      if (skip_mask[i] || skip_mask[i + 1]) {
        if (watch_segment) {
          log_with_loop("skipped_by_mask\n");
        }
        continue;
      }
      candidate_segments++;
      Vec3 intersection;
      if (!SegmentPlaneIntersection(segment.a, segment.b, surface.plane, options.eps_plane,
                                    &intersection)) {
        if (surface.kind == LoopKind::kMulti && i == 46) {
          double d_a = Dot(Sub(segment.a, surface.plane.c), surface.plane.n_hat);
          double d_b = Dot(Sub(segment.b, surface.plane.c), surface.plane.n_hat);
          std::cerr << "[debug] multi loop=" << surface.loop_id
                    << " segment=46 plane_miss d_a=" << d_a
                    << " d_b=" << d_b << "\n";
        }
        if (watch_segment) {
          log_with_loop("plane_miss\n");
        }
        continue;
      }
      if (surface.kind == LoopKind::kMulti && i == 46) {
        std::cerr << "[debug] multi loop=" << surface.loop_id
                  << " segment=46 plane_hit"
                  << " point=(" << intersection.x << "," << intersection.y << ","
                  << intersection.z << ")\n";
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
        log_with_loop("plane_hit");
        std::cerr << " in_polygon=" << in_poly
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
