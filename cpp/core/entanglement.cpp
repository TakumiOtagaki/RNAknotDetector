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

bool IsTargetMultiloop_debug(const std::vector<BasePair> &pairs) {
  if (pairs.size() != 3) {
    return false;
  }
  std::vector<std::pair<int, int>> sorted;
  sorted.reserve(pairs.size());
  for (const auto &pair : pairs) {
    sorted.push_back(SortedPair(pair));
  }
  std::sort(sorted.begin(), sorted.end());
  return sorted[0] == std::make_pair(63, 121) &&
         sorted[1] == std::make_pair(70, 96) &&
         sorted[2] == std::make_pair(98, 105);
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
    bool watch_target_multiloop =
        (surface.kind == LoopKind::kMulti &&
         IsTargetMultiloop_debug(surface.closing_pairs));
    std::vector<char> skip_mask(map.n_res + 1, 0);
    for (int idx : surface.skip_residues) {
      if (idx > 0 && idx <= map.n_res) {
        skip_mask[idx] = 1;
      }
    }
    int candidate_segments = 0;
    int plane_hits = 0;
    int triangle_hits = 0;
    bool use_triangles = !surface.triangles.empty();
    if (!use_triangles && (!surface.plane.valid || !surface.polygon.valid)) {
      continue;
    }
    if (watch_target_multiloop) {
      std::cerr << "[debug] target_multiloop triangles=" << surface.triangles.size()
                << "\n";
    }
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
        if (watch_target_multiloop && i == 46) {
          std::cerr << "[debug] target_multiloop loop=" << surface.loop_id
                    << " segment=46 skipped_by_mask\n";
        }
        if (watch_segment) {
          log_with_loop("skipped_by_mask\n");
        }
        continue;
      }
      if (watch_target_multiloop && i == 46) {
        std::cerr << "[debug] target_multiloop segment46 a=("
                  << segment.a.x << "," << segment.a.y << "," << segment.a.z
                  << ") b=(" << segment.b.x << "," << segment.b.y << ","
                  << segment.b.z << ")\n";
        if (!surface.triangles.empty()) {
          const auto &tri = surface.triangles.front();
          std::cerr << "[debug] target_multiloop tri0 a=(" << tri.a.x << ","
                    << tri.a.y << "," << tri.a.z << ") b=(" << tri.b.x << ","
                    << tri.b.y << "," << tri.b.z << ") c=(" << tri.c.x << ","
                    << tri.c.y << "," << tri.c.z << ")\n";
        }
      }
      candidate_segments++;
      Vec3 intersection;
      bool hit = false;
      if (use_triangles) {
        int triangle_tests = 0;
        for (const auto &tri : surface.triangles) {
          triangle_tests++;
          if (SegmentIntersectsTriangle(segment.a, segment.b, tri,
                                        options.eps_triangle, &intersection)) {
            hit = true;
            triangle_hits++;
            break;
          }
        }
        if (watch_target_multiloop && i == 46) {
          std::cerr << "[debug] target_multiloop loop=" << surface.loop_id
                    << " segment=46 triangle_"
                    << (hit ? "hit" : "miss")
                    << " tests=" << triangle_tests << "\n";
        }
      } else {
        if (!SegmentPlaneIntersection(segment.a, segment.b, surface.plane, options.eps_plane,
                                      &intersection)) {
          if (watch_target_multiloop && i == 46) {
            double d_a = Dot(Sub(segment.a, surface.plane.c), surface.plane.n_hat);
            double d_b = Dot(Sub(segment.b, surface.plane.c), surface.plane.n_hat);
            std::cerr << "[debug] target_multiloop loop=" << surface.loop_id
                      << " segment=46 plane_miss d_a=" << d_a
                      << " d_b=" << d_b << "\n";
          }
          if (watch_segment) {
            log_with_loop("plane_miss\n");
          }
          continue;
        }
        if (watch_target_multiloop && i == 46) {
          std::cerr << "[debug] target_multiloop loop=" << surface.loop_id
                    << " segment=46 plane_hit"
                    << " point=(" << intersection.x << "," << intersection.y << ","
                    << intersection.z << ")\n";
        }
        plane_hits++;
        Vec3 d = Sub(intersection, surface.plane.c);
        Vec2 q{Dot(d, surface.plane.e1), Dot(d, surface.plane.e2)};
        bool in_poly = PointInPolygon2D(q, surface.polygon, options.eps_polygon);
      if (watch_target_multiloop && i == 46) {
        std::cerr << "[debug] target_multiloop loop=" << surface.loop_id
                  << " segment=46 in_polygon=" << in_poly
                  << " q=(" << q.x << "," << q.y << ")\n";
      }
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
        hit = true;
      }
      if (!hit) {
        continue;
      }
      int64_t key = HitKey(surface.loop_id, segment.id);
      if (hit_keys.insert(key).second) {
        result.hits.push_back(HitInfo{surface.loop_id, segment.id, intersection});
      }
    }
    (void)candidate_segments;
    (void)plane_hits;
    (void)triangle_hits;
  }
  result.K = static_cast<int>(result.hits.size());
  return result;
}

}  // namespace rna
