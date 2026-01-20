#include "surface_builder.h"

#include <algorithm>
#include <cmath>
#include <iostream>

#include "coord_utils.h"
#include "geometry2d.h"
#include "geometry3d.h"
#include "loop_utils.h"
#include "pair_utils.h"

namespace rna {
namespace {

struct OrderedPoint {
  Vec3 p;
  Vec2 q;
  double angle;
};

std::vector<OrderedPoint> OrderPointsByAngle(const std::vector<Vec3> &points,
                                             const Plane &plane) {
  std::vector<OrderedPoint> ordered;
  if (!plane.valid || points.size() < 3) {
    return ordered;
  }
  ordered.reserve(points.size());
  Vec2 center{0.0, 0.0};
  for (const auto &p : points) {
    Vec3 d = Sub(p, plane.c);
    double x = Dot(d, plane.e1);
    double y = Dot(d, plane.e2);
    center.x += x;
    center.y += y;
  }
  center.x /= static_cast<double>(points.size());
  center.y /= static_cast<double>(points.size());
  for (const auto &p : points) {
    Vec3 d = Sub(p, plane.c);
    double x = Dot(d, plane.e1);
    double y = Dot(d, plane.e2);
    Vec3 proj = Add(plane.c, Add(Scale(plane.e1, x), Scale(plane.e2, y)));
    double angle = std::atan2(y - center.y, x - center.x);
    ordered.push_back(OrderedPoint{proj, Vec2{x, y}, angle});
  }
  std::sort(ordered.begin(), ordered.end(),
            [](const OrderedPoint &a, const OrderedPoint &b) {
              return a.angle < b.angle;
            });
  return ordered;
}

std::vector<int> BuildBoundaryIndices(const Loop &loop, int n_res) {
  std::vector<int> boundary_indices;
  boundary_indices.reserve(loop.boundary_residues.size() +
                           loop.closing_pairs.size() * 2);
  std::vector<char> seen(n_res + 1, 0);
  auto add_index = [&](int res_index) {
    if (res_index <= 0 || res_index > n_res) {
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
    if (loop.closing_pairs.size() == 3) {
      std::vector<std::pair<int, int>> check = pairs;
      if (check.size() == 3 &&
          check[0] == std::make_pair(63, 121) &&
          check[1] == std::make_pair(70, 96) &&
          check[2] == std::make_pair(98, 105)) {
        std::cerr << "[debug] target_multiloop boundary_indices:";
        for (int idx : boundary_indices) {
          std::cerr << " " << idx;
        }
        std::cerr << "\n";
      }
    }
    return boundary_indices;
  }

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
    for (int res_index : loop.boundary_residues) {
      add_index(res_index);
    }
    for (const auto &pair : loop.closing_pairs) {
      add_index(pair.i);
      add_index(pair.j);
    }
  }
  return boundary_indices;
}

}  // namespace

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

    std::vector<int> boundary_indices = BuildBoundaryIndices(loop, map.n_res);

    std::vector<Vec3> boundary_points;
    boundary_points.reserve(boundary_indices.size());
    for (int res_index : boundary_indices) {
      if (!map.has_coord[res_index]) {
        continue;
      }
      boundary_points.push_back(map.coords[res_index]);
    }
    if (options.surface_mode == SurfaceMode::kBestFitPlane) {
      surface.plane = FitPlane(boundary_points, options.eps_collinear);
      surface.polygon = ProjectPolygon(boundary_points, surface.plane);
    } else {
      surface.plane = FitPlane(boundary_points, options.eps_collinear);
      surface.polygon.valid = false;
      surface.polygon.vertices.clear();
      if (surface.plane.valid && boundary_points.size() >= 3) {
        std::vector<OrderedPoint> ordered =
            OrderPointsByAngle(boundary_points, surface.plane);
        surface.polygon.vertices.clear();
        surface.polygon.vertices.reserve(ordered.size());
        for (const auto &point : ordered) {
          surface.polygon.vertices.push_back(point.q);
        }
        surface.polygon.valid = surface.polygon.vertices.size() >= 3;
        const Vec3 &p0 = ordered[0].p;
        for (size_t i = 1; i + 1 < ordered.size(); ++i) {
          Triangle tri{p0, ordered[i].p, ordered[i + 1].p};
          Vec3 ab = Sub(tri.b, tri.a);
          Vec3 ac = Sub(tri.c, tri.a);
          double area = Norm(Cross(ab, ac));
          if (area <= options.eps_collinear) {
            continue;
          }
          surface.triangles.push_back(tri);
        }
      }
    }
    surfaces.push_back(std::move(surface));
  }
  return surfaces;
}

}  // namespace rna
