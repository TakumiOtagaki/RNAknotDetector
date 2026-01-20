#include "surface_builder.h"

#include <algorithm>
#include <array>
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

double SignedArea2D(const std::vector<Vec2> &poly) {
  if (poly.size() < 3) {
    return 0.0;
  }
  double area = 0.0;
  for (size_t i = 0; i < poly.size(); ++i) {
    const Vec2 &a = poly[i];
    const Vec2 &b = poly[(i + 1) % poly.size()];
    area += a.x * b.y - a.y * b.x;
  }
  return 0.5 * area;
}

double Cross2D(const Vec2 &a, const Vec2 &b, const Vec2 &c) {
  double abx = b.x - a.x;
  double aby = b.y - a.y;
  double acx = c.x - a.x;
  double acy = c.y - a.y;
  return abx * acy - aby * acx;
}

bool PointInTriangle2D(const Vec2 &p,
                       const Vec2 &a,
                       const Vec2 &b,
                       const Vec2 &c,
                       double eps) {
  double c1 = Cross2D(a, b, p);
  double c2 = Cross2D(b, c, p);
  double c3 = Cross2D(c, a, p);
  bool has_neg = (c1 < -eps) || (c2 < -eps) || (c3 < -eps);
  bool has_pos = (c1 > eps) || (c2 > eps) || (c3 > eps);
  return !(has_neg && has_pos);
}

std::vector<std::array<int, 3>> EarClipTriangulate(
    const std::vector<Vec2> &poly,
    double eps) {
  std::vector<std::array<int, 3>> tris;
  if (poly.size() < 3) {
    return tris;
  }
  double area = SignedArea2D(poly);
  if (std::abs(area) <= eps) {
    return tris;
  }
  int orientation = (area > 0.0) ? 1 : -1;
  std::vector<int> indices;
  indices.reserve(poly.size());
  for (size_t i = 0; i < poly.size(); ++i) {
    indices.push_back(static_cast<int>(i));
  }

  int guard = 0;
  while (indices.size() > 3 && guard < 10000) {
    bool ear_found = false;
    size_t n = indices.size();
    for (size_t i = 0; i < n; ++i) {
      int i_prev = indices[(i + n - 1) % n];
      int i_curr = indices[i];
      int i_next = indices[(i + 1) % n];
      const Vec2 &a = poly[i_prev];
      const Vec2 &b = poly[i_curr];
      const Vec2 &c = poly[i_next];
      double cross = Cross2D(a, b, c);
      if (orientation * cross <= eps) {
        continue;
      }
      bool has_inside = false;
      for (size_t k = 0; k < n; ++k) {
        int idx = indices[k];
        if (idx == i_prev || idx == i_curr || idx == i_next) {
          continue;
        }
        if (PointInTriangle2D(poly[idx], a, b, c, eps)) {
          has_inside = true;
          break;
        }
      }
      if (has_inside) {
        continue;
      }
      tris.push_back({i_prev, i_curr, i_next});
      indices.erase(indices.begin() + static_cast<long>(i));
      ear_found = true;
      break;
    }
    if (!ear_found) {
      tris.clear();
      return tris;
    }
    guard++;
  }
  if (indices.size() == 3) {
    tris.push_back({indices[0], indices[1], indices[2]});
  }
  return tris;
}

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
    for (int res_index : loop.boundary_residues) {
      add_index(res_index);
    }
    for (const auto &pair : loop.closing_pairs) {
      add_index(pair.i);
      add_index(pair.j);
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
        std::vector<Vec2> poly2d;
        std::vector<Vec3> poly3d;
        poly2d.reserve(boundary_points.size());
        poly3d.reserve(boundary_points.size());
        for (const auto &p : boundary_points) {
          Vec3 d = Sub(p, surface.plane.c);
          double x = Dot(d, surface.plane.e1);
          double y = Dot(d, surface.plane.e2);
          Vec3 proj = Add(surface.plane.c,
                          Add(Scale(surface.plane.e1, x),
                              Scale(surface.plane.e2, y)));
          poly2d.push_back(Vec2{x, y});
          poly3d.push_back(proj);
        }
        surface.polygon.vertices = poly2d;
        surface.polygon.valid = surface.polygon.vertices.size() >= 3;
        std::vector<std::array<int, 3>> tris =
            EarClipTriangulate(poly2d, 1e-12);
        for (const auto &tri : tris) {
          Triangle t{poly3d[tri[0]], poly3d[tri[1]], poly3d[tri[2]]};
          Vec3 ab = Sub(t.b, t.a);
          Vec3 ac = Sub(t.c, t.a);
          double area = Norm(Cross(ab, ac));
          if (area <= options.eps_collinear) {
            continue;
          }
          surface.triangles.push_back(t);
        }
      }
    }
    surfaces.push_back(std::move(surface));
  }
  return surfaces;
}

}  // namespace rna
