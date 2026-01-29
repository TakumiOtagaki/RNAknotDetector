#pragma once

#include <algorithm>
#include <cmath>
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
  int res_a = 0;
  int res_b = 0;
  AtomKind atom_a = AtomKind::kSingle;
  AtomKind atom_b = AtomKind::kSingle;
  Vec3 a;
  Vec3 b;

  Segment() = default;
  Segment(int id_in,
          int res_a_in,
          int res_b_in,
          AtomKind atom_a_in,
          AtomKind atom_b_in,
          const Vec3 &a_in,
          const Vec3 &b_in)
      : id(id_in),
        res_a(res_a_in),
        res_b(res_b_in),
        atom_a(atom_a_in),
        atom_b(atom_b_in),
        a(a_in),
        b(b_in) {}
};

struct PolylinePoint {
  int res_index = 0;
  AtomKind atom_kind = AtomKind::kSingle;
  Vec3 point;

  PolylinePoint() = default;
  PolylinePoint(int res_index_in, AtomKind atom_kind_in, const Vec3 &point_in)
      : res_index(res_index_in), atom_kind(atom_kind_in), point(point_in) {}
};

inline bool IsFiniteCoord(const Vec3 &v) {
  return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z);
}

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
    const Vec3 &v = res.atoms[atom_index];
    if (!IsFiniteCoord(v)) {
      continue;
    }
    map.coords[res.res_index] = v;
    map.has_coord[res.res_index] = 1;
  }
  return map;
}

inline std::vector<PolylinePoint> BuildPolylinePoints(
    const std::vector<ResidueCoord> &coords,
    int atom_index_p,
    int atom_index_c4) {
  std::vector<PolylinePoint> points;
  CoordMap map_p = BuildCoordMap(coords, atom_index_p);
  CoordMap map_c4 = BuildCoordMap(coords, atom_index_c4);
  int n_res = std::max(map_p.n_res, map_c4.n_res);
  points.reserve(n_res * 2);
  for (int i = 1; i <= n_res; ++i) {
    if (i <= map_p.n_res && map_p.has_coord[i]) {
      points.push_back(PolylinePoint{i, AtomKind::kP, map_p.coords[i]});
    }
    if (i <= map_c4.n_res && map_c4.has_coord[i]) {
      points.push_back(PolylinePoint{i, AtomKind::kC4, map_c4.coords[i]});
    }
  }
  return points;
}

inline std::vector<Segment> BuildSegmentsFromPolyline(
    const std::vector<PolylinePoint> &points) {
  std::vector<Segment> segments;
  if (points.size() < 2) {
    return segments;
  }
  segments.reserve(points.size() - 1);
  for (size_t i = 0; i + 1 < points.size(); ++i) {
    const auto &a = points[i];
    const auto &b = points[i + 1];
    segments.push_back(
        Segment{static_cast<int>(i + 1),
                a.res_index,
                b.res_index,
                a.atom_kind,
                b.atom_kind,
                a.point,
                b.point});
  }
  return segments;
}

}  // namespace rna
