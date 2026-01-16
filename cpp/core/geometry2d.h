#pragma once

#include <vector>

#include "geometry3d.h"
#include "geometry_types.h"

namespace rna {

struct Polygon2D {
  std::vector<Vec2> vertices;
  bool valid = false;
};

Polygon2D ProjectPolygon(const std::vector<Vec3> &points, const Plane &plane);
bool PointInPolygon2D(const Vec2 &q, const Polygon2D &poly, double eps_polygon);

}  // namespace rna
