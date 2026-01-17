#pragma once

#include <vector>

#include "geometry_types.h"

namespace rna {

struct Plane {
  Vec3 c;
  Vec3 n_hat;
  Vec3 e1;
  Vec3 e2;
  bool valid = false;
};

struct Triangle {
  Vec3 a;
  Vec3 b;
  Vec3 c;
};

Vec3 Add(const Vec3 &a, const Vec3 &b);
Vec3 Sub(const Vec3 &a, const Vec3 &b);
Vec3 Scale(const Vec3 &v, double s);
double Dot(const Vec3 &a, const Vec3 &b);
Vec3 Cross(const Vec3 &a, const Vec3 &b);
double Norm(const Vec3 &v);
Vec3 Normalize(const Vec3 &v);

Plane FitPlane(const std::vector<Vec3> &points, double eps_collinear);
bool SegmentPlaneIntersection(const Vec3 &a,
                              const Vec3 &b,
                              const Plane &plane,
                              double eps_plane,
                              Vec3 *out_point);
bool SegmentIntersectsTriangle(const Vec3 &a,
                               const Vec3 &b,
                               const Triangle &tri,
                               double eps,
                               Vec3 *out_point);

}  // namespace rna
