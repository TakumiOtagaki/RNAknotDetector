#include "geometry2d.h"

#include <algorithm>
#include <cmath>

namespace rna {
namespace {

double DistancePointSegmentSquared(const Vec2 &p, const Vec2 &a, const Vec2 &b) {
  double vx = b.x - a.x;
  double vy = b.y - a.y;
  double wx = p.x - a.x;
  double wy = p.y - a.y;
  double vv = vx * vx + vy * vy;
  if (vv <= 0.0) {
    double dx = p.x - a.x;
    double dy = p.y - a.y;
    return dx * dx + dy * dy;
  }
  double t = (wx * vx + wy * vy) / vv;
  if (t < 0.0) {
    double dx = p.x - a.x;
    double dy = p.y - a.y;
    return dx * dx + dy * dy;
  }
  if (t > 1.0) {
    double dx = p.x - b.x;
    double dy = p.y - b.y;
    return dx * dx + dy * dy;
  }
  double projx = a.x + t * vx;
  double projy = a.y + t * vy;
  double dx = p.x - projx;
  double dy = p.y - projy;
  return dx * dx + dy * dy;
}

double Cross2D(const Vec2 &a, const Vec2 &b, const Vec2 &c) {
  double abx = b.x - a.x;
  double aby = b.y - a.y;
  double acx = c.x - a.x;
  double acy = c.y - a.y;
  return abx * acy - aby * acx;
}

std::vector<Vec2> ConvexHull(std::vector<Vec2> points) {
  if (points.size() < 3) {
    return points;
  }
  std::sort(points.begin(), points.end(), [](const Vec2 &a, const Vec2 &b) {
    if (a.x == b.x) {
      return a.y < b.y;
    }
    return a.x < b.x;
  });
  std::vector<Vec2> hull;
  hull.reserve(points.size() * 2);
  for (const auto &p : points) {
    while (hull.size() >= 2) {
      size_t n = hull.size();
      if (Cross2D(hull[n - 2], hull[n - 1], p) <= 0.0) {
        hull.pop_back();
      } else {
        break;
      }
    }
    hull.push_back(p);
  }
  size_t lower_size = hull.size();
  for (size_t i = points.size(); i-- > 0;) {
    const auto &p = points[i];
    while (hull.size() > lower_size) {
      size_t n = hull.size();
      if (Cross2D(hull[n - 2], hull[n - 1], p) <= 0.0) {
        hull.pop_back();
      } else {
        break;
      }
    }
    hull.push_back(p);
  }
  if (!hull.empty()) {
    hull.pop_back();
  }
  return hull;
}

}  // namespace

Polygon2D ProjectPolygon(const std::vector<Vec3> &points, const Plane &plane) {
  Polygon2D poly;
  if (!plane.valid) {
    return poly;
  }
  std::vector<Vec2> vertices;
  vertices.reserve(points.size());
  for (const auto &p : points) {
    Vec3 d = Sub(p, plane.c);
    vertices.push_back(Vec2{Dot(d, plane.e1), Dot(d, plane.e2)});
  }
  if (vertices.size() < 3) {
    return poly;
  }
  poly.vertices = ConvexHull(std::move(vertices));
  poly.valid = poly.vertices.size() >= 3;
  return poly;
}

bool PointInPolygon2D(const Vec2 &q, const Polygon2D &poly, double eps_polygon) {
  if (!poly.valid || poly.vertices.size() < 3) {
    return false;
  }
  double eps2 = eps_polygon * eps_polygon;
  for (size_t i = 0; i < poly.vertices.size(); ++i) {
    const Vec2 &a = poly.vertices[i];
    const Vec2 &b = poly.vertices[(i + 1) % poly.vertices.size()];
    if (DistancePointSegmentSquared(q, a, b) <= eps2) {
      return true;
    }
  }
  bool inside = false;
  for (size_t i = 0, j = poly.vertices.size() - 1; i < poly.vertices.size(); j = i++) {
    const Vec2 &pi = poly.vertices[i];
    const Vec2 &pj = poly.vertices[j];
    bool intersect = ((pi.y > q.y) != (pj.y > q.y)) &&
                     (q.x < (pj.x - pi.x) * (q.y - pi.y) / (pj.y - pi.y + 1e-12) + pi.x);
    if (intersect) {
      inside = !inside;
    }
  }
  return inside;
}

}  // namespace rna
