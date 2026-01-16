#include "geometry3d.h"

#include <cmath>

namespace rna {

Vec3 Add(const Vec3 &a, const Vec3 &b) {
  return Vec3{a.x + b.x, a.y + b.y, a.z + b.z};
}

Vec3 Sub(const Vec3 &a, const Vec3 &b) {
  return Vec3{a.x - b.x, a.y - b.y, a.z - b.z};
}

Vec3 Scale(const Vec3 &v, double s) {
  return Vec3{v.x * s, v.y * s, v.z * s};
}

double Dot(const Vec3 &a, const Vec3 &b) {
  return a.x * b.x + a.y * b.y + a.z * b.z;
}

Vec3 Cross(const Vec3 &a, const Vec3 &b) {
  return Vec3{
      a.y * b.z - a.z * b.y,
      a.z * b.x - a.x * b.z,
      a.x * b.y - a.y * b.x,
  };
}

double Norm(const Vec3 &v) {
  return std::sqrt(Dot(v, v));
}

Vec3 Normalize(const Vec3 &v) {
  double n = Norm(v);
  if (n <= 0.0) {
    return Vec3{0.0, 0.0, 0.0};
  }
  return Scale(v, 1.0 / n);
}

namespace {

void JacobiEigenSymmetric3x3(double a[3][3],
                             double eigenvalues[3],
                             double eigenvectors[3][3]) {
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      eigenvectors[i][j] = (i == j) ? 1.0 : 0.0;
    }
  }
  for (int iter = 0; iter < 50; ++iter) {
    int p = 0;
    int q = 1;
    double max_offdiag = std::abs(a[p][q]);
    for (int i = 0; i < 3; ++i) {
      for (int j = i + 1; j < 3; ++j) {
        double val = std::abs(a[i][j]);
        if (val > max_offdiag) {
          max_offdiag = val;
          p = i;
          q = j;
        }
      }
    }
    if (max_offdiag < 1e-12) {
      break;
    }
    double phi = 0.5 * std::atan2(2.0 * a[p][q], a[q][q] - a[p][p]);
    double c = std::cos(phi);
    double s = std::sin(phi);

    double app = c * c * a[p][p] - 2.0 * s * c * a[p][q] + s * s * a[q][q];
    double aqq = s * s * a[p][p] + 2.0 * s * c * a[p][q] + c * c * a[q][q];
    a[p][p] = app;
    a[q][q] = aqq;
    a[p][q] = 0.0;
    a[q][p] = 0.0;

    for (int k = 0; k < 3; ++k) {
      if (k == p || k == q) {
        continue;
      }
      double akp = c * a[k][p] - s * a[k][q];
      double akq = s * a[k][p] + c * a[k][q];
      a[k][p] = akp;
      a[p][k] = akp;
      a[k][q] = akq;
      a[q][k] = akq;
    }

    for (int k = 0; k < 3; ++k) {
      double vip = c * eigenvectors[k][p] - s * eigenvectors[k][q];
      double viq = s * eigenvectors[k][p] + c * eigenvectors[k][q];
      eigenvectors[k][p] = vip;
      eigenvectors[k][q] = viq;
    }
  }
  for (int i = 0; i < 3; ++i) {
    eigenvalues[i] = a[i][i];
  }
}

}  // namespace

Plane FitPlane(const std::vector<Vec3> &points, double eps_collinear) {
  Plane plane;
  if (points.size() < 3) {
    return plane;
  }
  Vec3 c{0.0, 0.0, 0.0};
  for (const auto &p : points) {
    c = Add(c, p);
  }
  c = Scale(c, 1.0 / static_cast<double>(points.size()));

  double cov[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
  for (const auto &p : points) {
    Vec3 d = Sub(p, c);
    cov[0][0] += d.x * d.x;
    cov[0][1] += d.x * d.y;
    cov[0][2] += d.x * d.z;
    cov[1][0] += d.y * d.x;
    cov[1][1] += d.y * d.y;
    cov[1][2] += d.y * d.z;
    cov[2][0] += d.z * d.x;
    cov[2][1] += d.z * d.y;
    cov[2][2] += d.z * d.z;
  }

  double evals[3] = {0.0, 0.0, 0.0};
  double evecs[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
  JacobiEigenSymmetric3x3(cov, evals, evecs);

  int min_idx = 0;
  int max_idx = 0;
  for (int i = 1; i < 3; ++i) {
    if (evals[i] < evals[min_idx]) {
      min_idx = i;
    }
    if (evals[i] > evals[max_idx]) {
      max_idx = i;
    }
  }
  double max_eval = evals[max_idx];
  if (max_eval <= 0.0) {
    return plane;
  }
  double ratio = evals[min_idx] / max_eval;
  if (ratio < eps_collinear) {
    return plane;
  }
  Vec3 n_hat{evecs[0][min_idx], evecs[1][min_idx], evecs[2][min_idx]};
  n_hat = Normalize(n_hat);
  if (Norm(n_hat) <= 0.0) {
    return plane;
  }

  Vec3 ref = (std::abs(n_hat.x) < 0.9) ? Vec3{1.0, 0.0, 0.0} : Vec3{0.0, 1.0, 0.0};
  Vec3 e1 = Normalize(Cross(ref, n_hat));
  Vec3 e2 = Cross(n_hat, e1);

  plane.c = c;
  plane.n_hat = n_hat;
  plane.e1 = e1;
  plane.e2 = e2;
  plane.valid = true;
  return plane;
}

bool SegmentPlaneIntersection(const Vec3 &a,
                              const Vec3 &b,
                              const Plane &plane,
                              double eps_plane,
                              Vec3 *out_point) {
  if (!plane.valid) {
    return false;
  }
  double d_a = Dot(Sub(a, plane.c), plane.n_hat);
  double d_b = Dot(Sub(b, plane.c), plane.n_hat);
  if (d_a * d_b > 0.0) {
    return false;
  }
  if (std::abs(d_a) < eps_plane || std::abs(d_b) < eps_plane) {
    return false;
  }
  double denom = d_a - d_b;
  if (std::abs(denom) <= 0.0) {
    return false;
  }
  double t = d_a / denom;
  if (t <= 0.0 || t >= 1.0) {
    return false;
  }
  if (out_point != nullptr) {
    Vec3 diff = Sub(b, a);
    *out_point = Add(a, Scale(diff, t));
  }
  return true;
}

}  // namespace rna
