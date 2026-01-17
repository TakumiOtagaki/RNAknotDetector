#pragma once

#include <vector>

#include "geometry2d.h"
#include "geometry3d.h"

namespace rna {

struct ResidueCoord {
  int res_index = 0;  // 1-based
  std::vector<Vec3> atoms;
};

struct BasePair {
  int i = 0;  // 1-based
  int j = 0;  // 1-based
  enum class Type {
    kUnclassified,
    kCanonical,
    kNonCanonical,
  };
  Type bp_type = Type::kUnclassified;
};

enum class LoopKind {
  kHairpin,
  kInternal,
  kMulti,
  kUnknown,
};

struct Loop {
  int id = 0;
  LoopKind kind = LoopKind::kUnknown;
  std::vector<BasePair> closing_pairs;
  std::vector<int> boundary_residues;  // unpaired residues on loop boundary
};

struct LoopBuildOptions {
  bool main_layer_only = false;
};

enum class SurfaceMode {
  kBestFitPlane,
  kTrianglePlanes,
};

struct Surface {
  int loop_id = 0;
  LoopKind kind = LoopKind::kUnknown;
  std::vector<BasePair> closing_pairs;
  Plane plane;
  Polygon2D polygon;
  std::vector<Triangle> triangles;
  std::vector<int> skip_residues;
};

struct HitInfo {
  int loop_id = 0;
  int segment_id = 0;  // backbone segment index i for (i, i+1)
  Vec3 point;
};

struct Result {
  int K = 0;
  std::vector<HitInfo> hits;
};

struct SurfaceBuildOptions {
  int atom_index = 0;         // which atom in ResidueCoord::atoms to use
  double eps_collinear = 1e-6;  // ratio threshold for near-collinearity
  SurfaceMode surface_mode = SurfaceMode::kTrianglePlanes;
};

struct EvaluateOptions {
  int atom_index = 0;
  double eps_plane = 1e-2;
  double eps_polygon = 1e-2;
  double eps_triangle = 1e-8;
};

std::vector<Loop> BuildLoops(const std::vector<BasePair> &base_pairs,
                             int n_res,
                             const LoopBuildOptions &options = {});

std::vector<BasePair> ExtractMainLayer(const std::vector<BasePair> &base_pairs);

std::vector<Surface> BuildSurfaces(const std::vector<ResidueCoord> &coords,
                                   const std::vector<Loop> &loops,
                                   const SurfaceBuildOptions &options = {});

Result EvaluateEntanglement(const std::vector<ResidueCoord> &coords,
                            const std::vector<Surface> &surfaces,
                            const EvaluateOptions &options = {});

}  // namespace rna
