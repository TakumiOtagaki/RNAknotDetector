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

enum class AtomKind {
  kSingle,
  kP,
  kC4,
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
  int res_a = 0;
  int res_b = 0;
  AtomKind atom_a = AtomKind::kSingle;
  AtomKind atom_b = AtomKind::kSingle;
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
  int multi_chunk = 12;       // chunk size factor for multiloop boundary splitting
};

struct EvaluateOptions {
  int atom_index = 0;
  int atom_index_p = 0;
  int atom_index_c4 = 1;
  enum class PolylineMode {
    kSingleAtom,
    kPC4Alternating,
  };
  PolylineMode polyline_mode = PolylineMode::kSingleAtom;
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
