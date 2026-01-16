#pragma once

#include <string>
#include <vector>

namespace rna {

struct Vec3 {
  double x;
  double y;
  double z;
};

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
  bool include_multi = false;
};

// Build closed elements (loops) from base pairs.
// Each loop corresponds to an outer closing pair (i, j) and immediate child pairs.
// boundary_residues holds unpaired residues on the loop boundary (minimal set).
// Placeholder: no surface generation and no validation for crossing pairs.
std::vector<Loop> BuildLoops(const std::vector<BasePair> &base_pairs,
                             int n_res,
                             const LoopBuildOptions &options = {});

// Collect closing pairs that belong to multi-loops.
// Placeholder: used for PyMOL debug coloring.
std::vector<BasePair> CollectMultiLoopPairs(const std::vector<BasePair> &base_pairs,
                                            int n_res,
                                            const LoopBuildOptions &options = {});

}  // namespace rna
