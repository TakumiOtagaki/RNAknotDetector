#include "entanglement_core.h"
#include "pseudoknot_decomposition.h"

#include <algorithm>
#include <stdexcept>

namespace rna {
namespace {

std::vector<int> BuildPairMap(const std::vector<BasePair> &base_pairs, int n_res) {
  std::vector<int> pair_map(n_res + 1, 0);
  for (const auto &bp : base_pairs) {
    if (bp.i <= 0 || bp.j <= 0 || bp.i > n_res || bp.j > n_res) {
      throw std::invalid_argument("Base pair index out of range");
    }
    if (bp.i == bp.j) {
      throw std::invalid_argument("Base pair cannot be self-paired");
    }
    int i = std::min(bp.i, bp.j);
    int j = std::max(bp.i, bp.j);
    if (pair_map[i] != 0 || pair_map[j] != 0) {
      throw std::invalid_argument("Residue paired multiple times");
    }
    pair_map[i] = j;
    pair_map[j] = i;
  }
  return pair_map;
}

bool IsPaired(const std::vector<int> &pair_map, int idx) {
  return pair_map[idx] != 0;
}


std::vector<int> CollectUnpaired(const std::vector<int> &pair_map, int start, int end) {
  std::vector<int> residues;
  for (int k = start; k <= end; ++k) {
    if (!IsPaired(pair_map, k)) {
      residues.push_back(k);
    }
  }
  return residues;
}

// Find immediate child base pairs inside (i, j).
// A child pair is the first paired region encountered when scanning the
// interval; nested pairs inside that region are ignored.
// Used to count how many stems close the loop bounded by (i, j).
std::vector<BasePair> FindChildPairs(const std::vector<int> &pair_map, int i, int j) {
  std::vector<BasePair> child_pairs;
  int depth = 0;
  for (int idx = i + 1; idx <= j - 1; ++idx) {
    if (!IsPaired(pair_map, idx)) {
      continue;
    }
    int partner = pair_map[idx];
    if (idx < partner) {
      if (depth == 0) {
        child_pairs.push_back(BasePair{idx, partner, BasePair::Type::kUnclassified});
      }
      depth++;
    } else if (idx > partner) {
      depth--;
    }
  }
  return child_pairs;
}

// Classify loop by counting immediate child pairs within (i, j).
// closing_pairs includes the outer pair (i, j) and each immediate child pair.
// boundary collects unpaired residues that lie on the loop boundary.
// Rules: 0 child -> hairpin, 1 child -> internal/bulge/stacking, 2+ -> multi.
LoopKind ClassifyLoop(const std::vector<int> &pair_map,
                      int i,
                      int j,
                      std::vector<int> *boundary,
                      std::vector<BasePair> *closing_pairs) {
  closing_pairs->clear();
  closing_pairs->push_back(BasePair{i, j, BasePair::Type::kUnclassified});

  std::vector<BasePair> child_pairs = FindChildPairs(pair_map, i, j);
  closing_pairs->insert(closing_pairs->end(), child_pairs.begin(), child_pairs.end());

  if (child_pairs.empty()) {
    *boundary = CollectUnpaired(pair_map, i + 1, j - 1);
    return LoopKind::kHairpin;
  }
  if (child_pairs.size() == 1) {
    int k = std::min(child_pairs[0].i, child_pairs[0].j);
    int l = std::max(child_pairs[0].i, child_pairs[0].j);
    auto left = CollectUnpaired(pair_map, i + 1, k - 1);
    auto right = CollectUnpaired(pair_map, l + 1, j - 1);
    boundary->clear();
    boundary->reserve(left.size() + right.size());
    boundary->insert(boundary->end(), left.begin(), left.end());
    boundary->insert(boundary->end(), right.begin(), right.end());
    return LoopKind::kInternal;
  }

  *boundary = CollectUnpaired(pair_map, i + 1, j - 1);
  return LoopKind::kMulti;
}

}  // namespace


// Build closed elements (loops) from base pairs.
// Each loop corresponds to an outer closing pair (i, j) and immediate child pairs.
// boundary_residues holds unpaired residues on the loop boundary (minimal set).
// Pseudoknots are not supported.
std::vector<Loop> BuildLoops(const std::vector<BasePair> &base_pairs,
                             int n_res,
                             const LoopBuildOptions &options) {
  if (n_res <= 0) {
    throw std::invalid_argument("n_res must be positive");
  }
  // Placeholder: assumes pseudoknot-free input; no validation for crossing pairs.
  std::vector<BasePair> filtered_pairs = base_pairs;
  if (options.main_layer_only) {
    filtered_pairs = ExtractMainLayer(base_pairs);
  }
  std::vector<int> pair_map = BuildPairMap(filtered_pairs, n_res);

  std::vector<Loop> loops;
  int loop_id = 1;
  for (int i = 1; i <= n_res; ++i) {
    int j = pair_map[i];
    if (j == 0 || i > j) {
      continue;
    }
    std::vector<int> boundary;
    std::vector<BasePair> closing_pairs;
    LoopKind kind = ClassifyLoop(pair_map, i, j, &boundary, &closing_pairs);
    if (kind == LoopKind::kMulti && !options.include_multi) {
      continue;
    }
    Loop loop;
    loop.id = loop_id++;
    loop.kind = kind;
    loop.closing_pairs = std::move(closing_pairs);
    loop.boundary_residues = std::move(boundary);
    loops.push_back(std::move(loop));
  }
  return loops;
}

std::vector<BasePair> ExtractMainLayer(const std::vector<BasePair> &base_pairs) {
  if (base_pairs.empty()) {
    return {};
  }
  return ExtractMainLayerFromBasePairs(base_pairs);
}

}  // namespace rna
