#include "entanglement_core.h"

#include <algorithm>
#include <stdexcept>
#include <unordered_set>

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

bool HasPairedInRange(const std::vector<int> &pair_map, int start, int end) {
  for (int k = start; k <= end; ++k) {
    if (IsPaired(pair_map, k)) {
      return true;
    }
  }
  return false;
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

LoopKind ClassifyLoop(const std::vector<int> &pair_map, int i, int j,
                      std::vector<int> *boundary) {
  if (i + 1 > j - 1) {
    boundary->clear();
    return LoopKind::kHairpin;
  }
  bool has_paired_inside = HasPairedInRange(pair_map, i + 1, j - 1);
  if (!has_paired_inside) {
    *boundary = CollectUnpaired(pair_map, i + 1, j - 1);
    return LoopKind::kHairpin;
  }

  int k = 0;
  for (int idx = i + 1; idx <= j - 1; ++idx) {
    if (IsPaired(pair_map, idx)) {
      k = idx;
      break;
    }
  }
  if (k == 0) {
    *boundary = CollectUnpaired(pair_map, i + 1, j - 1);
    return LoopKind::kHairpin;
  }

  int l = pair_map[k];
  if (k > l) {
    std::swap(k, l);
  }
  if (k <= i || l >= j) {
    *boundary = CollectUnpaired(pair_map, i + 1, j - 1);
    return LoopKind::kUnknown;
  }

  bool left_branch = HasPairedInRange(pair_map, i + 1, k - 1);
  bool right_branch = HasPairedInRange(pair_map, l + 1, j - 1);
  if (left_branch || right_branch) {
    *boundary = CollectUnpaired(pair_map, i + 1, j - 1);
    return LoopKind::kMulti;
  }

  boundary->clear();
  auto left = CollectUnpaired(pair_map, i + 1, k - 1);
  auto right = CollectUnpaired(pair_map, l + 1, j - 1);
  boundary->reserve(left.size() + right.size());
  boundary->insert(boundary->end(), left.begin(), left.end());
  boundary->insert(boundary->end(), right.begin(), right.end());
  return LoopKind::kInternal;
}

}  // namespace

std::vector<Loop> BuildLoops(const std::vector<BasePair> &base_pairs,
                             int n_res,
                             const LoopBuildOptions &options) {
  if (n_res <= 0) {
    throw std::invalid_argument("n_res must be positive");
  }
  std::vector<int> pair_map = BuildPairMap(base_pairs, n_res);

  std::vector<Loop> loops;
  int loop_id = 1;
  for (int i = 1; i <= n_res; ++i) {
    int j = pair_map[i];
    if (j == 0 || i > j) {
      continue;
    }
    std::vector<int> boundary;
    LoopKind kind = ClassifyLoop(pair_map, i, j, &boundary);
    if (kind == LoopKind::kMulti && !options.include_multi) {
      continue;
    }
    Loop loop;
    loop.id = loop_id++;
    loop.i = i;
    loop.j = j;
    loop.kind = kind;
    loop.boundary_residues = std::move(boundary);
    loops.push_back(std::move(loop));
  }
  return loops;
}

}  // namespace rna
