#pragma once

#include <algorithm>
#include <limits>
#include <stdexcept>
#include <vector>

#include "entanglement.h"
#include "pair_utils.h"

namespace rna {

inline std::vector<int> BuildPairMap(const std::vector<BasePair> &base_pairs, int n_res) {
  std::vector<int> pair_map(n_res + 1, 0);
  for (const auto &bp : base_pairs) {
    if (bp.i <= 0 || bp.j <= 0 || bp.i > n_res || bp.j > n_res) {
      throw std::invalid_argument("Base pair index out of range");
    }
    if (bp.i == bp.j) {
      throw std::invalid_argument("Base pair cannot be self-paired");
    }
    auto [i, j] = SortedPair(bp);
    if (pair_map[i] != 0 || pair_map[j] != 0) {
      throw std::invalid_argument("Residue paired multiple times");
    }
    pair_map[i] = j;
    pair_map[j] = i;
  }
  return pair_map;
}

inline bool IsPaired(const std::vector<int> &pair_map, int idx) {
  return pair_map[idx] != 0;
}

// CollectUnpaired は、pair_map を使って 指定された区間 [start, end] にある「未ペアの残基インデックス」だけを集める関数です。
inline std::vector<int> CollectUnpaired(const std::vector<int> &pair_map, int start, int end) {
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
inline std::vector<BasePair> FindChildPairs(const std::vector<int> &pair_map, int i, int j) {
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
inline LoopKind ClassifyLoop(const std::vector<int> &pair_map,
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
    auto [k, l] = SortedPair(child_pairs[0]);
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

inline std::vector<int> BuildSkipResidues(const Loop &loop) {
  std::vector<int> skip;
  if (loop.closing_pairs.empty()) {
    return skip;
  }
  if (loop.kind == LoopKind::kHairpin) {
    auto [i, j] = SortedPair(loop.closing_pairs[0]);
    for (int k = i; k <= j; ++k) {
      skip.push_back(k);
    }
    return skip;
  }
  if (loop.kind == LoopKind::kInternal) {
    if (loop.closing_pairs.size() < 2) {
      auto [i, j] = SortedPair(loop.closing_pairs[0]);
      for (int k = i; k <= j; ++k) {
        skip.push_back(k);
      }
      return skip;
    }
    auto [i, j] = SortedPair(loop.closing_pairs[0]);
    auto [k, l] = SortedPair(loop.closing_pairs[1]);
    for (int idx = i; idx <= k; ++idx) {
      skip.push_back(idx);
    }
    for (int idx = l; idx <= j; ++idx) {
      skip.push_back(idx);
    }
    return skip;
  }
  if (loop.kind == LoopKind::kMulti) {
    int min_res = std::numeric_limits<int>::max();
    int max_res = std::numeric_limits<int>::min();
    for (const auto &pair : loop.closing_pairs) {
      auto [i, j] = SortedPair(pair);
      min_res = std::min(min_res, i);
      max_res = std::max(max_res, j);
      skip.push_back(i);
      skip.push_back(j);
    }
    if (min_res <= max_res) {
      for (int idx = min_res; idx <= max_res; ++idx) {
        skip.push_back(idx);
      }
    }
    return skip;
  }
  return skip;
}

}  // namespace rna
