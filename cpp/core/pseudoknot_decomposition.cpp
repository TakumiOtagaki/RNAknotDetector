#include "pseudoknot_decomposition.h"

#include <algorithm>
#include <cstdint>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

namespace rna {
namespace {

using Pair = std::pair<int, int>;

int64_t PairKey(int i, int j) {
  int a = std::min(i, j);
  int b = std::max(i, j);
  return (static_cast<int64_t>(a) << 32) | static_cast<uint32_t>(b);
}

std::vector<int> UniqueSortedResidues(const std::vector<Pair> &pairs) {
  std::vector<int> residues;
  residues.reserve(pairs.size() * 2);
  for (const auto &pair : pairs) {
    residues.push_back(pair.first);
    residues.push_back(pair.second);
  }
  std::sort(residues.begin(), residues.end());
  residues.erase(std::unique(residues.begin(), residues.end()), residues.end());
  return residues;
}

std::vector<Pair> CompressPairs(const std::vector<Pair> &pairs,
                                std::vector<int> *inv_hash) {
  *inv_hash = UniqueSortedResidues(pairs);
  std::unordered_map<int, int> hash;
  hash.reserve(inv_hash->size());
  for (size_t idx = 0; idx < inv_hash->size(); ++idx) {
    hash.emplace((*inv_hash)[idx], static_cast<int>(idx));
  }
  std::vector<Pair> compressed;
  compressed.reserve(pairs.size());
  for (const auto &pair : pairs) {
    compressed.emplace_back(hash.at(pair.first), hash.at(pair.second));
  }
  return compressed;
}

std::vector<Pair> DecompressPairs(const std::vector<Pair> &pairs,
                                  const std::vector<int> &inv_hash) {
  std::vector<Pair> decompressed;
  decompressed.reserve(pairs.size());
  for (const auto &pair : pairs) {
    decompressed.emplace_back(inv_hash.at(pair.first), inv_hash.at(pair.second));
  }
  return decompressed;
}

std::vector<Pair> ExtractMainLayerPairs(const std::vector<Pair> &pairs) {
  if (pairs.empty()) {
    return {};
  }

  std::vector<int> inv_hash;
  std::vector<Pair> compressed = CompressPairs(pairs, &inv_hash);
  int L = static_cast<int>(inv_hash.size());

  std::vector<std::vector<int>> gamma(L, std::vector<int>(L, 0));
  std::unordered_set<int64_t> pair_set;
  pair_set.reserve(compressed.size());
  for (const auto &pair : compressed) {
    pair_set.insert(PairKey(pair.first, pair.second));
  }

  auto gamma_value = [&gamma](int i, int j) {
    if (i < 0 || j < 0 || i >= static_cast<int>(gamma.size()) ||
        j >= static_cast<int>(gamma.size())) {
      return 0;
    }
    if (i > j) {
      return 0;
    }
    return gamma[i][j];
  };

  for (int d = 1; d < L; ++d) {
    for (int i = 0; i + d < L; ++i) {
      int j = i + d;
      int best = std::max(gamma_value(i + 1, j), gamma_value(i, j - 1));
      int diag = gamma_value(i + 1, j - 1);
      if (pair_set.count(PairKey(i, j)) != 0) {
        best = std::max(best, diag + 1);
      } else {
        best = std::max(best, diag);
      }
      for (int k = i; k < j; ++k) {
        best = std::max(best, gamma_value(i, k) + gamma_value(k + 1, j));
      }
      gamma[i][j] = best;
    }
  }

  std::vector<Pair> layer;
  std::vector<char> route(L, '.');
  std::vector<Pair> trace_stack;
  trace_stack.emplace_back(0, L - 1);

  while (!trace_stack.empty()) {
    Pair node = trace_stack.back();
    trace_stack.pop_back();
    int i = node.first;
    int j = node.second;
    if (i >= j) {
      continue;
    }
    if (gamma_value(i + 1, j) == gamma[i][j]) {
      trace_stack.emplace_back(i + 1, j);
      continue;
    }
    if (gamma_value(i, j - 1) == gamma[i][j]) {
      trace_stack.emplace_back(i, j - 1);
      continue;
    }
    if (pair_set.count(PairKey(i, j)) != 0 &&
        gamma_value(i + 1, j - 1) + 1 == gamma[i][j]) {
      if (route[i] == '.' && route[j] == '.') {
        route[i] = '(';
        route[j] = ')';
        layer.emplace_back(i, j);
        trace_stack.emplace_back(i + 1, j - 1);
        continue;
      }
    }
    for (int k = i; k < j; ++k) {
      if (gamma_value(i, k) + gamma_value(k + 1, j) == gamma[i][j]) {
        trace_stack.emplace_back(k + 1, j);
        trace_stack.emplace_back(i, k);
        break;
      }
    }
  }

  return DecompressPairs(layer, inv_hash);
}

}  // namespace

std::vector<BasePair> ExtractMainLayerFromBasePairs(
    const std::vector<BasePair> &base_pairs) {
  if (base_pairs.empty()) {
    return {};
  }

  std::vector<Pair> pairs;
  pairs.reserve(base_pairs.size());
  std::unordered_map<int64_t, BasePair::Type> type_map;
  type_map.reserve(base_pairs.size());

  for (const auto &bp : base_pairs) {
    if (bp.i == bp.j) {
      throw std::invalid_argument("Base pair cannot be self-paired");
    }
    int i = std::min(bp.i, bp.j);
    int j = std::max(bp.i, bp.j);
    pairs.emplace_back(i, j);
    type_map.emplace(PairKey(i, j), bp.bp_type);
  }

  std::vector<Pair> layer_pairs = ExtractMainLayerPairs(pairs);
  std::vector<BasePair> result;
  result.reserve(layer_pairs.size());
  for (const auto &pair : layer_pairs) {
    BasePair::Type bp_type = BasePair::Type::kUnclassified;
    auto it = type_map.find(PairKey(pair.first, pair.second));
    if (it != type_map.end()) {
      bp_type = it->second;
    }
    result.push_back(BasePair{pair.first, pair.second, bp_type});
  }
  return result;
}

}  // namespace rna
