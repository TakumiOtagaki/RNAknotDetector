#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "entanglement_core.h"

namespace py = pybind11;

namespace {

std::vector<rna::BasePair> ToBasePairs(
    const std::vector<std::pair<int, int>> &bp_list,
    rna::BasePair::Type bp_type) {
  std::vector<rna::BasePair> result;
  result.reserve(bp_list.size());
  for (const auto &pair : bp_list) {
    result.push_back(rna::BasePair{pair.first, pair.second, bp_type});
  }
  return result;
}

}  // namespace

PYBIND11_MODULE(rnaknotdetector_core, m) {
  m.doc() = "Minimal bindings for RNAknotDetector core.";

  py::enum_<rna::BasePair::Type>(m, "BasePairType")
      .value("UNCLASSIFIED", rna::BasePair::Type::kUnclassified)
      .value("CANONICAL", rna::BasePair::Type::kCanonical)
      .value("NON_CANONICAL", rna::BasePair::Type::kNonCanonical);

  m.def(
      "get_multiloop_pairs",
      [](const std::vector<std::pair<int, int>> &bp_list, int n_res) {
        auto pairs = ToBasePairs(bp_list, rna::BasePair::Type::kCanonical);
        auto multi_pairs = rna::CollectMultiLoopPairs(pairs, n_res);
        std::vector<std::pair<int, int>> result;
        result.reserve(multi_pairs.size());
        for (const auto &pair : multi_pairs) {
          result.emplace_back(pair.i, pair.j);
        }
        return result;
      },
      py::arg("bp_list"),
      py::arg("n_res"),
      "Return closing pairs that belong to multi-loops.");

  m.def(
      "get_main_layer_pairs",
      [](const std::vector<std::pair<int, int>> &bp_list) {
        auto pairs = ToBasePairs(bp_list, rna::BasePair::Type::kCanonical);
        auto main_pairs = rna::ExtractMainLayer(pairs);
        std::vector<std::pair<int, int>> result;
        result.reserve(main_pairs.size());
        for (const auto &pair : main_pairs) {
          result.emplace_back(pair.i, pair.j);
        }
        return result;
      },
      py::arg("bp_list"),
      "Return base pairs in the main (maximum) pseudoknot-free layer.");
}
