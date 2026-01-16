#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "entanglement.h"

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

  py::class_<rna::Vec2>(m, "Vec2")
      .def(py::init<double, double>(), py::arg("x") = 0.0, py::arg("y") = 0.0)
      .def_readwrite("x", &rna::Vec2::x)
      .def_readwrite("y", &rna::Vec2::y);

  py::class_<rna::Vec3>(m, "Vec3")
      .def(py::init<double, double, double>(),
           py::arg("x") = 0.0,
           py::arg("y") = 0.0,
           py::arg("z") = 0.0)
      .def_readwrite("x", &rna::Vec3::x)
      .def_readwrite("y", &rna::Vec3::y)
      .def_readwrite("z", &rna::Vec3::z);

  py::class_<rna::ResidueCoord>(m, "ResidueCoord")
      .def(py::init<int, std::vector<rna::Vec3>>(),
           py::arg("res_index"),
           py::arg("atoms"))
      .def_readwrite("res_index", &rna::ResidueCoord::res_index)
      .def_readwrite("atoms", &rna::ResidueCoord::atoms);

  py::enum_<rna::BasePair::Type>(m, "BasePairType")
      .value("UNCLASSIFIED", rna::BasePair::Type::kUnclassified)
      .value("CANONICAL", rna::BasePair::Type::kCanonical)
      .value("NON_CANONICAL", rna::BasePair::Type::kNonCanonical);

  py::enum_<rna::LoopKind>(m, "LoopKind")
      .value("HAIRPIN", rna::LoopKind::kHairpin)
      .value("INTERNAL", rna::LoopKind::kInternal)
      .value("MULTI", rna::LoopKind::kMulti)
      .value("UNKNOWN", rna::LoopKind::kUnknown);

  py::class_<rna::BasePair>(m, "BasePair")
      .def(py::init<int, int, rna::BasePair::Type>(),
           py::arg("i"),
           py::arg("j"),
           py::arg("bp_type") = rna::BasePair::Type::kUnclassified)
      .def_readwrite("i", &rna::BasePair::i)
      .def_readwrite("j", &rna::BasePair::j)
      .def_readwrite("bp_type", &rna::BasePair::bp_type);

  py::class_<rna::Loop>(m, "Loop")
      .def(py::init<>())
      .def_readwrite("id", &rna::Loop::id)
      .def_readwrite("kind", &rna::Loop::kind)
      .def_readwrite("closing_pairs", &rna::Loop::closing_pairs)
      .def_readwrite("boundary_residues", &rna::Loop::boundary_residues);

  py::class_<rna::Plane>(m, "Plane")
      .def(py::init<>())
      .def_readwrite("c", &rna::Plane::c)
      .def_readwrite("n_hat", &rna::Plane::n_hat)
      .def_readwrite("e1", &rna::Plane::e1)
      .def_readwrite("e2", &rna::Plane::e2)
      .def_readwrite("valid", &rna::Plane::valid);

  py::class_<rna::Polygon2D>(m, "Polygon2D")
      .def(py::init<>())
      .def_readwrite("vertices", &rna::Polygon2D::vertices)
      .def_readwrite("valid", &rna::Polygon2D::valid);

  py::class_<rna::Surface>(m, "Surface")
      .def(py::init<>())
      .def_readwrite("loop_id", &rna::Surface::loop_id)
      .def_readwrite("kind", &rna::Surface::kind)
      .def_readwrite("plane", &rna::Surface::plane)
      .def_readwrite("polygon", &rna::Surface::polygon)
      .def_readwrite("skip_residues", &rna::Surface::skip_residues);

  py::class_<rna::HitInfo>(m, "HitInfo")
      .def(py::init<>())
      .def_readwrite("loop_id", &rna::HitInfo::loop_id)
      .def_readwrite("segment_id", &rna::HitInfo::segment_id)
      .def_readwrite("point", &rna::HitInfo::point);

  py::class_<rna::Result>(m, "Result")
      .def(py::init<>())
      .def_readwrite("K", &rna::Result::K)
      .def_readwrite("hits", &rna::Result::hits);

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

  m.def(
      "build_loops",
      [](const std::vector<std::pair<int, int>> &bp_list,
         int n_res,
         bool include_multi,
         bool main_layer_only) {
        auto pairs = ToBasePairs(bp_list, rna::BasePair::Type::kCanonical);
        rna::LoopBuildOptions options;
        options.include_multi = include_multi;
        options.main_layer_only = main_layer_only;
        return rna::BuildLoops(pairs, n_res, options);
      },
      py::arg("bp_list"),
      py::arg("n_res"),
      py::arg("include_multi") = false,
      py::arg("main_layer_only") = false,
      "Build loops from base pairs.");

  m.def(
      "build_surfaces",
      [](const std::vector<rna::ResidueCoord> &coords,
         const std::vector<rna::Loop> &loops,
         int atom_index,
         double eps_collinear) {
        rna::SurfaceBuildOptions options;
        options.atom_index = atom_index;
        options.eps_collinear = eps_collinear;
        return rna::BuildSurfaces(coords, loops, options);
      },
      py::arg("coords"),
      py::arg("loops"),
      py::arg("atom_index") = 0,
      py::arg("eps_collinear") = 1e-6,
      "Build surfaces from loops and coordinates.");

  m.def(
      "evaluate_entanglement",
      [](const std::vector<rna::ResidueCoord> &coords,
         const std::vector<rna::Surface> &surfaces,
         int atom_index,
         double eps_plane,
         double eps_polygon) {
        rna::EvaluateOptions options;
        options.atom_index = atom_index;
        options.eps_plane = eps_plane;
        options.eps_polygon = eps_polygon;
        return rna::EvaluateEntanglement(coords, surfaces, options);
      },
      py::arg("coords"),
      py::arg("surfaces"),
      py::arg("atom_index") = 0,
      py::arg("eps_plane") = 1e-2,
      py::arg("eps_polygon") = 1e-2,
      "Evaluate entanglement for surfaces.");
}
