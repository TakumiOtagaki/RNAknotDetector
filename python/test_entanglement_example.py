from __future__ import annotations

import argparse
from typing import Iterable, List, Tuple

import rnaknotdetector_core as core

from input_layer import ResidueCoord as PyResidueCoord
from input_layer import load_coords
from secstruct2bpseq import parse_secstruct, read_secstruct_file


def _to_cpp_coords(py_coords: Iterable[PyResidueCoord]) -> List[core.ResidueCoord]:
    coords: List[core.ResidueCoord] = []
    for res in py_coords:
        atoms = [core.Vec3(*atom) for atom in res.atoms]
        coords.append(core.ResidueCoord(res.res_index, atoms))
    return coords


def _pair_map_to_list(pair_map: List[int]) -> List[Tuple[int, int]]:
    pairs: List[Tuple[int, int]] = []
    for i in range(1, len(pair_map)):
        j = pair_map[i]
        if j > i:
            pairs.append((i, j))
    return pairs


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Evaluate entanglement on an example PDB + secstruct."
    )
    parser.add_argument("pdb_path", help="Input PDB/mmCIF path.")
    parser.add_argument("secstruct_path", help="Secstruct file path.")
    parser.add_argument("--chain", default=None, help="Chain ID (if multiple chains).")
    parser.add_argument("--eps-plane", type=float, default=1e-2)
    parser.add_argument("--eps-polygon", type=float, default=1e-2)
    parser.add_argument("--eps-collinear", type=float, default=1e-6)
    parser.add_argument(
        "--include-multi",
        action="store_true",
        help="Include multiloop surfaces when building loops.",
    )
    parser.add_argument(
        "--main-layer-only",
        action="store_true",
        help="Build loops from main layer pairs only (default).",
    )
    args = parser.parse_args()

    coords_py = load_coords(args.pdb_path, chain_id=args.chain)
    coords_cpp = _to_cpp_coords(coords_py)

    _, secstruct = read_secstruct_file(args.secstruct_path)
    pair_map = parse_secstruct(secstruct)
    bp_list = _pair_map_to_list(pair_map)
    print(f"[debug] input bp count = {len(bp_list)}")

    if args.main_layer_only or not bp_list:
        main_pairs = core.get_main_layer_pairs(bp_list)
        print(f"[debug] main layer bp count = {len(main_pairs)}")
        loop_pairs = main_pairs
    else:
        loop_pairs = bp_list
    loops = core.build_loops(
        loop_pairs,
        len(pair_map) - 1,
        include_multi=args.include_multi,
        main_layer_only=False,
    )
    print(f"[debug] loops built = {len(loops)}")
    surfaces = core.build_surfaces(coords_cpp, loops, eps_collinear=args.eps_collinear)
    valid_surfaces = sum(1 for s in surfaces if s.plane.valid and s.polygon.valid)
    print(f"[debug] surfaces built = {len(surfaces)} valid = {valid_surfaces}")
    result = core.evaluate_entanglement(
        coords_cpp,
        surfaces,
        eps_plane=args.eps_plane,
        eps_polygon=args.eps_polygon,
    )

    print(f"K = {result.K}")
    for hit in result.hits:
        p = hit.point
        print(
            f"hit loop={hit.loop_id} segment={hit.segment_id} "
            f"point=({p.x:.3f},{p.y:.3f},{p.z:.3f})"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
