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
        "--surface-mode",
        type=int,
        default=1,
        help="Surface mode: 0=best-fit plane, 1=triangle planes.",
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
    loops = core.build_loops(loop_pairs, len(pair_map) - 1, main_layer_only=False)
    print(f"[debug] loops built = {len(loops)}")
    surfaces = core.build_surfaces(
        coords_cpp,
        loops,
        eps_collinear=args.eps_collinear,
        surface_mode=args.surface_mode,
    )
    valid_surfaces = sum(1 for s in surfaces if s.plane.valid and s.polygon.valid)
    print(f"[debug] surfaces built = {len(surfaces)} valid = {valid_surfaces}")
    result = core.evaluate_entanglement(
        coords_cpp,
        surfaces,
        eps_plane=args.eps_plane,
        eps_polygon=args.eps_polygon,
    )

    print("[debug] indices are 1-based (residue indices and segment ids)")
    print(f"K = {result.K}")
    loop_map = {loop.id: loop for loop in loops}
    res_id_map = {res.res_index: res.pdb_res_id for res in coords_py}
    for hit in result.hits:
        loop = loop_map.get(hit.loop_id)
        if loop is None:
            continue
        p = hit.point
        loop_type = str(loop.kind)
        closing_pairs = [(bp.i, bp.j) for bp in loop.closing_pairs]
        i = hit.segment_id
        j = hit.segment_id + 1
        i_label = res_id_map.get(i, str(i))
        j_label = res_id_map.get(j, str(j))
        print(
            f"hit loop={hit.loop_id} type={loop_type} "
            f"pairs={closing_pairs} segment=({i_label},{j_label}) "
            f"point=({p.x:.3f},{p.y:.3f},{p.z:.3f})"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
