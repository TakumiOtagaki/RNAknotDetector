from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from pymol import cmd
from pymol.cgo import BEGIN, COLOR, CYLINDER, END, SPHERE, TRIANGLES, VERTEX

import rnaknotdetector_core as core

from input_layer import parse_bpseq
from secstruct2bpseq import parse_secstruct, read_secstruct_file


@dataclass(frozen=True)
class ResidueKey:
    chain: str
    resi: str
    ins_code: str
    segi: str


def rnaknot_hit(
    model_name: str,
    ss_path: str,
    chain: str = "A",
    surface_mode: int = 1,
    polyline_mode: int = 1,
    main_layer_only: bool = True,
    eps_plane: float = 1e-2,
    eps_polygon: float = 1e-2,
    eps_collinear: float = 1e-6,
    cmd_obj=None,
) -> List[int]:
    cmd_handle = cmd_obj or cmd
    if cmd_handle is None:
        raise RuntimeError("pymol.cmd is not available")

    coords_cpp, res_id_map, atom_map = _extract_coords_from_pymol(
        cmd_handle, model_name, chain
    )
    if not coords_cpp:
        raise RuntimeError("No residues found for the requested model/chain")

    bp_list = _load_base_pairs(ss_path)
    if main_layer_only or not bp_list:
        bp_list = core.get_main_layer_pairs(bp_list)

    loops = core.build_loops(bp_list, len(res_id_map), main_layer_only=False)
    surfaces = core.build_surfaces(
        coords_cpp,
        loops,
        eps_collinear=eps_collinear,
        surface_mode=surface_mode,
    )
    result = core.evaluate_entanglement(
        coords_cpp,
        surfaces,
        polyline_mode=polyline_mode,
        eps_plane=eps_plane,
        eps_polygon=eps_polygon,
    )

    loop_map = {loop.id: loop for loop in loops}
    surface_map = {surface.loop_id: surface for surface in surfaces}
    hit_indices: List[int] = []

    for idx, hit in enumerate(result.hits, start=1):
        loop = loop_map.get(hit.loop_id)
        if loop is None:
            continue
        surface = surface_map.get(hit.loop_id)
        if surface is None:
            continue
        hit_indices.append(idx)
        loop_type = _format_loop_type(loop)
        closing_pairs = [(bp.i, bp.j) for bp in loop.closing_pairs]
        a_label = _format_endpoint(res_id_map, hit.res_a, hit.atom_a)
        b_label = _format_endpoint(res_id_map, hit.res_b, hit.atom_b)
        print(
            f"hit loop={hit.loop_id} type={loop_type} "
            f"pairs={closing_pairs} segment=({a_label},{b_label})"
        )
        selection_name = f"hit_{idx}"
        selection = _build_hit_selection(loop, hit, res_id_map, chain)
        cmd_handle.select(selection_name, f"({model_name} and {selection})")
        tri_count = len(surface.triangles) if surface.triangles is not None else 0
        print(f"[debug] hit={idx} loop={hit.loop_id} tri_n={tri_count}")
        _draw_hit_objects(
            cmd_handle,
            model_name,
            idx,
            hit,
            surface,
            atom_map,
        )

    return hit_indices


def _load_base_pairs(ss_path: str) -> List[Tuple[int, int]]:
    path = Path(ss_path)
    suffix = path.suffix.lower()
    if suffix == ".secstruct":
        _, secstruct = read_secstruct_file(str(path))
        pair_map = parse_secstruct(secstruct)
        return [(i, j) for i, j in enumerate(pair_map) if i > 0 and j > i]
    if suffix == ".bpseq":
        lines = path.read_text(encoding="utf-8").splitlines()
        pairs = parse_bpseq(lines)
        return [(bp.i, bp.j) for bp in pairs]
    raise ValueError(f"Unsupported ss_path extension: {path.suffix}")


def _extract_coords_from_pymol(
    cmd_handle,
    model_name: str,
    chain: str,
) -> Tuple[List[core.ResidueCoord], Dict[int, str], Dict[int, Dict[str, Tuple[float, float, float]]]]:
    model = cmd_handle.get_model(f"{model_name} and chain {chain}")
    residue_order: List[ResidueKey] = []
    residue_atoms: Dict[ResidueKey, Dict[str, Tuple[float, float, float]]] = {}

    for atom in model.atom:
        key = ResidueKey(atom.chain, atom.resi, atom.ins_code, atom.segi)
        if key not in residue_atoms:
            residue_atoms[key] = {}
            residue_order.append(key)
        residue_atoms[key][atom.name] = (float(atom.coord[0]), float(atom.coord[1]), float(atom.coord[2]))

    coords_cpp: List[core.ResidueCoord] = []
    res_id_map: Dict[int, str] = {}
    atom_map: Dict[int, Dict[str, Tuple[float, float, float]]] = {}
    for idx, key in enumerate(residue_order, start=1):
        atoms = residue_atoms.get(key, {})
        res_id_map[idx] = _format_residue_id(key)
        atom_map[idx] = atoms
        p = atoms.get("P")
        c4 = atoms.get("C4'")
        coords = [
            core.Vec3(*p) if p else core.Vec3(math.nan, math.nan, math.nan),
            core.Vec3(*c4) if c4 else core.Vec3(math.nan, math.nan, math.nan),
        ]
        coords_cpp.append(core.ResidueCoord(idx, coords))
    return coords_cpp, res_id_map, atom_map


def _format_residue_id(key: ResidueKey) -> str:
    ins = key.ins_code.strip() if key.ins_code else ""
    return f"{key.chain}:{key.resi}{ins}"


def _format_loop_type(loop) -> str:
    if loop.kind == core.LoopKind.INTERNAL and not loop.boundary_residues:
        return "LoopKind.STACKING"
    return str(loop.kind)


def _atom_kind_label(kind) -> str:
    if kind == core.AtomKind.P:
        return "P"
    if kind == core.AtomKind.C4:
        return "C4'"
    return "X"


def _format_endpoint(res_id_map: Dict[int, str], res_index: int, atom_kind) -> str:
    label = res_id_map.get(res_index, str(res_index))
    atom_label = _atom_kind_label(atom_kind)
    return f"{label}:{atom_label}"


def _build_hit_selection(
    loop,
    hit,
    res_id_map: Dict[int, str],
    chain: str,
) -> str:
    residues = set(loop.boundary_residues)
    for bp in loop.closing_pairs:
        residues.add(bp.i)
        residues.add(bp.j)
    residues.add(hit.res_a)
    residues.add(hit.res_b)
    resi = "+".join(str(idx) for idx in sorted(residues))
    return f"(chain {chain} and resi {resi})"


def _draw_hit_objects(
    cmd_handle,
    model_name: str,
    hit_index: int,
    hit,
    surface,
    atom_map: Dict[int, Dict[str, Tuple[float, float, float]]],
) -> None:
    segment = _segment_coords(hit, atom_map)
    if segment is None:
        return
    seg_a, seg_b = segment
    cgo = []
    cgo.extend([COLOR, 1.0, 0.0, 0.0])
    cgo.extend([SPHERE, hit.point.x, hit.point.y, hit.point.z, 0.4])
    cgo.extend([COLOR, 0.2, 0.6, 0.2])
    tri = _find_hit_triangle(seg_a, seg_b, surface.triangles)
    if tri is not None:
        cgo.extend([BEGIN, TRIANGLES])
        cgo.extend([VERTEX, tri[0][0], tri[0][1], tri[0][2]])
        cgo.extend([VERTEX, tri[1][0], tri[1][1], tri[1][2]])
        cgo.extend([VERTEX, tri[2][0], tri[2][1], tri[2][2]])
        cgo.extend([END])
    else:
        _append_polygon_fan(cgo, surface)
    cgo.extend([COLOR, 0.1, 0.3, 0.8])
    cgo.extend(
        [
            CYLINDER,
            seg_a[0],
            seg_a[1],
            seg_a[2],
            seg_b[0],
            seg_b[1],
            seg_b[2],
            0.2,
            0.1,
            0.3,
            0.8,
            0.1,
            0.3,
            0.8,
        ]
    )
    cmd_handle.load_cgo(cgo, f"{model_name}_hit_{hit_index}_geom")


def _segment_coords(
    hit,
    atom_map: Dict[int, Dict[str, Tuple[float, float, float]]],
) -> Optional[Tuple[Tuple[float, float, float], Tuple[float, float, float]]]:
    a = _atom_coord(hit.res_a, hit.atom_a, atom_map)
    b = _atom_coord(hit.res_b, hit.atom_b, atom_map)
    if a is None or b is None:
        return None
    return a, b


def _atom_coord(
    res_index: int,
    atom_kind,
    atom_map: Dict[int, Dict[str, Tuple[float, float, float]]],
) -> Optional[Tuple[float, float, float]]:
    atoms = atom_map.get(res_index, {})
    if atom_kind == core.AtomKind.P:
        return atoms.get("P")
    if atom_kind == core.AtomKind.C4:
        return atoms.get("C4'")
    return None


def _find_hit_triangle(
    a: Tuple[float, float, float],
    b: Tuple[float, float, float],
    triangles: Iterable[core.Triangle],
) -> Optional[Tuple[Tuple[float, float, float], Tuple[float, float, float], Tuple[float, float, float]]]:
    for tri in triangles:
        if _segment_intersects_triangle(a, b, tri):
            return (
                (tri.a.x, tri.a.y, tri.a.z),
                (tri.b.x, tri.b.y, tri.b.z),
                (tri.c.x, tri.c.y, tri.c.z),
            )
    return None


def _append_polygon_fan(cgo, surface) -> None:
    if not surface.plane.valid or not surface.polygon.valid:
        return
    vertices = surface.polygon.vertices
    if vertices is None or len(vertices) < 3:
        return
    center_x = sum(v.x for v in vertices) / len(vertices)
    center_y = sum(v.y for v in vertices) / len(vertices)
    center = _plane_point(surface.plane, center_x, center_y)
    cgo.extend([BEGIN, TRIANGLES])
    for i in range(len(vertices)):
        a2 = vertices[i]
        b2 = vertices[(i + 1) % len(vertices)]
        a3 = _plane_point(surface.plane, a2.x, a2.y)
        b3 = _plane_point(surface.plane, b2.x, b2.y)
        cgo.extend([VERTEX, center[0], center[1], center[2]])
        cgo.extend([VERTEX, a3[0], a3[1], a3[2]])
        cgo.extend([VERTEX, b3[0], b3[1], b3[2]])
    cgo.extend([END])


def _plane_point(plane, x: float, y: float) -> Tuple[float, float, float]:
    return (
        plane.c.x + plane.e1.x * x + plane.e2.x * y,
        plane.c.y + plane.e1.y * x + plane.e2.y * y,
        plane.c.z + plane.e1.z * x + plane.e2.z * y,
    )


def _segment_intersects_triangle(
    a: Tuple[float, float, float],
    b: Tuple[float, float, float],
    tri: core.Triangle,
    eps: float = 1e-8,
) -> bool:
    ax, ay, az = a
    bx, by, bz = b
    dirx = bx - ax
    diry = by - ay
    dirz = bz - az
    e1x = tri.b.x - tri.a.x
    e1y = tri.b.y - tri.a.y
    e1z = tri.b.z - tri.a.z
    e2x = tri.c.x - tri.a.x
    e2y = tri.c.y - tri.a.y
    e2z = tri.c.z - tri.a.z
    pvecx = diry * e2z - dirz * e2y
    pvecy = dirz * e2x - dirx * e2z
    pvecz = dirx * e2y - diry * e2x
    det = e1x * pvecx + e1y * pvecy + e1z * pvecz
    if abs(det) < eps:
        return False
    inv_det = 1.0 / det
    tvecx = ax - tri.a.x
    tvecy = ay - tri.a.y
    tvecz = az - tri.a.z
    u = (tvecx * pvecx + tvecy * pvecy + tvecz * pvecz) * inv_det
    if u < 0.0 or u > 1.0:
        return False
    qvecx = tvecy * e1z - tvecz * e1y
    qvecy = tvecz * e1x - tvecx * e1z
    qvecz = tvecx * e1y - tvecy * e1x
    v = (dirx * qvecx + diry * qvecy + dirz * qvecz) * inv_det
    if v < 0.0 or u + v > 1.0:
        return False
    t = (e2x * qvecx + e2y * qvecy + e2z * qvecz) * inv_det
    return 0.0 < t < 1.0


cmd.extend("rnaknot_hit", rnaknot_hit)
