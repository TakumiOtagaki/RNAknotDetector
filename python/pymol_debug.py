from __future__ import annotations

from typing import Iterable, List, Optional, Sequence, Tuple

try:
    from pymol import cmd
except ImportError:  # pragma: no cover - used only inside PyMOL
    cmd = None

try:
    import rnaknotdetector_core as core
except ImportError:  # pragma: no cover - bindings not built yet
    core = None


def color_multiloop_pairs(
    bp_list: Sequence[Tuple[int, int]],
    n_res: int,
    chain_id: Optional[str] = None,
    color: str = "red",
    cmd_obj=None,
) -> List[int]:
    if core is None:
        raise RuntimeError("rnaknotdetector_core module is not available")
    cmd_handle = cmd_obj or cmd
    if cmd_handle is None:
        raise RuntimeError("pymol.cmd is not available")

    pairs = core.get_multiloop_pairs(bp_list, n_res)
    residues = _unique_residues_from_pairs(pairs)
    if not residues:
        return []

    selection = _selection_from_residues(residues, chain_id=chain_id)
    cmd_handle.color(color, selection)
    return residues


def _unique_residues_from_pairs(pairs: Iterable[Tuple[int, int]]) -> List[int]:
    residues = {idx for pair in pairs for idx in pair}
    return sorted(residues)


def _selection_from_residues(residues: Sequence[int], chain_id: Optional[str]) -> str:
    resi = "+".join(str(idx) for idx in residues)
    if chain_id:
        return f"(chain {chain_id} and resi {resi})"
    return f"resi {resi}"
