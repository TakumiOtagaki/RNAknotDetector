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


def _unique_residues_from_pairs(pairs: Iterable[Tuple[int, int]]) -> List[int]:
    residues = {idx for pair in pairs for idx in pair}
    return sorted(residues)


def _selection_from_residues(residues: Sequence[int], chain_id: Optional[str]) -> str:
    resi = "+".join(str(idx) for idx in residues)
    if chain_id:
        return f"(chain {chain_id} and resi {resi})"
    return f"resi {resi}"
