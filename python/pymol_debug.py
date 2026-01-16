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

try:
    from secstruct2bpseq import parse_secstruct, read_secstruct_file
except ImportError:  # pragma: no cover - used only when helper is available
    parse_secstruct = None
    read_secstruct_file = None


def print_main_layer_pairs(
    model_name: str,
    chain_id: str,
    secstruct_path: Optional[str] = None,
    cmd_obj=None,
) -> List[Tuple[int, int]]:
    if core is None:
        raise RuntimeError("rnaknotdetector_core module is not available")
    cmd_handle = cmd_obj or cmd
    if cmd_handle is None:
        raise RuntimeError("pymol.cmd is not available")
    if read_secstruct_file is None or parse_secstruct is None:
        raise RuntimeError("secstruct2bpseq helpers are not available")

    if secstruct_path is None:
        secstruct_path = f"{model_name}.secstruct"

    sequence, secstruct = read_secstruct_file(secstruct_path)
    pair_map = parse_secstruct(secstruct)
    bp_list = [(i, j) for i, j in enumerate(pair_map) if i > 0 and j > i]
    main_pairs = core.get_main_layer_pairs(bp_list)

    residues = _extract_residues(cmd_handle, model_name, chain_id)
    _print_pairs(main_pairs, residues, len(secstruct))
    return main_pairs


def _unique_residues_from_pairs(pairs: Iterable[Tuple[int, int]]) -> List[int]:
    residues = {idx for pair in pairs for idx in pair}
    return sorted(residues)


def _selection_from_residues(residues: Sequence[int], chain_id: Optional[str]) -> str:
    resi = "+".join(str(idx) for idx in residues)
    if chain_id:
        return f"(chain {chain_id} and resi {resi})"
    return f"resi {resi}"


def _extract_residues(cmd_handle, model_name: str, chain_id: str) -> List[str]:
    model = cmd_handle.get_model(f"{model_name} and chain {chain_id}")
    seen = set()
    residues: List[str] = []
    for atom in model.atom:
        key = (atom.chain, atom.resi, atom.ins_code, atom.segi)
        if key in seen:
            continue
        seen.add(key)
        residues.append(atom.resi)
    return residues


def _print_pairs(
    pairs: Sequence[Tuple[int, int]],
    residues: Sequence[str],
    expected_length: int,
) -> None:
    if residues and len(residues) != expected_length:
        print(
            f"[warning] residue count ({len(residues)}) does not match "
            f"secstruct length ({expected_length})"
        )
    if not pairs:
        print("[info] main layer pairs: none")
        return

    print("[info] main layer pairs:")
    for i, j in pairs:
        if residues and i <= len(residues) and j <= len(residues):
            print(f"  {i}({residues[i-1]}) - {j}({residues[j-1]})")
        else:
            print(f"  {i} - {j}")


if cmd is not None:  # pragma: no cover - PyMOL runtime only
    cmd.extend("print_main_layer_pairs", print_main_layer_pairs)
