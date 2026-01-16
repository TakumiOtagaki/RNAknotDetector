from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, List, Optional, Sequence, Set, Tuple

from Bio.PDB import MMCIFParser, PDBParser


@dataclass
class ResidueCoord:
    res_index: int
    atoms: List[Tuple[float, float, float]]
    pdb_res_id: str


@dataclass
class BasePair:
    i: int
    j: int
    bp_type: Optional[str] = None


def load_coords(
    path: str,
    atom_names: Sequence[str] = ("C4'",),
    chain_id: Optional[str] = None,
    model_index: int = 0,
    missing_policy: str = "skip",
    include_hetero: bool = False,
) -> List[ResidueCoord]:
    parser = _select_parser(path)
    structure = parser.get_structure("rna", path)
    return load_coords_from_structure(
        structure,
        atom_names=atom_names,
        chain_id=chain_id,
        model_index=model_index,
        missing_policy=missing_policy,
        include_hetero=include_hetero,
    )


def load_coords_from_structure(
    structure,
    atom_names: Sequence[str] = ("C4'",),
    chain_id: Optional[str] = None,
    model_index: int = 0,
    missing_policy: str = "skip",
    include_hetero: bool = False,
) -> List[ResidueCoord]:
    model = structure[model_index]
    chain = _select_chain(model, chain_id=chain_id)
    coords: List[ResidueCoord] = []
    seq_index = 1

    for residue in chain.get_residues():
        if _should_skip_residue(residue, include_hetero=include_hetero):
            continue
        atom_coords = _extract_atoms(
            residue, atom_names=atom_names, missing_policy=missing_policy
        )
        if atom_coords is None:
            continue
        pdb_res_id = _format_residue_id(chain.id, residue.get_id())
        coords.append(ResidueCoord(seq_index, atom_coords, pdb_res_id))
        seq_index += 1

    return coords


def parse_dot_bracket(dot_bracket: str) -> List[BasePair]:
    stack: List[int] = []
    pairs: List[BasePair] = []
    for idx, char in enumerate(dot_bracket, start=1):
        if char == "(":
            stack.append(idx)
        elif char == ")":
            if not stack:
                raise ValueError("Unbalanced dot-bracket: too many closing parens")
            i = stack.pop()
            pairs.append(BasePair(i=i, j=idx))
        elif char == ".":
            continue
        else:
            raise ValueError(f"Unsupported dot-bracket symbol: {char}")
    if stack:
        raise ValueError("Unbalanced dot-bracket: too many opening parens")
    return pairs


def parse_bpseq(lines: Iterable[str]) -> List[BasePair]:
    pairs: List[BasePair] = []
    for raw in lines:
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 3:
            raise ValueError(f"Invalid bpseq line: {line}")
        i = int(parts[0])
        j = int(parts[2])
        if j > i:
            pairs.append(BasePair(i=i, j=j))
    return pairs


def _select_parser(path: str):
    lower = path.lower()
    if lower.endswith(".cif") or lower.endswith(".mmcif"):
        return MMCIFParser(QUIET=True)
    return PDBParser(QUIET=True)


def _select_chain(model, chain_id: Optional[str]):
    chains = list(model.get_chains())
    if chain_id is None:
        if len(chains) != 1:
            raise ValueError("chain_id is required when multiple chains are present")
        return chains[0]
    for chain in chains:
        if chain.id == chain_id:
            return chain
    raise ValueError(f"chain_id not found: {chain_id}")


def _should_skip_residue(residue, include_hetero: bool) -> bool:
    hetflag, _, _ = residue.get_id()
    if not include_hetero and hetflag != " ":
        return True
    resname = residue.get_resname().upper()
    if resname in {"HOH", "WAT"}:
        return True
    if not _is_rna_residue(resname):
        return True
    return False


def _extract_atoms(residue, atom_names: Sequence[str], missing_policy: str):
    coords: List[Tuple[float, float, float]] = []
    for name in atom_names:
        if name not in residue:
            if missing_policy == "skip":
                return None
            if missing_policy == "nan":
                coords.append((float("nan"), float("nan"), float("nan")))
                continue
            raise ValueError(f"Missing atom {name} in residue {residue.get_id()}")
        atom = residue[name]
        x, y, z = atom.get_coord().tolist()
        coords.append((float(x), float(y), float(z)))
    return coords


def _format_residue_id(chain_id: str, res_id) -> str:
    hetflag, resseq, icode = res_id
    icode_str = icode.strip()
    if icode_str:
        return f"{chain_id}:{resseq}{icode_str}"
    return f"{chain_id}:{resseq}"


def _is_rna_residue(resname: str) -> bool:
    return resname in _RNA_RESIDUES


_RNA_RESIDUES: Set[str] = {
    "A",
    "C",
    "G",
    "U",
    "I",
    "URA",
    "URA3",
    "ADE",
    "ADE3",
    "CYT",
    "CYT3",
    "GUA",
    "GUA3",
    "H2U",
    "PSU",
    "OMG",
    "2MG",
    "7MG",
    "M2G",
    "M5C",
    "1MA",
    "5MC",
    "OMC",
    "A23",
    "C23",
    "G23",
    "U23",
}
