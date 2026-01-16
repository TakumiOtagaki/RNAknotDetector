from __future__ import annotations

import argparse
import sys
from typing import Dict, List, Sequence, Tuple


UNPAIRED_CHARS = {".", "-", "x", "X"}
OPEN_TO_CLOSE = {
    "(": ")",
    "[": "]",
    "{": "}",
    "<": ">",
}
CLOSE_TO_OPEN = {v: k for k, v in OPEN_TO_CLOSE.items()}
SEQUENCE_CHARS = set("ACGUTNacgutn")
SECSTRUCT_CHARS = set(OPEN_TO_CLOSE.keys()) | set(CLOSE_TO_OPEN.keys()) | UNPAIRED_CHARS


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Convert a dot-bracket-style secstruct file to BPSEQ."
    )
    parser.add_argument("secstruct_path", help="Input secstruct file path.")
    parser.add_argument(
        "-o",
        "--output",
        help="Output BPSEQ path (default: stdout).",
        default="-",
    )
    args = parser.parse_args()

    sequence, secstruct = read_secstruct_file(args.secstruct_path)
    pair_map = parse_secstruct(secstruct)
    lines = format_bpseq(sequence, pair_map)

    if args.output == "-":
        sys.stdout.write("".join(lines))
    else:
        with open(args.output, "w", encoding="utf-8") as handle:
            handle.write("".join(lines))
    return 0


def read_secstruct_file(path: str) -> Tuple[str, str]:
    raw_lines = _load_lines(path)
    seq_parts: List[str] = []
    ss_parts: List[str] = []

    for line in raw_lines:
        if _looks_like_secstruct(line):
            ss_parts.append(line)
            continue
        if _looks_like_sequence(line):
            seq_parts.append(line)
            continue
        raise ValueError(f"Unrecognized line in secstruct file: {line}")

    if not ss_parts:
        raise ValueError("No secstruct string found in input.")

    secstruct = "".join(ss_parts)
    sequence = "".join(seq_parts) if seq_parts else "N" * len(secstruct)

    if len(sequence) != len(secstruct):
        raise ValueError(
            f"Sequence length ({len(sequence)}) does not match secstruct length "
            f"({len(secstruct)})."
        )

    return sequence, secstruct


def _load_lines(path: str) -> List[str]:
    lines: List[str] = []
    with open(path, "r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            if line.startswith(">"):
                continue
            lines.append(line)
    return lines


def _looks_like_sequence(line: str) -> bool:
    return all(char in SEQUENCE_CHARS for char in line)


def _looks_like_secstruct(line: str) -> bool:
    return all(char in SECSTRUCT_CHARS for char in line)


def parse_secstruct(secstruct: str) -> List[int]:
    pair_map: List[int] = [0] * (len(secstruct) + 1)
    stacks: Dict[str, List[int]] = {ch: [] for ch in OPEN_TO_CLOSE}

    for idx, char in enumerate(secstruct, start=1):
        if char in OPEN_TO_CLOSE:
            stacks[char].append(idx)
            continue
        if char in CLOSE_TO_OPEN:
            open_char = CLOSE_TO_OPEN[char]
            if not stacks[open_char]:
                raise ValueError(f"Unbalanced secstruct: unexpected {char} at {idx}")
            i = stacks[open_char].pop()
            pair_map[i] = idx
            pair_map[idx] = i
            continue
        if char in UNPAIRED_CHARS:
            continue
        raise ValueError(f"Unsupported secstruct symbol: {char}")

    for open_char, stack in stacks.items():
        if stack:
            raise ValueError(
                f"Unbalanced secstruct: missing {OPEN_TO_CLOSE[open_char]}"
            )

    return pair_map


def format_bpseq(sequence: str, pair_map: Sequence[int]) -> List[str]:
    lines: List[str] = []
    for idx in range(1, len(pair_map)):
        base = sequence[idx - 1]
        partner = pair_map[idx]
        lines.append(f"{idx} {base} {partner}\n")
    return lines


if __name__ == "__main__":
    raise SystemExit(main())
