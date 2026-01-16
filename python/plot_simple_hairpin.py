from __future__ import annotations

import argparse
from typing import List, Tuple

import matplotlib.pyplot as plt


def load_c4_coords(path: str) -> List[Tuple[int, float, float, float]]:
    coords: List[Tuple[int, float, float, float]] = []
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith("ATOM"):
                continue
            atom_name = line[12:16].strip()
            if atom_name != "C4'":
                continue
            res_id = int(line[22:26])
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            coords.append((res_id, x, y, z))
    coords.sort(key=lambda t: t[0])
    return coords


def main() -> int:
    parser = argparse.ArgumentParser(description="Plot C4' coordinates in 3D.")
    parser.add_argument("pdb_path", help="Input PDB path.")
    args = parser.parse_args()

    coords = load_c4_coords(args.pdb_path)
    if not coords:
        raise SystemExit("No C4' atoms found.")

    xs = [c[1] for c in coords]
    ys = [c[2] for c in coords]
    zs = [c[3] for c in coords]

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    ax.plot(xs, ys, zs, marker="o")
    for res_id, x, y, z in coords:
        ax.text(x, y, z, str(res_id), size=8)

    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.set_title("C4' Trace")
    plt.tight_layout()
    plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
