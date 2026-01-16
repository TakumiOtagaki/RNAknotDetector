# RNAknotDetector

## Notes

### Loop boundary
In the current C++ core, `boundary` means the set of unpaired residues on the loop boundary.
It is used as a minimal placeholder for surface generation.

- hairpin: unpaired residues in (i+1 .. j-1)
- internal/bulge: unpaired residues between outer (i,j) and inner (k,l) pairs
- multi: unpaired residues in (i+1 .. j-1) (coarse placeholder)

### PyMOL debug
`python/pymol_debug.py` contains helper utilities for interactive inspection.
It expects a pybind11 module named `rnaknotdetector_core` to be built.

```sh
/opt/homebrew/Cellar/pymol/3.1.0_3/libexec/bin/python -m pip install pandas numpy pybind11 biopython
PYTHON=/opt/homebrew/Cellar/pymol/3.1.0_3/libexec/bin/python make clean
PYTHON=/opt/homebrew/Cellar/pymol/3.1.0_3/libexec/bin/python make
```
わざわざこの 方法でやらないといけない。
ここクソ面倒な現象の原因は uv で pymol を install できないという点にある。
ちょっとやってみようとしたけど無理だった。

```pymol console
print_main_layer_pairs 6t3r, A, /Users/ootagakitakumi/RNAknotDetector/examples/6t3r.secstruct
color_main_layer_pairs 6t3r, A, /Users/ootagakitakumi/RNAknotDetector/examples/6t3r.secstruct, pink
```


```pymol console
load /Users/ootagakitakumi/RNAknotDetector/examples/Example02.pdb
print_main_layer_pairs Example02, A, /Users/ootagakitakumi/RNAknotDetector/examples/Example02.secstruct
color_main_layer_pairs Example02, A, /Users/ootagakitakumi/RNAknotDetector/examples/Example02.secstruct, pink
# pkv Example02
```

pseudoknot visualizer の出力と一致することを 6t3r と Example02 で確認しておいた
