# RNAknotDetector

## Notes

### Loop boundary
In the current C++ core, `boundary` means the set of unpaired residues on the loop boundary.
It is used as a minimal placeholder for surface generation.

- hairpin: unpaired residues in (i+1 .. j-1)
- internal/bulge: unpaired residues between outer (i,j) and inner (k,l) pairs
- multi: unpaired residues in (i+1 .. j-1) (coarse placeholder)
