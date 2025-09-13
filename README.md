# PointGroup

## Overview

**PointGroup** is a Python package for analyzing point group and space group characters (irreducible representations) from tight-binding models (Wannier90/pythtb/respack format). It provides utilities for symmetry operation analysis, character calculation, and group symbol conversion, leveraging libraries such as spglib and seekpath.

---

## Main Features

- **Tight-binding Model Wrapper**  
  - Wraps Wannier90/pythtb/respack models, manages lattice/orbital/band/k-path information (`TBModel` class).
- **Automatic Point/Space Group Detection**  
  - Uses spglib to detect international and Schonflies symbols and all symmetry operations (rotation/translation).
- **Character Calculation**  
  - Computes the character (trace) for a given k-point, band, and symmetry operation (`GetCharacter` class).
- **Symmetry Operation Classification**  
  - Classifies symmetry operations (rotation, mirror, inversion, etc.) and converts to Schonflies/international symbols.
- **Visualization Support**  
  - Utilities for structure visualization (e.g., CIF export via ASE).
- **Extended Character Tables**  
  - Built-in character tables for major point groups (`pg_ch_extra`).

---

## Directory Structure

```
pointgroup/
    __init__.py
    get_character.py      # Main class for character calculation
    tbmodel.py            # Tight-binding model wrapper
    sym_op.py             # Symmetry operation utilities
    schonflies.py         # Schonflies symbol conversion
    utils.py              # Logging/visualization utilities
    pglib/                # Point group helper modules
        judge_pg_sym_op.py
        output_manager.py
        pg_ch_extra.py    # Character tables
        transform.py
requirements.txt
README.md
```

---

## Example Usage

```python
from pointgroup.tbmodel import TBModel
from pointgroup.get_character import GetCharacter

tb = TBModel(wannier_path="...", wannier_prefix="...")
tb.gen_pythtb()
gc = GetCharacter(tb)
char = gc.get_character(kidx=0, orb_idx=0, op_idx=0)
print("Character:", char)
```

---

## Dependencies
- numpy
- spglib
- seekpath
- ase
- pythtb, respack (optional)

---

## License
Under Consideration

---

## Notes
- This package is intended for research in quantum chemistry and solid-state physics.
- For details, see docstrings and test codes in each module.

---

