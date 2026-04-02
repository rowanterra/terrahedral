# Terrahedral

CIF -> coordination geometry analysis for metalloenzymes and minerals.

## Install

```bash
pip install -e .
```

Or with optional dependencies:

```bash
pip install -e ".[full]"   # adds gemmi
pip install -e ".[dev]"    # adds pytest, ruff
```

## Quick start

```python
from terrahedral import CoordinationSite

# From a local CIF
site = CoordinationSite.from_cif("structure.cif", metal="Cu1")

# From the PDB (fetches automatically)
site = CoordinationSite.from_pdb("6GTL", metal_asym="B")

# From COD
site = CoordinationSite.from_cod("4070506", metal="Pt9")

# Analyze
site.classify_geometry()
print(site.summary())
print(f"tau4 = {site.tau4()}")
```

## CLI

```bash
terrahedral info structure.cif --metal Cu1
terrahedral info --pdb 6GTL --metal Cu --asym B
terrahedral list-metals --pdb 6GTL
terrahedral fetch --pdb 6GTL
```

## Tests

```bash
pip install -e ".[dev]"
pytest
```

## License

Business Source License 1.1 — free for academic and non-commercial use.
Commercial use requires a paid license. Contact rowan.r.terra@gmail.com.

On April 2, 2030 (or 4 years after each version's release), the code
converts to GPL-3.0.
