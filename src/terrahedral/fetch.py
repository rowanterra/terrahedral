"""
Fetch CIF files from online databases.

Supported sources:
  - RCSB PDB  (mmCIF format)  → fetch_pdb("6GTL")
  - COD       (small-molecule) → fetch_cod("4070506")

Search:
  - search_pdb("lanmodulin")  → [{pdb_id, title, resolution, ...}, ...]
  - search_cod("quartz")      → [{cod_id, formula, ...}, ...]
"""

from __future__ import annotations

import json
import os
import re
import urllib.parse
import urllib.request
import urllib.error
from pathlib import Path
from typing import Optional

RCSB_URL = "https://files.rcsb.org/download/{pdb_id}.cif"
COD_URL = "https://www.crystallography.net/cod/{cod_id}.cif"

RCSB_SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
RCSB_DATA_URL = "https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
RCSB_GRAPHQL_URL = "https://data.rcsb.org/graphql"

COD_SEARCH_URL = "https://www.crystallography.net/cod/result"

_DEFAULT_CACHE = Path.home() / ".cache" / "terrahedral"


def _ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


# ─────────────────────────────────────────────────────────────────────
# PDB search
# ─────────────────────────────────────────────────────────────────────

# Known siderophore names — searching these without requiring an iron
# entity returns mostly biosynthetic enzymes / apo-binding proteins,
# not the actual iron–siderophore complex structures.
_SIDEROPHORE_TERMS = {
    "enterobactin", "ferric enterobactin",
    "desferrioxamine", "desferrioxamine b", "deferoxamine", "deferrioxamine",
    "rhizoferrin", "ferrichrome", "pyoverdine", "pyochelin",
    "vibriobactin", "bacillibactin", "petrobactin", "yersiniabactin",
    "aerobactin", "mycobactin", "staphyloferrin", "salmochelin",
    "achromobactin", "citrate siderophore",
    "ferrioxamine", "coprogen", "rhodotorulic acid",
}

# Metal comp_ids to require for siderophore searches (Fe²⁺ / Fe³⁺)
_IRON_COMP_IDS = ["FE", "FE2"]


def _is_siderophore_query(query: str) -> bool:
    """Return True if the query appears to be searching for a siderophore."""
    q = query.lower().strip()
    # Direct match or query contains a known siderophore term
    for term in _SIDEROPHORE_TERMS:
        if term in q or q in term:
            return True
    # Heuristic: query contains 'siderophore'
    if "siderophore" in q:
        return True
    return False


def _search_pdb_with_metal(
    query: str,
    metal_comp_ids: list[str],
    *,
    max_rows: int = 24,
    sort: list[dict] | None = None,
) -> list[str]:
    """Combined RCSB search: full-text AND require specific metal entities."""
    metal_nodes = [
        {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_nonpolymer_entity_container_identifiers.nonpolymer_comp_id",
                "operator": "exact_match",
                "value": mid,
            },
        }
        for mid in metal_comp_ids
    ]

    payload = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "full_text",
                    "parameters": {"value": query},
                },
                {
                    "type": "group",
                    "logical_operator": "or",
                    "nodes": metal_nodes,
                },
            ],
        },
        "return_type": "entry",
        "request_options": {
            "paginate": {"start": 0, "rows": max_rows},
            "sort": sort or [{"sort_by": "score", "direction": "desc"}],
            "scoring_strategy": "combined",
        },
    }

    try:
        req = urllib.request.Request(
            RCSB_SEARCH_URL,
            data=json.dumps(payload).encode("utf-8"),
            headers={"Content-Type": "application/json"},
            method="POST",
        )
        with urllib.request.urlopen(req, timeout=15) as resp:
            data = json.loads(resp.read())
    except (urllib.error.URLError, urllib.error.HTTPError, json.JSONDecodeError):
        return []

    return [r["identifier"] for r in data.get("result_set", [])]


def search_pdb(
    query: str,
    *,
    max_results: int = 8,
    sort_by: str = "resolution",
) -> list[dict]:
    """
    Search RCSB PDB by text query (protein name, enzyme, organism, etc.).

    Parameters
    ----------
    query : free-text search (e.g. "lanmodulin", "catalase", "iron oxygenase")
    max_results : max entries to return (default 8)
    sort_by : "resolution" (best first), "date" (newest first), or "score" (relevance)

    Returns
    -------
    List of dicts with: pdb_id, title, resolution, method, date, authors, metal_count
    Sorted by resolution (best first) by default.
    """
    if not query.strip():
        return []

    # RCSB Search API v2
    sort_map = {
        "resolution": [{"sort_by": "rcsb_entry_info.resolution_combined", "direction": "asc"}],
        "date": [{"sort_by": "rcsb_accession_info.deposit_date", "direction": "desc"}],
        "score": [{"sort_by": "score", "direction": "desc"}],
    }
    sort_spec = sort_map.get(sort_by, sort_map["resolution"])

    # ── Siderophore-aware search ──
    # Siderophore names return mostly biosynthetic enzymes / apo-receptors
    # unless we also require an iron entity in the structure.
    if _is_siderophore_query(query):
        iron_ids = _search_pdb_with_metal(
            query, _IRON_COMP_IDS,
            max_rows=max_results * 3,
            sort=sort_spec,
        )
        if iron_ids:
            return _fetch_pdb_metadata(iron_ids, max_results=max_results)
        # Fall through to plain search if combined query found nothing

    search_payload = {
        "query": {
            "type": "terminal",
            "service": "full_text",
            "parameters": {"value": query},
        },
        "return_type": "entry",
        "request_options": {
            "paginate": {"start": 0, "rows": max_results * 3},
            "sort": sort_spec,
            "scoring_strategy": "combined",
        },
    }

    try:
        req = urllib.request.Request(
            RCSB_SEARCH_URL,
            data=json.dumps(search_payload).encode("utf-8"),
            headers={"Content-Type": "application/json"},
            method="POST",
        )
        with urllib.request.urlopen(req, timeout=15) as resp:
            data = json.loads(resp.read())
    except (urllib.error.URLError, urllib.error.HTTPError, json.JSONDecodeError):
        return []

    result_set = data.get("result_set", [])
    if not result_set:
        return []

    pdb_ids = [r["identifier"] for r in result_set[:max_results * 3]]

    # Batch-fetch metadata via GraphQL (one call)
    return _fetch_pdb_metadata(pdb_ids, max_results=max_results)


def _fetch_pdb_metadata(
    pdb_ids: list[str],
    *,
    max_results: int = 8,
) -> list[dict]:
    """Fetch title, resolution, method, date, metal info for a batch of PDB IDs."""
    if not pdb_ids:
        return []

    ids_str = ", ".join(f'"{pid}"' for pid in pdb_ids)
    gql = f"""{{
      entries(entry_ids: [{ids_str}]) {{
        rcsb_id
        struct {{
          title
        }}
        rcsb_entry_info {{
          resolution_combined
          experimental_method
          deposited_atom_count
        }}
        rcsb_accession_info {{
          deposit_date
        }}
        audit_author {{
          name
        }}
        nonpolymer_entities {{
          nonpolymer_comp {{
            chem_comp {{
              id
              name
              formula
              type
            }}
          }}
        }}
      }}
    }}"""

    try:
        req = urllib.request.Request(
            RCSB_GRAPHQL_URL,
            data=json.dumps({"query": gql}).encode("utf-8"),
            headers={"Content-Type": "application/json"},
            method="POST",
        )
        with urllib.request.urlopen(req, timeout=15) as resp:
            result = json.loads(resp.read())
    except (urllib.error.URLError, urllib.error.HTTPError, json.JSONDecodeError):
        # Fallback: return bare IDs
        return [{"pdb_id": pid, "title": "", "resolution": None} for pid in pdb_ids[:max_results]]

    entries = result.get("data", {}).get("entries", [])
    out = []
    for entry in entries:
        if not entry:
            continue

        pdb_id = entry.get("rcsb_id", "")
        info = entry.get("rcsb_entry_info") or {}
        accession = entry.get("rcsb_accession_info") or {}
        struct = entry.get("struct") or {}
        authors_list = entry.get("audit_author") or []

        # Detect metals in nonpolymer entities
        metals_found = []
        metals_seen = set()
        from terrahedral.core import METAL_ELEMENTS
        for npe in (entry.get("nonpolymer_entities") or []):
            comp = (npe.get("nonpolymer_comp") or {}).get("chem_comp") or {}
            comp_id = comp.get("id", "")
            formula = comp.get("formula", "") or ""
            # Check if the component is a bare metal ion (e.g. FE, MN, CA)
            if comp_id in METAL_ELEMENTS:
                if comp_id not in metals_seen:
                    metals_found.append(comp_id)
                    metals_seen.add(comp_id)
            else:
                # Check formula for metals embedded in ligands
                # (e.g. siderophore–iron complexes: "C25 H45 Fe N6 O8")
                for token in formula.split():
                    el = re.match(r"([A-Z][a-z]?)", token)
                    if el and el.group(1) in METAL_ELEMENTS:
                        sym = el.group(1)
                        if sym not in metals_seen:
                            metals_found.append(sym)
                            metals_seen.add(sym)

        resolution = info.get("resolution_combined")
        if isinstance(resolution, list) and resolution:
            resolution = resolution[0]

        authors = [a.get("name", "") for a in authors_list[:3]]
        author_str = ", ".join(authors)
        if len(authors_list) > 3:
            author_str += f" (+{len(authors_list) - 3})"

        out.append({
            "pdb_id": pdb_id,
            "title": (struct.get("title") or "").strip(),
            "resolution": resolution,
            "method": info.get("experimental_method") or "",
            "date": (accession.get("deposit_date") or "")[:10],
            "authors": author_str,
            "metals": metals_found,
            "metal_count": len(metals_found),
        })

    # Sort: entries with metals first, then by resolution
    out.sort(key=lambda e: (
        0 if e["metal_count"] > 0 else 1,
        e["resolution"] if e["resolution"] else 999,
    ))

    return out[:max_results]


# ─────────────────────────────────────────────────────────────────────
# COD search
# ─────────────────────────────────────────────────────────────────────

# Common compound names → element sets for smart query parsing.
# Enables "iron nitrate" → search by elements Fe + N + O.
_COMPOUND_NAMES = {
    # Metals
    "iron": "Fe", "ferric": "Fe", "ferrous": "Fe",
    "copper": "Cu", "cupric": "Cu", "cuprous": "Cu",
    "zinc": "Zn", "nickel": "Ni", "cobalt": "Co",
    "manganese": "Mn", "chromium": "Cr", "vanadium": "V",
    "titanium": "Ti", "molybdenum": "Mo", "tungsten": "W",
    "platinum": "Pt", "palladium": "Pd", "gold": "Au",
    "silver": "Ag", "mercury": "Hg", "tin": "Sn",
    "lead": "Pb", "bismuth": "Bi", "aluminum": "Al",
    "aluminium": "Al", "magnesium": "Mg", "calcium": "Ca",
    "sodium": "Na", "potassium": "K", "barium": "Ba",
    "strontium": "Sr", "lithium": "Li", "ruthenium": "Ru",
    "rhodium": "Rh", "iridium": "Ir", "osmium": "Os",
    "cadmium": "Cd", "zirconium": "Zr",
    # Anions / ligands → elements they contribute
    "nitrate": ["N", "O"], "nitrite": ["N", "O"],
    "sulfate": ["S", "O"], "sulphate": ["S", "O"],
    "sulfide": ["S"], "sulphide": ["S"],
    "phosphate": ["P", "O"], "carbonate": ["C", "O"],
    "chloride": ["Cl"], "bromide": ["Br"], "iodide": ["I"],
    "fluoride": ["F"], "oxide": ["O"], "hydroxide": ["O", "H"],
    "cyanide": ["C", "N"], "thiocyanate": ["S", "C", "N"],
    "acetate": ["C", "O", "H"], "oxalate": ["C", "O"],
    "citrate": ["C", "O", "H"],
    "bipyridine": ["C", "N", "H"], "phenanthroline": ["C", "N", "H"],
    "hydrate": ["H", "O"], "aqua": ["H", "O"], "water": ["H", "O"],
}

# All valid element symbols for detection in queries
_ALL_ELEMENTS = {
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",
    "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr",
    "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    "In", "Sn", "Sb", "Te", "I", "Xe",
    "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb",
    "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
    "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn",
    "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am",
}


def _parse_query_elements(query: str) -> list[str]:
    """
    Extract element symbols from a natural-language compound query.

    Handles:
      - Direct element symbols: "Fe N O" → [Fe, N, O]
      - Compound names: "iron nitrate" → [Fe, N, O]
      - Formulas: "Fe(NO3)3" → [Fe, N, O]
      - Mixed: "copper sulfate" → [Cu, S, O]
    """
    elements = set()
    words = re.split(r'[\s,;()\[\]]+', query.strip().lower())

    for word in words:
        # Strip trailing numbers/charges: "fe3+" → "fe", "iron(iii)" → "iron"
        clean = re.sub(r'[\d()+\-iii]+$', '', word).strip()
        if not clean:
            continue

        # Check compound name lookup (also try stripping numeric prefixes
        # like "di", "tri", "tetra", "penta", "hexa", "mono", etc.)
        if clean in _COMPOUND_NAMES:
            val = _COMPOUND_NAMES[clean]
            if isinstance(val, list):
                elements.update(val)
            else:
                elements.add(val)
            continue

        # Strip common numeric prefixes and try again
        stripped = re.sub(
            r'^(mono|di|tri|tetra|penta|hexa|hepta|octa|nona|deca|bis|tris)',
            '', clean,
        )
        if stripped and stripped != clean and stripped in _COMPOUND_NAMES:
            val = _COMPOUND_NAMES[stripped]
            if isinstance(val, list):
                elements.update(val)
            else:
                elements.add(val)
            continue

        # Check if it's a raw element symbol (case-insensitive)
        # Try exact match first, then capitalized
        cap = clean.capitalize()
        if cap in _ALL_ELEMENTS:
            elements.add(cap)
            continue

    # Also try to parse formula fragments like "Fe(NO3)3", "CuSO4", "FeCl3"
    # Extract element symbols via regex: uppercase letter + optional lowercase
    # But skip Roman numerals (I, V, X used in oxidation states)
    formula_els = re.findall(r'([A-Z][a-z]?)', query)
    # Filter: don't count isolated "I" or "V" or "X" that appear inside
    # parenthesized Roman numerals like (III), (IV), (II)
    roman_context = set()
    for m in re.finditer(r'\(([IVX]+)\)', query):
        for ch in m.group(1):
            roman_context.add(ch)
    for el in formula_els:
        if el in _ALL_ELEMENTS:
            # Skip single-char elements that are part of Roman numerals
            if len(el) == 1 and el in roman_context:
                continue
            elements.add(el)

    return sorted(elements)


def _build_cod_params(
    query: str = "",
    *,
    elements: Optional[list[str]] = None,
    exclude_elements: Optional[list[str]] = None,
    formula: Optional[str] = None,
    year_min: Optional[int] = None,
    strict_min_elements: Optional[int] = None,
    strict_max_elements: Optional[int] = None,
) -> dict:
    """Build COD REST API query parameters."""
    params = {"format": "json"}

    if formula:
        params["formula"] = formula
    if query.strip():
        params["text"] = query.strip()
    if elements:
        for i, el in enumerate(elements[:8], 1):
            params[f"el{i}"] = el
    if exclude_elements:
        for i, el in enumerate(exclude_elements[:4], 1):
            params[f"nel{i}"] = el
    if year_min is not None:
        params["year"] = str(year_min)
    if strict_min_elements is not None:
        params["strictmin"] = str(strict_min_elements)
    if strict_max_elements is not None:
        params["strictmax"] = str(strict_max_elements)

    return params


def _cod_query(params: dict, max_results: int = 8) -> list[dict]:
    """Execute a COD REST API query and return parsed results."""
    try:
        url = f"{COD_SEARCH_URL}?{urllib.parse.urlencode(params)}"
        req = urllib.request.Request(url, headers={"User-Agent": "terrahedral/0.1"})
        with urllib.request.urlopen(req, timeout=15) as resp:
            raw = resp.read().decode("utf-8", errors="replace")
    except (urllib.error.URLError, urllib.error.HTTPError) as e:
        print(f"  [cod] Query failed: {e}")
        return []

    # Try JSON first (preferred), fall back to CSV, then HTML
    raw_stripped = raw.strip()
    if raw_stripped.startswith("{") or raw_stripped.startswith("["):
        return _parse_cod_json(raw, max_results=max_results)
    elif "," in raw.split("\n")[0] and not raw_stripped.startswith("<"):
        return _parse_cod_csv(raw, max_results=max_results)
    else:
        # HTML fallback — extract COD IDs from links
        cod_ids = re.findall(r'/cod/(\d+)\.cif', raw)
        cod_ids = list(dict.fromkeys(cod_ids))[:max_results * 2]
        if cod_ids:
            return _fetch_cod_metadata(cod_ids, max_results=max_results)
        return []


def search_cod(
    query: str,
    *,
    max_results: int = 8,
    elements: Optional[list[str]] = None,
    exclude_elements: Optional[list[str]] = None,
    formula: Optional[str] = None,
    year_min: Optional[int] = None,
    strict_min_elements: Optional[int] = None,
    strict_max_elements: Optional[int] = None,
) -> list[dict]:
    """
    Search the Crystallography Open Database.

    Supports natural-language queries like "iron nitrate", "copper sulfate",
    "FeCl3", or direct element lists. Automatically extracts elements from
    compound names and formulas.

    Parameters
    ----------
    query : search text — compound name, formula, mineral name, elements, etc.
    max_results : max entries to return
    elements : list of elements that MUST be present (e.g. ["Fe", "S"])
    exclude_elements : list of elements that must NOT be present
    formula : exact Hill-notation formula (e.g. "Fe N3 O9")
    year_min : only entries published on or after this year
    strict_min_elements : minimum distinct elements in formula
    strict_max_elements : maximum distinct elements in formula

    Returns
    -------
    List of dicts with: cod_id, formula, title, mineral, spacegroup,
    cell_volume, year, has_coordinates, source
    """
    if not query.strip() and not elements and not formula:
        return []

    results = []
    seen_ids: set[str] = set()

    def _dedup_add(new_results: list[dict]) -> None:
        for r in new_results:
            cid = r.get("cod_id", "")
            if cid and cid not in seen_ids:
                seen_ids.add(cid)
                results.append(r)

    # ── Strategy 1: If explicit formula given, search by formula ──
    if formula:
        params = _build_cod_params(
            formula=formula,
            elements=elements, exclude_elements=exclude_elements,
            year_min=year_min, strict_min_elements=strict_min_elements,
            strict_max_elements=strict_max_elements,
        )
        _dedup_add(_cod_query(params, max_results=max_results))
        if len(results) >= max_results:
            return results[:max_results]

    # ── Strategy 2: Parse elements from the query text ──
    # "iron nitrate" → [Fe, N, O]; "CuSO4" → [Cu, S, O]
    parsed_els = _parse_query_elements(query) if query.strip() else []
    search_elements = elements or parsed_els

    # Detect if query is "inorganic" — no carbon-containing words/elements.
    # If so, auto-exclude C to avoid drowning in organometallic hits.
    _ORGANIC_WORDS = {
        "bipyridine", "phenanthroline", "acetate", "oxalate", "citrate",
        "porphyrin", "salen", "edta", "amino", "methyl", "ethyl", "benzene",
        "phenyl", "organic", "carbonate", "cyanide", "thiocyanate",
    }
    query_lower = query.strip().lower()
    query_is_organic = (
        "C" in parsed_els
        or any(w in query_lower for w in _ORGANIC_WORDS)
    )
    auto_exclude = list(exclude_elements) if exclude_elements else []
    if not query_is_organic and parsed_els and "C" not in parsed_els:
        if "C" not in auto_exclude:
            auto_exclude.append("C")

    if search_elements and len(search_elements) >= 2:
        # Element-based search (most reliable for compounds)
        params = _build_cod_params(
            elements=search_elements, exclude_elements=auto_exclude or None,
            year_min=year_min, strict_min_elements=strict_min_elements,
            strict_max_elements=strict_max_elements,
        )
        # If we parsed exactly the right number of elements, constrain
        # to that exact count to avoid getting huge organic molecules
        # that happen to contain Fe+N+O.  But only if strictmin/max
        # weren't already set by the caller.
        if parsed_els and not strict_min_elements and not strict_max_elements:
            nel = len(search_elements)
            # Allow +1 for hydrogen (hydrates are common)
            params["strictmin"] = str(nel)
            params["strictmax"] = str(nel + 2)

        # Request extra results for ranking/filtering
        _dedup_add(_cod_query(params, max_results=max_results * 3))
        if len(results) >= max_results:
            results = _rank_cod_results(results, search_elements)
            return results[:max_results]

    # ── Strategy 3: Text search as supplement / fallback ──
    if query.strip():
        params = _build_cod_params(
            query=query, elements=elements, exclude_elements=auto_exclude or None,
            year_min=year_min, strict_min_elements=strict_min_elements,
            strict_max_elements=strict_max_elements,
        )
        _dedup_add(_cod_query(params, max_results=max_results * 2))

    results = _rank_cod_results(results, search_elements)
    return results[:max_results]


def _rank_cod_results(
    results: list[dict],
    target_elements: list[str],
) -> list[dict]:
    """
    Rank COD search results by relevance.

    Prefers:
      - Formulas with fewer total atoms (simple salts over organometallics)
      - Formulas with fewer distinct elements (closer match to query)
      - Newer structures (more recent year)
    """
    target_set = set(target_elements) if target_elements else set()

    def _score(r: dict) -> tuple:
        formula = r.get("formula", "")

        # Count distinct elements in formula
        formula_els = set(re.findall(r'([A-Z][a-z]?)', formula))
        n_elements = len(formula_els) if formula_els else 99

        # Count total atoms (sum of all numeric subscripts)
        total_atoms = 0
        for m in re.finditer(r'([A-Z][a-z]?)(\d*\.?\d*)', formula):
            n = m.group(2)
            total_atoms += float(n) if n else 1.0

        # Bonus: how well elements match the query (extra elements = penalty)
        extra_elements = len(formula_els - target_set) if target_set else 0

        # Year (prefer newer, 0 if unknown)
        try:
            year = int(r.get("year", "0"))
        except (ValueError, TypeError):
            year = 0

        # Sort key: fewer extra elements, fewer total atoms, newer year
        return (extra_elements, total_atoms, -year)

    results.sort(key=_score)
    return results


def _parse_cod_json(json_text: str, *, max_results: int = 8) -> list[dict]:
    """Parse COD JSON search results into structured dicts."""
    try:
        data = json.loads(json_text)
    except (json.JSONDecodeError, ValueError):
        return []

    results = []

    # COD JSON can be a dict keyed by COD ID, or a list
    if isinstance(data, dict):
        items = list(data.items())
    elif isinstance(data, list):
        items = [(str(i), d) for i, d in enumerate(data)]
    else:
        return []

    for cod_id, entry in items:
        if len(results) >= max_results:
            break

        if isinstance(entry, dict):
            rec = entry
            # The key itself is the COD ID in dict format
            cod_id_str = str(rec.get("file", cod_id)).strip()
        else:
            continue

        formula = str(rec.get("formula", "")).strip().strip("'\" ")
        mineral = str(rec.get("mineral", "")).strip().strip("'\" ")
        title = str(rec.get("compname", rec.get("chemname", ""))).strip().strip("'\" ")
        sg = str(rec.get("sg", rec.get("sgHall", ""))).strip().strip("'\" ")
        year = str(rec.get("year", "")).strip()
        vol = str(rec.get("vol", "")).strip()
        z = str(rec.get("Z", rec.get("Zprime", ""))).strip()

        has_coords = True
        try:
            if z and int(float(z)) == 0:
                has_coords = False
        except (ValueError, TypeError):
            pass

        if not title and mineral:
            title = mineral
        if not title and formula:
            title = formula

        results.append({
            "cod_id": cod_id_str,
            "formula": formula,
            "title": title,
            "mineral": mineral,
            "spacegroup": sg,
            "year": year,
            "cell_volume": vol,
            "has_coordinates": has_coords,
            "source": "COD",
        })

    return results


def _parse_cod_csv(csv_text: str, *, max_results: int = 8) -> list[dict]:
    """Parse COD CSV search results into structured dicts."""
    import csv
    import io

    reader = csv.DictReader(io.StringIO(csv_text))
    results = []

    for row in reader:
        if len(results) >= max_results:
            break

        cod_id = row.get("file", "").strip()
        if not cod_id:
            continue

        formula = row.get("formula", "").strip().strip("'\"")
        mineral = row.get("mineral", "").strip().strip("'\"")
        title = row.get("compname", row.get("title", "")).strip().strip("'\"")
        sg = row.get("sg", row.get("sgHall", "")).strip().strip("'\"")
        year = row.get("year", "").strip()
        vol = row.get("vol", "").strip()
        z = row.get("Z", row.get("Zprime", "")).strip()

        has_coords = True
        try:
            if z and int(float(z)) == 0:
                has_coords = False
        except (ValueError, TypeError):
            pass

        if not title and mineral:
            title = mineral

        results.append({
            "cod_id": cod_id,
            "formula": formula,
            "title": title,
            "mineral": mineral,
            "spacegroup": sg,
            "year": year,
            "cell_volume": vol,
            "has_coordinates": has_coords,
            "source": "COD",
        })

    return results


def search_amcsd(
    query: str,
    *,
    max_results: int = 8,
) -> list[dict]:
    """
    Search the American Mineralogist Crystal Structure Database.

    Uses a curated index of common minerals (instant, always works),
    then tries AMCSD directly, then falls back to COD 9xxxxxx filtering.
    """
    if not query.strip():
        return []

    q = query.strip().lower()
    results = []

    # Strategy 1: Curated mineral index — instant, always works
    # Includes simple fuzzy matching for common misspellings
    def _fuzzy_match(query_str: str, target_str: str, threshold: int = 2) -> bool:
        """Check if query is a substring match OR within edit distance threshold."""
        ql = query_str.lower()
        tl = target_str.lower()
        # Exact substring
        if ql in tl or tl in ql:
            return True
        # Edit distance for short queries (likely mineral names)
        if len(ql) >= 4 and len(tl) >= 4:
            # Simple Levenshtein for the mineral name comparison
            # Only compute if lengths are close enough
            if abs(len(ql) - len(tl)) > threshold:
                return False
            # DP edit distance
            m, n = len(ql), len(tl)
            dp = list(range(n + 1))
            for i in range(1, m + 1):
                prev = dp[0]
                dp[0] = i
                for j in range(1, n + 1):
                    tmp = dp[j]
                    if ql[i - 1] == tl[j - 1]:
                        dp[j] = prev
                    else:
                        dp[j] = 1 + min(prev, dp[j], dp[j - 1])
                    prev = tmp
            if dp[n] <= threshold:
                return True
        return False

    for entry in _MINERAL_INDEX:
        name_match = _fuzzy_match(q, entry["mineral"].lower())
        formula_match = q in entry.get("formula", "").lower()
        element_match = any(q.lower() == el.lower() for el in entry.get("elements", []))
        if name_match or formula_match or element_match:
            results.append({
                "cod_id": entry["cod_id"],
                "amcsd_id": entry.get("amcsd_id", ""),
                "formula": entry.get("formula", ""),
                "title": entry["mineral"],
                "mineral": entry["mineral"],
                "spacegroup": entry.get("spacegroup", ""),
                "year": entry.get("year", ""),
                "has_coordinates": True,
                "source": "AMCSD",
            })

    if len(results) >= max_results:
        return results[:max_results]

    # Strategy 2: Query AMCSD at rruff.geo.arizona.edu
    try:
        import ssl
        params = urllib.parse.urlencode({"mineral": query.strip()})
        url = f"https://rruff.geo.arizona.edu/AMS/result.php?{params}"
        req = urllib.request.Request(url, headers={"User-Agent": "terrahedral/0.1"})
        ctx = ssl.create_default_context()
        ctx.check_hostname = False
        ctx.verify_mode = ssl.CERT_NONE
        with urllib.request.urlopen(req, timeout=8, context=ctx) as resp:
            html = resp.read().decode("utf-8", errors="replace")

        # Parse AMCSD IDs from download links
        amcsd_ids = re.findall(r'download\.php\?id=(\d+)\.cif', html)
        amcsd_ids = list(dict.fromkeys(amcsd_ids))

        # Also try parsing database codes from the page
        db_codes = re.findall(r'database_code_amcsd\s+(\d+)', html)
        amcsd_ids = list(dict.fromkeys(amcsd_ids + db_codes))

        existing_ids = {r["cod_id"] for r in results}
        for amid in amcsd_ids:
            if amid in existing_ids:
                continue
            entry = {
                "cod_id": amid, "amcsd_id": amid,
                "formula": "", "title": query.strip().title(),
                "mineral": query.strip().title(), "spacegroup": "",
                "year": "", "has_coordinates": True, "source": "AMCSD",
            }
            # Quick CIF header fetch for metadata
            try:
                cif_url = f"https://rruff.geo.arizona.edu/AMS/download.php?id={amid}.cif&down=cif"
                cif_req = urllib.request.Request(cif_url, headers={"User-Agent": "terrahedral/0.1"})
                with urllib.request.urlopen(cif_req, timeout=4, context=ctx) as cr:
                    hdr = cr.read(2048).decode("utf-8", errors="replace")
                for line in hdr.split("\n"):
                    ls = line.strip()
                    if ls.startswith("_chemical_formula_sum"):
                        entry["formula"] = ls.split(None, 1)[-1].strip("' \"")
                    elif ls.startswith("_chemical_name_mineral"):
                        entry["mineral"] = ls.split(None, 1)[-1].strip("' \"")
                        entry["title"] = entry["mineral"]
                    elif ls.startswith("_symmetry_space_group_name"):
                        entry["spacegroup"] = ls.split(None, 1)[-1].strip("' \"")
                    elif ls.startswith("_journal_year"):
                        entry["year"] = ls.split(None, 1)[-1].strip()
                    elif ls.startswith("_diffrn_radiation_type"):
                        rad = ls.split(None, 1)[-1].strip("' \"").lower()
                        if "neutron" in rad:
                            entry["method"] = "Neutron"
                        elif "synchrotron" in rad or "sr" in rad:
                            entry["method"] = "Synchrotron XRD"
                        else:
                            entry["method"] = "XRD"
                    elif ls.startswith("_journal_name_full"):
                        entry["journal"] = ls.split(None, 1)[-1].strip("' \"")
            except Exception:
                pass
            results.append(entry)
            existing_ids.add(amid)
            if len(results) >= max_results:
                break
    except Exception as e:
        print(f"  [amcsd] Direct search failed: {e}")

    if len(results) >= max_results:
        return results[:max_results]

    # Strategy 3: COD text search filtered to 9xxxxxx
    try:
        params = urllib.parse.urlencode({"text": query.strip(), "format": "urls"})
        url = f"{COD_SEARCH_URL}?{params}"
        req = urllib.request.Request(url, headers={"User-Agent": "terrahedral/0.1"})
        with urllib.request.urlopen(req, timeout=12) as resp:
            text = resp.read().decode("utf-8", errors="replace")
        cod_ids = re.findall(r"/cod/(\d+)\.cif", text)
        cod_ids = list(dict.fromkeys(cod_ids))
        amcsd_cods = [cid for cid in cod_ids if cid.startswith("9")]
        existing_ids = {r["cod_id"] for r in results}
        need = max_results - len(results)
        new_ids = [cid for cid in amcsd_cods if cid not in existing_ids][:need]
        if new_ids:
            extra = _fetch_cod_metadata(new_ids, max_results=need)
            for r in extra:
                r["has_coordinates"] = True
                r["source"] = "AMCSD (via COD)"
            results.extend(extra)
    except Exception as e:
        print(f"  [amcsd] COD fallback failed: {e}")

    return results[:max_results]


# ── Curated mineral index ──
# Common minerals with known-good AMCSD/COD IDs that have full coordinates.
# AMCSD IDs (short, like "00594") are for direct AMCSD download.
# COD IDs (7-digit starting with 9, like "9000594") are the same entry in COD.
_MINERAL_INDEX = [
    {"mineral": "Pyrite", "formula": "Fe S2", "cod_id": "9000594", "amcsd_id": "00594",
     "spacegroup": "P a -3", "year": "1976", "elements": ["Fe", "S"], "method": "XRD"},
    {"mineral": "Magnetite", "formula": "Fe3 O4", "cod_id": "9000926", "amcsd_id": "00926",
     "spacegroup": "F d -3 m", "year": "1915", "elements": ["Fe", "O"], "method": "XRD"},
    {"mineral": "Hematite", "formula": "Fe2 O3", "cod_id": "9000139", "amcsd_id": "00139",
     "spacegroup": "R -3 c", "year": "1970", "elements": ["Fe", "O"], "method": "XRD"},
    {"mineral": "Goethite", "formula": "Fe O (O H)", "cod_id": "9003078", "amcsd_id": "03078",
     "spacegroup": "P b n m", "year": "1970", "elements": ["Fe", "O", "H"], "method": "XRD"},
    {"mineral": "Quartz", "formula": "Si O2", "cod_id": "9000775", "amcsd_id": "00775",
     "spacegroup": "P 32 2 1", "year": "1980", "elements": ["Si", "O"], "method": "XRD"},
    {"mineral": "Calcite", "formula": "Ca C O3", "cod_id": "9000095", "amcsd_id": "00095",
     "spacegroup": "R -3 c", "year": "1993", "elements": ["Ca", "C", "O"], "method": "XRD"},
    {"mineral": "Dolomite", "formula": "Ca Mg (C O3)2", "cod_id": "9000200", "amcsd_id": "00200",
     "spacegroup": "R -3", "year": "1977", "elements": ["Ca", "Mg", "C", "O"], "method": "XRD"},
    {"mineral": "Galena", "formula": "Pb S", "cod_id": "9000030", "amcsd_id": "00030",
     "spacegroup": "F m -3 m", "year": "1962", "elements": ["Pb", "S"], "method": "XRD"},
    {"mineral": "Sphalerite", "formula": "Zn S", "cod_id": "9000104", "amcsd_id": "00104",
     "spacegroup": "F -4 3 m", "year": "1962", "elements": ["Zn", "S"], "method": "XRD"},
    {"mineral": "Chalcopyrite", "formula": "Cu Fe S2", "cod_id": "9000375", "amcsd_id": "00375",
     "spacegroup": "I -4 2 d", "year": "1973", "elements": ["Cu", "Fe", "S"], "method": "XRD"},
    {"mineral": "Arsenopyrite", "formula": "Fe As S", "cod_id": "9000232", "amcsd_id": "00232",
     "spacegroup": "P 21/c", "year": "1981", "elements": ["Fe", "As", "S"], "method": "XRD"},
    {"mineral": "Olivine", "formula": "Mg2 Si O4", "cod_id": "9000319", "amcsd_id": "00319",
     "spacegroup": "P b n m", "year": "1982", "elements": ["Mg", "Si", "O"], "method": "XRD"},
    {"mineral": "Forsterite", "formula": "Mg2 Si O4", "cod_id": "9000319", "amcsd_id": "00319",
     "spacegroup": "P b n m", "year": "1982", "elements": ["Mg", "Si", "O"], "method": "XRD"},
    {"mineral": "Fayalite", "formula": "Fe2 Si O4", "cod_id": "9000176", "amcsd_id": "00176",
     "spacegroup": "P b n m", "year": "1982", "elements": ["Fe", "Si", "O"], "method": "XRD"},
    {"mineral": "Garnet", "formula": "Mg3 Al2 Si3 O12", "cod_id": "9000232", "amcsd_id": "00232",
     "spacegroup": "I a -3 d", "year": "1976", "elements": ["Mg", "Al", "Si", "O"], "method": "XRD"},
    {"mineral": "Pyrope", "formula": "Mg3 Al2 Si3 O12", "cod_id": "9000095", "amcsd_id": "00095",
     "spacegroup": "I a -3 d", "year": "1976", "elements": ["Mg", "Al", "Si", "O"], "method": "XRD"},
    {"mineral": "Rutile", "formula": "Ti O2", "cod_id": "9000318", "amcsd_id": "00318",
     "spacegroup": "P 42/m n m", "year": "1956", "elements": ["Ti", "O"], "method": "XRD"},
    {"mineral": "Corundum", "formula": "Al2 O3", "cod_id": "9000141", "amcsd_id": "00141",
     "spacegroup": "R -3 c", "year": "1962", "elements": ["Al", "O"], "method": "XRD"},
    {"mineral": "Spinel", "formula": "Mg Al2 O4", "cod_id": "9000070", "amcsd_id": "00070",
     "spacegroup": "F d -3 m", "year": "1976", "elements": ["Mg", "Al", "O"], "method": "XRD"},
    {"mineral": "Cuprite", "formula": "Cu2 O", "cod_id": "9000529", "amcsd_id": "00529",
     "spacegroup": "P n -3 m", "year": "1962", "elements": ["Cu", "O"], "method": "XRD"},
    {"mineral": "Cassiterite", "formula": "Sn O2", "cod_id": "9000457", "amcsd_id": "00457",
     "spacegroup": "P 42/m n m", "year": "1965", "elements": ["Sn", "O"], "method": "XRD"},
    {"mineral": "Halite", "formula": "Na Cl", "cod_id": "9000596", "amcsd_id": "00596",
     "spacegroup": "F m -3 m", "year": "1981", "elements": ["Na", "Cl"], "method": "XRD"},
    {"mineral": "Fluorite", "formula": "Ca F2", "cod_id": "9000056", "amcsd_id": "00056",
     "spacegroup": "F m -3 m", "year": "1981", "elements": ["Ca", "F"], "method": "XRD"},
    {"mineral": "Rhodochrosite", "formula": "Mn C O3", "cod_id": "9000597", "amcsd_id": "00597",
     "spacegroup": "R -3 c", "year": "1981", "elements": ["Mn", "C", "O"], "method": "XRD"},
    {"mineral": "Siderite", "formula": "Fe C O3", "cod_id": "9000598", "amcsd_id": "00598",
     "spacegroup": "R -3 c", "year": "1981", "elements": ["Fe", "C", "O"], "method": "XRD"},
    {"mineral": "Birnessite", "formula": "Mn O2", "cod_id": "9001658", "amcsd_id": "01658",
     "spacegroup": "C 2/m", "year": "1992", "elements": ["Mn", "O"], "method": "XRD"},
    {"mineral": "Todorokite", "formula": "Mn O2", "cod_id": "9002530", "amcsd_id": "02530",
     "spacegroup": "P 2/m", "year": "2003", "elements": ["Mn", "O"], "method": "XRD"},
    {"mineral": "Bornite", "formula": "Cu5 Fe S4", "cod_id": "9001136", "amcsd_id": "01136",
     "spacegroup": "F -4 3 m", "year": "1988", "elements": ["Cu", "Fe", "S"], "method": "XRD"},
    {"mineral": "Molybdenite", "formula": "Mo S2", "cod_id": "9000683", "amcsd_id": "00683",
     "spacegroup": "P 63/m m c", "year": "1983", "elements": ["Mo", "S"], "method": "XRD"},
    {"mineral": "Cinnabar", "formula": "Hg S", "cod_id": "9000680", "amcsd_id": "00680",
     "spacegroup": "P 32 2 1", "year": "1965", "elements": ["Hg", "S"], "method": "XRD"},
    {"mineral": "Ilmenite", "formula": "Fe Ti O3", "cod_id": "9000904", "amcsd_id": "00904",
     "spacegroup": "R -3", "year": "1981", "elements": ["Fe", "Ti", "O"], "method": "XRD"},
    {"mineral": "Chromite", "formula": "Fe Cr2 O4", "cod_id": "9000063", "amcsd_id": "00063",
     "spacegroup": "F d -3 m", "year": "1984", "elements": ["Fe", "Cr", "O"], "method": "XRD"},
    {"mineral": "Malachite", "formula": "Cu2 (C O3) (O H)2", "cod_id": "9000694", "amcsd_id": "00694",
     "spacegroup": "P 21/a", "year": "1981", "elements": ["Cu", "C", "O", "H"], "method": "XRD"},
    {"mineral": "Azurite", "formula": "Cu3 (C O3)2 (O H)2", "cod_id": "9000250", "amcsd_id": "00250",
     "spacegroup": "P 21/c", "year": "1972", "elements": ["Cu", "C", "O", "H"], "method": "XRD"},
    {"mineral": "Barite", "formula": "Ba S O4", "cod_id": "9001635", "amcsd_id": "01635",
     "spacegroup": "P n m a", "year": "1992", "elements": ["Ba", "S", "O"], "method": "XRD"},
    {"mineral": "Gypsum", "formula": "Ca S O4 (H2 O)2", "cod_id": "9000140", "amcsd_id": "00140",
     "spacegroup": "A 2/a", "year": "1982", "elements": ["Ca", "S", "O", "H"], "method": "XRD"},
    {"mineral": "Apatite", "formula": "Ca5 (P O4)3 F", "cod_id": "9002232", "amcsd_id": "02232",
     "spacegroup": "P 63/m", "year": "1969", "elements": ["Ca", "P", "O", "F"], "method": "XRD"},
    {"mineral": "Tourmaline", "formula": "Na Mg3 Al6 Si6 O18 (B O3)3 (O H)4",
     "cod_id": "9001505", "amcsd_id": "01505",
     "spacegroup": "R 3 m", "year": "1986", "elements": ["Na", "Mg", "Al", "Si", "O", "B", "H"], "method": "XRD"},
    {"mineral": "Muscovite", "formula": "K Al3 Si3 O10 (O H)2", "cod_id": "9000834",
     "amcsd_id": "00834", "spacegroup": "C 2/c", "year": "1984",
     "elements": ["K", "Al", "Si", "O", "H"], "method": "XRD"},
    {"mineral": "Biotite", "formula": "K Mg3 Al Si3 O10 (O H)2", "cod_id": "9000307",
     "amcsd_id": "00307", "spacegroup": "C 2/m", "year": "1975",
     "elements": ["K", "Mg", "Al", "Si", "O", "H"], "method": "XRD"},
    {"mineral": "Wustite", "formula": "Fe O", "cod_id": "9002659", "amcsd_id": "02659",
     "spacegroup": "F m -3 m", "year": "1984", "elements": ["Fe", "O"], "method": "XRD"},
    {"mineral": "Periclase", "formula": "Mg O", "cod_id": "9000501", "amcsd_id": "00501",
     "spacegroup": "F m -3 m", "year": "1976", "elements": ["Mg", "O"], "method": "XRD"},
    {"mineral": "Zincite", "formula": "Zn O", "cod_id": "9004178", "amcsd_id": "04178",
     "spacegroup": "P 63 m c", "year": "1969", "elements": ["Zn", "O"], "method": "XRD"},
    {"mineral": "Tenorite", "formula": "Cu O", "cod_id": "9000484", "amcsd_id": "00484",
     "spacegroup": "C 2/c", "year": "1970", "elements": ["Cu", "O"], "method": "XRD"},
    {"mineral": "Manganosite", "formula": "Mn O", "cod_id": "9006395", "amcsd_id": "06395",
     "spacegroup": "F m -3 m", "year": "2008", "elements": ["Mn", "O"], "method": "XRD"},
    {"mineral": "Pyrolusite", "formula": "Mn O2", "cod_id": "9000510", "amcsd_id": "00510",
     "spacegroup": "P 42/m n m", "year": "1965", "elements": ["Mn", "O"], "method": "XRD"},
    {"mineral": "Hausmannite", "formula": "Mn3 O4", "cod_id": "9000605", "amcsd_id": "00605",
     "spacegroup": "I 41/a m d", "year": "1968", "elements": ["Mn", "O"], "method": "XRD"},
    {"mineral": "Braunite", "formula": "Mn2+ Mn3+6 Si O12", "cod_id": "9001238",
     "amcsd_id": "01238", "spacegroup": "I 41/a c d", "year": "1986", "elements": ["Mn", "Si", "O"], "method": "XRD"},
]


def _fetch_cod_metadata(
    cod_ids: list[str],
    *,
    max_results: int = 8,
) -> list[dict]:
    """Fetch basic metadata for COD entries by reading the first few lines of each CIF."""
    out = []
    for cod_id in cod_ids[:max_results]:
        entry = {"cod_id": cod_id, "formula": "", "title": "", "mineral": ""}
        try:
            url = f"https://www.crystallography.net/cod/{cod_id}.cif"
            req = urllib.request.Request(url, headers={"User-Agent": "terrahedral/0.1"})
            with urllib.request.urlopen(req, timeout=8) as resp:
                # Read just the header (first 4KB) to extract metadata
                header = resp.read(4096).decode("utf-8", errors="replace")

            for line in header.split("\n"):
                line_s = line.strip()
                if line_s.startswith("_chemical_formula_sum"):
                    entry["formula"] = line_s.split(None, 1)[-1].strip("' \"")
                elif line_s.startswith("_chemical_name_mineral"):
                    entry["mineral"] = line_s.split(None, 1)[-1].strip("' \"")
                elif line_s.startswith("_chemical_name_systematic"):
                    entry["title"] = line_s.split(None, 1)[-1].strip("' \"")
                elif line_s.startswith("_publ_author_name") or line_s.startswith("loop_"):
                    # Past the header fields we care about
                    if entry["formula"] or entry["title"]:
                        break
        except Exception:
            pass

        # Use mineral name as title if no systematic name
        if not entry["title"] and entry["mineral"]:
            entry["title"] = entry["mineral"]

        out.append(entry)

    return out


def fetch_pdb(
    pdb_id: str,
    *,
    cache_dir: Optional[str | Path] = None,
    force: bool = False,
) -> Path:
    """
    Download an mmCIF file from RCSB PDB.

    Parameters
    ----------
    pdb_id : 4-character PDB identifier (e.g. "6GTL", "1W56")
    cache_dir : directory to store downloaded files (default: ~/.cache/terrahedral)
    force : re-download even if cached

    Returns
    -------
    Path to the local .cif file.
    """
    pdb_id = pdb_id.strip().upper()
    if len(pdb_id) != 4:
        raise ValueError(f"PDB ID must be 4 characters, got: {pdb_id!r}")

    cache = Path(cache_dir) if cache_dir else _DEFAULT_CACHE
    _ensure_dir(cache)

    local_path = cache / f"{pdb_id}.cif"

    if local_path.exists() and not force:
        return local_path

    url = RCSB_URL.format(pdb_id=pdb_id)
    try:
        urllib.request.urlretrieve(url, local_path)
    except urllib.error.HTTPError as e:
        if e.code == 404:
            raise ValueError(f"PDB entry '{pdb_id}' not found on RCSB") from e
        raise ConnectionError(f"RCSB download failed ({e.code}): {url}") from e
    except urllib.error.URLError as e:
        raise ConnectionError(f"Could not reach RCSB: {e.reason}") from e

    return local_path


def fetch_cod(
    cod_id: str,
    *,
    cache_dir: Optional[str | Path] = None,
    force: bool = False,
) -> Path:
    """
    Download a small-molecule CIF from the Crystallography Open Database.

    Parameters
    ----------
    cod_id : numeric COD identifier (e.g. "4070506")
    cache_dir : directory to store downloaded files (default: ~/.cache/terrahedral)
    force : re-download even if cached

    Returns
    -------
    Path to the local .cif file.
    """
    cod_id = cod_id.strip()
    if not cod_id.isdigit():
        raise ValueError(f"COD ID must be numeric, got: {cod_id!r}")

    cache = Path(cache_dir) if cache_dir else _DEFAULT_CACHE
    _ensure_dir(cache)

    local_path = cache / f"COD_{cod_id}.cif"

    if local_path.exists() and not force:
        return local_path

    url = COD_URL.format(cod_id=cod_id)
    try:
        urllib.request.urlretrieve(url, local_path)
    except urllib.error.HTTPError as e:
        if e.code == 404:
            raise ValueError(f"COD entry '{cod_id}' not found") from e
        raise ConnectionError(f"COD download failed ({e.code}): {url}") from e
    except urllib.error.URLError as e:
        raise ConnectionError(f"Could not reach COD: {e.reason}") from e

    return local_path


def fetch_amcsd(
    amcsd_id: str,
    *,
    cache_dir: Optional[str | Path] = None,
    force: bool = False,
) -> Path:
    """
    Download a CIF from the American Mineralogist Crystal Structure Database.

    Parameters
    ----------
    amcsd_id : AMCSD numeric identifier (e.g. "00594")
    cache_dir : directory to store downloaded files
    force : re-download even if cached

    Returns
    -------
    Path to the local .cif file.
    """
    import ssl

    amcsd_id = amcsd_id.strip()
    if not amcsd_id.isdigit():
        raise ValueError(f"AMCSD ID must be numeric, got: {amcsd_id!r}")

    cache = Path(cache_dir) if cache_dir else _DEFAULT_CACHE
    _ensure_dir(cache)

    local_path = cache / f"AMCSD_{amcsd_id}.cif"

    if local_path.exists() and not force:
        return local_path

    url = f"https://rruff.geo.arizona.edu/AMS/download.php?id={amcsd_id}.cif&down=cif"
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "terrahedral/0.1"})
        ctx = ssl.create_default_context()
        ctx.check_hostname = False
        ctx.verify_mode = ssl.CERT_NONE
        with urllib.request.urlopen(req, timeout=15, context=ctx) as resp:
            data = resp.read()
        local_path.write_bytes(data)
    except urllib.error.HTTPError as e:
        if e.code == 404:
            raise ValueError(f"AMCSD entry '{amcsd_id}' not found") from e
        raise ConnectionError(f"AMCSD download failed ({e.code}): {url}") from e
    except urllib.error.URLError as e:
        raise ConnectionError(f"Could not reach AMCSD: {e.reason}") from e

    return local_path
