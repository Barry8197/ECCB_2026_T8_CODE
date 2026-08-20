"""
Build the Session 1 knowledge-graph input files.

This is a DEVELOPER script, not something workshop participants run. It downloads
public reference data, cuts it down to a small breast-cancer-focused subset, and
writes a handful of CSVs that are small enough to commit to the repository.
Participants just read those CSVs.

Data sources
------------
Open Targets Platform  (CC0 1.0 - public domain, redistributable)
    https://platform-docs.opentargets.org/licence
    - disease.parquet                 disease entities + ontology hierarchy
    - target.parquet                  gene entities (Ensembl ID -> symbol)
    - association_overall_direct/     gene <-> disease association scores
    - association_by_datatype_direct/ the same scores split by evidence type

MONDO Disease Ontology  (CC BY 4.0 - redistributable with attribution)
    https://mondo.monarchinitiative.org/
    - mondo.sssom.tsv                 MONDO -> ICD-10 cross references

Why not DisGeNET: as of 2024 it is commercially licensed. The free academic
licence does not permit downloading the full database (search results only) and
requires per-user registration with ~7 day approval, so it cannot be shipped in
a public teaching repository.

Why these two together: Open Targets keys genes by Ensembl ID (ENSG...), which is
the same identifier space as the TCGA-BRCA transcriptomics matrix used elsewhere
in the workshop, so the knowledge graph joins to the omics data natively. Open
Targets' own ICD-10 coverage is sparse, so MONDO supplies the ICD-10 mapping.

Usage
-----
    python build_kg_data.py                  # default: breast cancer focus
    python build_kg_data.py --top-genes 50   # more genes per disease
    python build_kg_data.py --keep-cache     # keep the multi-GB downloads

Outputs (written to ../data/)
-----------------------------
    kg_nodes.csv        id, type, name, extra
    kg_edges.csv        source, target, type, weight, evidence
    kg_evidence.csv     source, target, datatype, weight, evidence
    icd10_map.csv       mondo_id, icd10_code, icd10_label
    README.md           provenance + licence note for the generated files
"""

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path
from urllib.request import urlopen
from urllib.error import URLError, HTTPError

import pandas as pd

try:
    import pyarrow.parquet as pq
except ImportError:  # pragma: no cover
    sys.exit("pyarrow is required:  pip install pyarrow")


# --------------------------------------------------------------------------
# Configuration
# --------------------------------------------------------------------------

OT_RELEASE = "latest"  # pin to e.g. "26.06" for a fully reproducible build
OT_BASE = f"https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/{OT_RELEASE}/output"
MONDO_SSSOM = (
    "https://raw.githubusercontent.com/monarch-initiative/mondo/master/"
    "src/ontology/mappings/mondo.sssom.tsv"
)

# The disease the workshop cohort is about. Everything is built outwards from here.
SEED_DISEASE = "MONDO_0007254"  # breast cancer -> ICD10CM:C50

# A few comparison diseases so the graph is not a single star and students can ask
# "which genes are shared between these?". Chosen to be well-annotated and to span
# both other cancers and non-cancer conditions.
# Grouped so the projection in Part 3 has real community structure to find:
# without non-cancer diseases in here, every disease shares genes with every
# other one and "which diseases resemble each other" has only one answer.
COMPARISON_DISEASES = {
    # --- other cancers ---
    "MONDO_0008170": "ovarian cancer",
    "MONDO_0005061": "lung adenocarcinoma",  # a well-annotated non-breast cancer
    "MONDO_0008315": "prostate cancer",
    "MONDO_0005575": "colorectal cancer",
    "MONDO_0005184": "pancreatic ductal adenocarcinoma",
    "MONDO_0005105": "melanoma",
    "MONDO_0018177": "glioblastoma",
    "MONDO_0004950": "gastric carcinoma",
    "MONDO_0018874": "acute myeloid leukemia",
    # --- autoimmune / inflammatory ---
    "MONDO_0008383": "rheumatoid arthritis",
    "MONDO_0005011": "Crohn disease",
    "MONDO_0005101": "ulcerative colitis",
    "MONDO_0005301": "multiple sclerosis",
    "MONDO_0007915": "systemic lupus erythematosus",
    "MONDO_0004979": "asthma",
    "MONDO_0005083": "psoriasis",
    "MONDO_0005147": "type 1 diabetes mellitus",
    # --- neurological ---
    "MONDO_0004975": "Alzheimer disease",
    "MONDO_0005180": "Parkinson disease",
    "MONDO_0005090": "schizophrenia",
    # --- cardiometabolic / other ---
    "MONDO_0005148": "type 2 diabetes mellitus",
    "MONDO_0005010": "coronary artery disorder",
    "MONDO_0005300": "chronic kidney disease",
    "MONDO_0005298": "osteoporosis",
}

SCORE_THRESHOLD = 0.05  # drop very weak associations; keeps the graph readable


# --------------------------------------------------------------------------
# Download helpers
# --------------------------------------------------------------------------

def as_list(value) -> list:
    """
    Coerce a parquet list-column cell to a plain list.

    Columns like `parents` and `descendants` come back as numpy arrays, for which
    the usual `value or []` idiom raises "truth value of an array is ambiguous".
    """
    if value is None:
        return []
    return list(value)


def download(url: str, dest: Path, desc: str) -> Path:
    """Download `url` to `dest`, skipping if the file is already present."""
    if dest.exists() and dest.stat().st_size > 0:
        print(f"  [cached] {desc}")
        return dest

    dest.parent.mkdir(parents=True, exist_ok=True)
    print(f"  [get]    {desc} ...", end="", flush=True)
    try:
        with urlopen(url) as response, open(dest, "wb") as handle:
            shutil.copyfileobj(response, handle)
    except (URLError, HTTPError) as exc:
        dest.unlink(missing_ok=True)
        sys.exit(f"\nFailed to download {url}\n  {exc}")
    print(f" {dest.stat().st_size / 1e6:.1f} MB")
    return dest


def list_parquet_parts(directory_url: str) -> list[str]:
    """Scrape an Open Targets FTP directory listing for its parquet part files."""
    import re

    with urlopen(directory_url) as response:
        html = response.read().decode("utf-8", errors="replace")
    parts = sorted(set(re.findall(r'part-\d+[^"<>\s]*\.parquet', html)))
    if not parts:
        sys.exit(f"No parquet parts found at {directory_url}")
    return parts


# --------------------------------------------------------------------------
# Build steps
# --------------------------------------------------------------------------

def load_diseases(cache: Path) -> pd.DataFrame:
    """Open Targets disease entities, including the ontology hierarchy."""
    path = download(f"{OT_BASE}/disease/disease.parquet", cache / "disease.parquet",
                    "Open Targets disease entities")
    df = pq.read_table(path).to_pandas()
    print(f"           {len(df):,} diseases")
    return df


def load_targets(cache: Path) -> pd.DataFrame:
    """Open Targets gene entities -> Ensembl ID to symbol lookup."""
    parts = list_parquet_parts(f"{OT_BASE}/target/")
    frames = []
    for part in parts:
        path = download(f"{OT_BASE}/target/{part}", cache / "target" / part,
                        f"Open Targets targets {part[:11]}")
        table = pq.read_table(path, columns=["id", "approvedSymbol", "approvedName",
                                             "biotype"])
        frames.append(table.to_pandas())
    df = pd.concat(frames, ignore_index=True)
    print(f"           {len(df):,} genes")
    return df


def load_associations(cache: Path, wanted: set[str], keep_cache: bool) -> pd.DataFrame:
    """
    Gene <-> disease association scores, filtered to the diseases we care about.

    The full dataset is ~1 GB across 14 parquet parts. Each part is downloaded,
    filtered down to `wanted`, and then discarded unless --keep-cache is set, so
    peak disk usage stays modest.
    """
    parts = list_parquet_parts(f"{OT_BASE}/association_overall_direct/")
    frames = []

    for i, part in enumerate(parts, 1):
        path = download(
            f"{OT_BASE}/association_overall_direct/{part}",
            cache / "assoc" / part,
            f"associations part {i}/{len(parts)}",
        )
        table = pq.read_table(
            path, columns=["diseaseId", "targetId", "associationScore", "evidenceCount"]
        )
        chunk = table.to_pandas()
        frames.append(chunk[chunk["diseaseId"].isin(wanted)])

        if not keep_cache:
            path.unlink(missing_ok=True)  # reclaim ~75 MB immediately

    df = pd.concat(frames, ignore_index=True)
    print(f"           {len(df):,} associations for {df['diseaseId'].nunique()} diseases")
    return df


def load_datatype_evidence(cache: Path, wanted: set[str], keep_cache: bool) -> pd.DataFrame:
    """
    Association scores broken down by *type* of evidence.

    `association_overall_direct` collapses everything into one number, which
    cannot distinguish "this gene causes the disease" from "this gene is hit by
    a drug used to treat the disease". This dataset keeps them apart:

        genetic_association, somatic_mutation, known_drug, literature,
        rna_expression, animal_model, affected_pathway, genetic_literature

    That distinction is the whole point of the tubulin exercise in Part 3 - the
    difference between BRCA2 and TUBB is invisible in the overall score and
    obvious here.

    The `timeseries` column is a nested per-year struct and accounts for most of
    the file size, so it is never read.
    """
    parts = list_parquet_parts(f"{OT_BASE}/association_by_datatype_direct/")
    frames = []

    for i, part in enumerate(parts, 1):
        path = download(
            f"{OT_BASE}/association_by_datatype_direct/{part}",
            cache / "datatype" / part,
            f"evidence by datatype part {i}/{len(parts)}",
        )
        table = pq.read_table(path, columns=["diseaseId", "targetId", "aggregationType",
                                             "aggregationValue", "associationScore",
                                             "evidenceCount"])
        chunk = table.to_pandas()
        chunk = chunk[(chunk["diseaseId"].isin(wanted))
                      & (chunk["aggregationType"] == "datatypeId")]
        frames.append(chunk)

        if not keep_cache:
            path.unlink(missing_ok=True)

    df = pd.concat(frames, ignore_index=True)
    df = df.rename(columns={"aggregationValue": "datatype"})
    print(f"           {len(df):,} evidence rows across "
          f"{df['datatype'].nunique()} evidence types")
    return df


def load_icd10_map(cache: Path) -> pd.DataFrame:
    """MONDO -> ICD-10 cross references, from the SSSOM mapping file."""
    path = download(MONDO_SSSOM, cache / "mondo.sssom.tsv", "MONDO SSSOM mappings")
    df = pd.read_csv(path, sep="\t", comment="#", dtype=str, low_memory=False)

    icd = df[df["object_id"].astype(str).str.startswith("ICD10")].copy()
    # Open Targets writes MONDO ids with an underscore; SSSOM uses a colon.
    icd["mondo_id"] = icd["subject_id"].str.replace(":", "_", regex=False)
    icd = icd.rename(columns={"object_id": "icd10_code", "object_label": "icd10_label"})

    out = icd[["mondo_id", "icd10_code", "icd10_label"]].drop_duplicates()
    print(f"           {len(out):,} MONDO->ICD-10 mappings")
    return out


def build_graph_tables(diseases, targets, associations, icd10, top_genes):
    """Assemble typed node and edge tables from the filtered source data."""
    nodes: list[dict] = []
    edges: list[dict] = []

    disease_by_id = diseases.set_index("id")

    # --- Disease nodes + the ontology hierarchy that connects them -----------
    disease_ids = set(associations["diseaseId"])
    for disease_id in sorted(disease_ids):
        if disease_id not in disease_by_id.index:
            continue
        row = disease_by_id.loc[disease_id]
        nodes.append({
            "id": disease_id,
            "type": "disease",
            "name": row["name"],
            "extra": "; ".join(as_list(row["therapeuticAreas"])),
        })

        # is_a edges, kept inside our subset so the graph stays connected
        for parent in as_list(row["parents"]):
            if parent in disease_ids:
                edges.append({
                    "source": disease_id, "target": parent,
                    "type": "is_a", "weight": "", "evidence": "",
                })

    # --- Gene nodes + gene->disease association edges -----------------------
    symbol_by_id = targets.set_index("id")
    kept_genes: set[str] = set()

    # Take the strongest associations per disease so every disease keeps its own
    # top genes rather than being swamped by the best-studied disease.
    ranked = associations[associations["associationScore"] >= SCORE_THRESHOLD]
    ranked = (
        ranked.sort_values("associationScore", ascending=False)
        .groupby("diseaseId", group_keys=False)
        .head(top_genes)
    )

    for row in ranked.itertuples(index=False):
        kept_genes.add(row.targetId)
        edges.append({
            "source": row.targetId,
            "target": row.diseaseId,
            "type": "associated_with",
            "weight": round(float(row.associationScore), 4),
            "evidence": int(row.evidenceCount) if pd.notna(row.evidenceCount) else 0,
        })

    for gene_id in sorted(kept_genes):
        if gene_id in symbol_by_id.index:
            gene = symbol_by_id.loc[gene_id]
            symbol, name = gene["approvedSymbol"], gene["approvedName"]
        else:
            symbol, name = gene_id, ""
        nodes.append({"id": gene_id, "type": "gene", "name": symbol, "extra": name})

    # --- ICD-10 nodes: the clinical anchor for each disease -----------------
    relevant_icd = icd10[icd10["mondo_id"].isin(disease_ids)]
    for row in relevant_icd.itertuples(index=False):
        nodes.append({
            "id": row.icd10_code,
            "type": "icd10",
            "name": row.icd10_label,
            "extra": "",
        })
        edges.append({
            "source": row.mondo_id, "target": row.icd10_code,
            "type": "maps_to", "weight": "", "evidence": "",
        })

    node_df = pd.DataFrame(nodes).drop_duplicates(subset="id").reset_index(drop=True)
    edge_df = pd.DataFrame(edges).drop_duplicates(
        subset=["source", "target", "type"]
    ).reset_index(drop=True)
    return node_df, edge_df


PROVENANCE = """\
# Session 1 - generated knowledge graph data

These files are **generated**. Do not edit them by hand - rerun
`../data-prep/build_kg_data.py` instead.

| File | Contents |
|------|----------|
| `kg_nodes.csv` | `id`, `type` (gene / disease / icd10), `name`, `extra` |
| `kg_edges.csv` | `source`, `target`, `type` (associated_with / is_a / maps_to), `weight`, `evidence` |
| `kg_evidence.csv` | `source`, `target`, `datatype`, `weight`, `evidence` - the same gene-disease edges split by *kind* of evidence |
| `icd10_map.csv` | MONDO id to ICD-10 code and label |

## Sources and licences

- **Open Targets Platform** - gene entities, disease entities and ontology
  hierarchy, gene-disease association scores.
  Licensed [CC0 1.0](https://platform-docs.opentargets.org/licence) (public
  domain), so redistribution here is unrestricted.
- **MONDO Disease Ontology** - MONDO to ICD-10 cross references.
  Licensed [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/);
  attribution is given here.

Genes are keyed by **Ensembl gene ID (`ENSG...`)**, matching the identifier space
of the TCGA-BRCA transcriptomics matrix used in Session 2, so gene lists produced
by the multi-omics work can be looked up in this graph directly.
"""


# --------------------------------------------------------------------------
# Entry point
# --------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--top-genes", type=int, default=30,
                        help="strongest gene associations to keep per disease (default: 30)")
    parser.add_argument("--out", type=Path, default=Path(__file__).parent.parent / "data",
                        help="output directory for the generated CSVs")
    parser.add_argument("--cache", type=Path, default=Path(__file__).parent / ".cache",
                        help="where to stash raw downloads")
    parser.add_argument("--keep-cache", action="store_true",
                        help="keep the ~1 GB association downloads instead of deleting them")
    args = parser.parse_args()

    print("Building Session 1 knowledge graph data\n")

    print("1. Reference data")
    diseases = load_diseases(args.cache)
    icd10 = load_icd10_map(args.cache)

    # Expand the seed disease to its subtypes so the hierarchy has something in it.
    seed_row = diseases[diseases["id"] == SEED_DISEASE]
    if seed_row.empty:
        sys.exit(f"Seed disease {SEED_DISEASE} not present in this Open Targets release.")
    descendants = as_list(seed_row.iloc[0]["descendants"])
    wanted = {SEED_DISEASE, *descendants, *COMPARISON_DISEASES}
    print(f"\n2. Disease scope: {SEED_DISEASE} + {len(descendants)} subtypes "
          f"+ {len(COMPARISON_DISEASES)} comparators = {len(wanted)} diseases")

    print("\n3. Gene annotations")
    targets = load_targets(args.cache)

    print("\n4. Associations (large - filtered while streaming)")
    associations = load_associations(args.cache, wanted, args.keep_cache)

    print("\n5. Evidence broken down by type (large - filtered while streaming)")
    evidence = load_datatype_evidence(args.cache, wanted, args.keep_cache)

    print("\n6. Assembling graph tables")
    nodes, edges = build_graph_tables(diseases, targets, associations, icd10,
                                      args.top_genes)

    # Keep evidence only for gene-disease pairs that actually survived into the
    # graph, so the side table cannot describe edges that are not there.
    kept_pairs = set(
        edges.loc[edges["type"] == "associated_with", ["source", "target"]]
        .itertuples(index=False, name=None)
    )
    evidence = evidence[
        [(t, d) in kept_pairs
         for t, d in zip(evidence["targetId"], evidence["diseaseId"])]
    ].copy()
    evidence = evidence.rename(columns={"targetId": "source", "diseaseId": "target",
                                        "associationScore": "weight",
                                        "evidenceCount": "evidence"})
    evidence["weight"] = evidence["weight"].round(4)
    evidence["evidence"] = evidence["evidence"].fillna(0).astype(int)
    evidence = evidence[["source", "target", "datatype", "weight", "evidence"]]
    evidence = evidence.sort_values(["target", "source", "datatype"])

    args.out.mkdir(parents=True, exist_ok=True)
    nodes.to_csv(args.out / "kg_nodes.csv", index=False)
    edges.to_csv(args.out / "kg_edges.csv", index=False)
    evidence.to_csv(args.out / "kg_evidence.csv", index=False)
    icd10[icd10["mondo_id"].isin(set(nodes["id"]))].to_csv(
        args.out / "icd10_map.csv", index=False
    )
    (args.out / "README.md").write_text(PROVENANCE)

    print(f"\nWrote to {args.out}/")
    print(f"  kg_nodes.csv   {len(nodes):,} nodes")
    for node_type, count in nodes["type"].value_counts().items():
        print(f"                   {count:>6,}  {node_type}")
    print(f"  kg_edges.csv   {len(edges):,} edges")
    for edge_type, count in edges["type"].value_counts().items():
        print(f"                   {count:>6,}  {edge_type}")
    print(f"  kg_evidence.csv  {len(evidence):,} rows")
    for datatype, count in evidence["datatype"].value_counts().items():
        print(f"                   {count:>6,}  {datatype}")

    if not args.keep_cache:
        print("\n(association downloads discarded; pass --keep-cache to retain them)")


if __name__ == "__main__":
    main()
