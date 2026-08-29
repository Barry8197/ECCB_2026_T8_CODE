"""
Build the Session 1 co-expression input file.

This is a DEVELOPER script, not something workshop participants run. It reads the
full workshop multi-omics bundle (~900 MB, not in this repository), takes the
transcriptomics view, and cuts it down to just the genes that are already nodes
in the knowledge graph. The result is small enough to commit, so Part 1 can build
a real inferred network without anyone downloading anything.

Why this subset
---------------
Part 1 contrasts an *inferred* network (computed from measurements) with a
*curated* knowledge graph (read from recorded facts). That contrast is much
sharper when both networks describe the **same genes**: any difference is then a
difference in how the edges were obtained, not a difference in which genes were
looked at.

So we intersect the transcriptomics matrix with `kg_nodes.csv`, keeping the gene
nodes that appear in both. 737 of the graph's 760 genes have expression measured.

The identifier gotcha
---------------------
TCGA carries an annotation-release suffix on its Ensembl IDs
(`ENSG00000012048.23`); Open Targets does not (`ENSG00000012048`). Joining the
two without stripping the suffix matches exactly nothing, silently. We strip it
here, and Part 3 makes participants hit the same wall deliberately.

Usage
-----
    python build_coexpression_data.py
    python build_coexpression_data.py --omics /path/to/omics.pkl

Output (written to <repo>/data/session-1-data/)
-----------------------------------------------
    coexpr_expression.csv.gz    patients x genes, log2 normalised expression
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
DATA_DIR = HERE.parents[2] / "data" / "session-1-data"
REPO_ROOT = HERE.parents[2]

# The bundle lives outside the repository because it is ~900 MB. This is where it
# sits in the workshop's own layout; override with --omics if yours differs.
DEFAULT_OMICS = REPO_ROOT.parent / "large_data" / "omics.pkl"

VIEW = "transcriptomics"
OUT_NAME = "coexpr_expression.csv.gz"

# Expression is log2-scale, so four decimals is far below measurement precision
# and roughly halves the file compared with full float64 repr.
DECIMALS = 4


def strip_version(ids) -> list[str]:
    """ENSG00000012048.23 -> ENSG00000012048."""
    return [str(i).split(".")[0] for i in ids]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--omics", type=Path, default=DEFAULT_OMICS,
                        help=f"path to the multi-omics pickle (default: {DEFAULT_OMICS})")
    parser.add_argument("--data-dir", type=Path, default=DATA_DIR,
                        help="where to write the output (default: <repo>/data/session-1-data)")
    args = parser.parse_args()

    if not args.omics.exists():
        raise SystemExit(
            f"No omics bundle at {args.omics}\n"
            f"It is ~900 MB and deliberately not committed. Pass --omics with its "
            f"real location."
        )

    nodes_path = args.data_dir / "kg_nodes.csv"
    if not nodes_path.exists():
        raise SystemExit(
            f"No {nodes_path}. Run build_kg_data.py first - this script subsets "
            f"the expression matrix to genes that are already graph nodes."
        )

    nodes = pd.read_csv(nodes_path)
    kg_genes = list(nodes.loc[nodes["type"] == "gene", "id"])

    blob = pd.read_pickle(args.omics)
    if VIEW not in blob:
        raise SystemExit(f"{args.omics} has no '{VIEW}' view. Keys: {list(blob)}")

    expression = blob[VIEW]
    versioned = list(expression.columns)
    expression.columns = strip_version(versioned)

    # Stripping the suffix can collide: pseudo-autosomal genes appear once per sex
    # chromosome with the same base ID. Duplicate columns would silently corrupt
    # every correlation computed from them, so fail loudly rather than guess.
    duplicated = expression.columns[expression.columns.duplicated()].unique()
    if len(duplicated):
        raise SystemExit(
            f"{len(duplicated)} duplicate gene ids after stripping the version "
            f"suffix, e.g. {list(duplicated[:5])}. Resolve before continuing."
        )

    keep = [g for g in kg_genes if g in expression.columns]
    missing = len(kg_genes) - len(keep)

    subset = expression.loc[:, keep].round(DECIMALS)
    subset.index.name = "patient_id"

    out_path = args.data_dir / OUT_NAME
    subset.to_csv(out_path, compression="gzip")

    print(f"knowledge graph gene nodes : {len(kg_genes)}")
    print(f"  measured in {VIEW:<14}: {len(keep)}")
    print(f"  not measured               : {missing}")
    print(f"patients                     : {subset.shape[0]}")
    print(f"\nwrote {out_path}  ({out_path.stat().st_size / 1e6:.1f} MB)")


if __name__ == "__main__":
    main()
