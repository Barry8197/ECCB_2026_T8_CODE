# Session 1: Introduction to Biological Networks and Knowledge Graphs

**New here?** [`session-story.md`](session-story.md) is a one-page narrative of what
this session does and what it finds — read that first, then the notebooks.

## Learning Goals
- Understand core network concepts and graph terminology
- Understand what a biomedical knowledge graph is, and how it differs from an
  inferred network — by building one of each over the same genes
- Build and query a knowledge graph using NetworkX
- Recognise annotation sparsity and identifier mismatches — the two failure modes
  that matter most in practice
- Combine curated and inferred evidence to separate corroboration from candidates

## Notebooks

| Notebook | Agenda slot | Contents |
|---|---|---|
| [`part-1-networks-and-knowledge-graphs.ipynb`](part-1-networks-and-knowledge-graphs.ipynb) | 09:10–09:40 | Graph vocabulary, inferred vs curated networks, the schema, data sources and licences, and building a co-expression network alongside the graph (§7) |
| [`part-2-building-a-knowledge-graph.ipynb`](part-2-building-a-knowledge-graph.ipynb) | 09:40–10:10 | Building the graph in NetworkX, degree and hubs, association scores, sparsification, annotation sparsity |
| [`part-3-querying-the-knowledge-graph.ipynb`](part-3-querying-the-knowledge-graph.ipynb) | 10:10–10:45 | **Practical.** ICD-10 coverage and ontology climbing (§A), shared-gene projection, separating cause from treatment, community detection (§B), bridging to the Session 2 omics data (§C, stretch), curated and inferred networks together (§D, closing) |

These are the **participant** notebooks: the interesting lines are replaced with
`### YOUR CODE HERE ###` for you to fill in. They are generated from the solution
notebooks in the instructors' repository — if something needs changing, it has to
change there and be regenerated, so edits made here will be overwritten.

## Layout

```text
session-story.md    the one-page narrative of the session - start here
part-*.ipynb        the three notebooks
s1_helpers.py       graph utilities, mostly adapted from co-expression material
concept-glossary.md every term used in the session, defined
data-prep/          developer scripts (build_kg_data.py,
                    build_coexpression_data.py, make_images.py) - not run by
                    participants; they write to ../../data/session-1-data/
images/             diagrams used by the notebooks
```

The data itself lives at the repository root, in `data/session-1-data/`.

## Data

**The knowledge graph.** Breast-cancer-focused, with 24 comparison diseases
spanning cancers, autoimmune, neurological and cardiometabolic conditions:
**881 nodes** (760 gene, 90 disease, 31 ICD-10) and **1,768 edges**
(`associated_with`, `is_a`, `maps_to`), plus **3,422 evidence rows** splitting
each gene–disease edge by kind of evidence. Built from
[Open Targets](https://platform.opentargets.org/) (CC0) and the
[MONDO Disease Ontology](https://mondo.monarchinitiative.org/) (CC BY 4.0).

**The expression matrix.** `coexpr_expression.csv.gz` — 500 patients × 737 genes
of TCGA-BRCA transcriptomics, subset to genes that are already graph nodes. Part 1
§7 builds an inferred co-expression network from it, and Part 3 §D puts the two
networks together.

All five files live in [`../../data/session-1-data/`](../../data/session-1-data/)
and are committed to this repository — together they are 1.3 MB, so Session 1
needs no separate download. See
[`../../data/session-1-data/README.md`](../../data/session-1-data/README.md)
for provenance and licences.

Genes are keyed by Ensembl gene ID, matching the TCGA-BRCA matrix used in
Session 2, so the two join natively.

## Requirements

`pandas`, `numpy`, `networkx`, `matplotlib`, `seaborn`, `scipy` (see
`environment.yml` at the repository root).

Everything runs from the committed data **except Section C of Part 3**, which
needs the full ~900 MB Session 2 omics pickle. Most participants will not have
it: the notebook catches the missing file, prints a skip message and carries on,
and Section D — the closing section that puts the curated and inferred networks
together — deliberately uses the committed 1.1 MB subset instead so that it runs
for everyone.
