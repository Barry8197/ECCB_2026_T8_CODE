# Session 1: Introduction to Biological Networks and Knowledge Graphs

**New here?** [`session-story.md`](session-story.md) is a one-page narrative of what
this session does and what it finds — read that first, then the notebooks.

## Learning Goals
- Understand core network concepts and graph terminology
- Understand what a biomedical knowledge graph is, and how it differs from an
  computed network — by building one of each over the same genes
- Build and query a knowledge graph using NetworkX
- Recognise annotation sparsity and identifier mismatches — the two failure modes
  that matter most in practice
- Combine curated and computed evidence to separate corroboration from candidates

## Notebooks

| Notebook | Agenda slot | Contents |
|---|---|---|
| [`part-1-networks-and-knowledge-graphs.ipynb`](part-1-networks-and-knowledge-graphs.ipynb) | 09:10–09:40 | Graph vocabulary, computed vs curated networks, the schema, data sources and licences, and building a co-expression network alongside the graph (§7) |
| [`part-2-building-a-knowledge-graph.ipynb`](part-2-building-a-knowledge-graph.ipynb) | 09:40–10:10 | Building the graph in NetworkX, degree and hubs, association scores, sparsification, annotation sparsity |
| [`part-3-querying-the-knowledge-graph.ipynb`](part-3-querying-the-knowledge-graph.ipynb) | 10:10–10:45 | **Practical.** ICD-10 coverage and ontology climbing (§A), shared-gene projection, separating cause from treatment, community detection (§B), bridging to the Session 2 omics data (§C, stretch), curated and computed networks together (§D, closing) |

Participant ("fill in the blanks") copies live in `participant/` and are generated
from the solution notebooks — never edit them by hand:

```bash
python ../../tools/make_participant_version.py "part-*.ipynb"
```

## Layout

```text
session-story.md the one-page narrative of the session
data/            generated data files (committed, small)
data-prep/       developer scripts - build_kg_data.py,
                 build_coexpression_data.py, make_images.py
images/          diagrams used by the notebooks
participant/     generated fill-in-the-blanks notebooks
s1_helpers.py    graph utilities, mostly adapted from co-expression material
concept-glossary.md
```

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
§7 builds a computed co-expression network from it, and Part 3 §D puts the two
networks together. 1.1 MB, so it needs no external download.

See [`data/README.md`](data/README.md) for provenance and licences.

Genes are keyed by Ensembl gene ID, matching the TCGA-BRCA matrix used in
Session 2, so the two join natively.

## Requirements

`pandas`, `numpy`, `networkx`, `matplotlib`, `seaborn`, `scipy` (see the
repository `requirements.txt`). Everything runs from the committed data except
Section C of Part 3, which needs the full Session 2 omics pickle; it is optional
and the notebook skips it cleanly if the file is absent.
