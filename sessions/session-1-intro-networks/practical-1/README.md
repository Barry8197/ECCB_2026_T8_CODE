# Practical 1: Querying a Phenotype Knowledge Graph

The practical is delivered as
[`../part-3-querying-the-knowledge-graph.ipynb`](../part-3-querying-the-knowledge-graph.ipynb).
Work through the participant copy in
[`../participant/`](../participant/) and fill in the `### YOUR CODE HERE ###` blocks.

## Objective

Start from ICD-10 codes — the vocabulary a hospital actually uses — and get from
there to genes. Discover why that only partly works, recover what we can by
traversing the disease ontology, and then ask the question the graph is genuinely
good at: which diseases share genes?

## What we will do

- [ ] Measure ICD-10 coverage across the graph's 71 diseases (spoiler: 7 of them)
- [ ] Diagnose *why* — a granularity mismatch, not a missing data source
- [ ] Climb the `is_a` hierarchy to raise coverage from 7 to 70
- [ ] Project the bipartite gene–disease graph onto diseases alone
- [ ] Recover hereditary breast-ovarian cancer syndrome from shared genes
- [ ] Spot the shared genes that are pharmacological rather than biological
- [ ] Join the graph to the Session 2 transcriptomics matrix without silently
      matching nothing

## Stretch goals

- [ ] Compare Jaccard against raw shared-gene counts — which diseases change rank, and why?
- [ ] Find the disease pairs whose similarity is carried entirely by promiscuous genes
- [ ] Export a subgraph to GraphML and open it in [Cytoscape](https://cytoscape.org/)
