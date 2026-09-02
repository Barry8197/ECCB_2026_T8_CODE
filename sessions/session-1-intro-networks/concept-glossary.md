# Glossary

# Session 1 : Introduction to Biological Networks and Knowledge Graphs

## Global concepts across all parts

| Term | Definition |
| --- | --- |
| **Network / graph** | A set of **nodes** (objects) and **edges** (pairs of objects that are related). Everything else is derived from those two sets. |
| **Node / vertex** | One object in the graph. Here: a gene, a disease, or an ICD-10 code. |
| **Edge** | A relationship between two nodes. |
| **Degree** | The number of edges attached to a node. A high-degree node is a **hub**. In a curated graph, degree measures *research attention* as much as biological importance. |
| **Path** | A sequence of edges leading from one node to another. Paths let a graph answer questions nobody stored the answer to. |
| **Connected component** | A group of nodes all reachable from one another. |
| **Density** | Fraction of possible edges that actually exist. Real biological graphs are sparse (ours: 0.0046). |
| **Clustering coefficient** | How often a node's neighbours are themselves connected. Near zero here, because the graph is close to bipartite. |
| **Directed / undirected** | Whether edges have a direction. `is_a` is genuinely directed; we work undirected but keep the direction as an edge attribute. |
| **Self-loop** | An edge from a node to itself. None in this graph. |
| --- | --- |
| **Knowledge graph** | A graph whose nodes and edges are **typed** and carry **provenance** — we can tell what kind of entity each node is, what kind of relationship each edge is, and on what evidence. |
| **Computed network** | A network whose edges are *calculated from measurements* — correlate every pair, keep what clears a threshold. Contrast with curated. |
| **Curator / curation** | A biologist employed by a database to read published papers and record their findings as structured entries (gene X, disease Y, evidence PMID Z). Also covers automated contributors — text mining, GWAS pipelines. Every edge in a knowledge graph exists because a curator or pipeline put it there. |
| **Curated network** | A network whose edges are *asserted facts* read out of a database. Can only contain what somebody already recorded. |
| **Node type** | The kind of entity a node represents: `gene`, `disease`, or `icd10`. |
| **Edge type** | The kind of relationship: `associated_with`, `is_a`, or `maps_to`. |
| **Provenance** | The evidence attached to an edge — here an association score and an evidence count. |
| **Bipartite graph** | A graph whose nodes split into two groups with edges only *between* groups, never within. Genes and diseases here are almost perfectly bipartite. |
| **Projection** | Collapsing a bipartite graph onto one node type, joining two nodes when they share neighbours. Turns "which diseases resemble each other?" into a structural question. |
| **Jaccard index** | Size of the intersection divided by size of the union. Used to score shared genes so that well-studied diseases do not automatically rank highest. |
| **Sparsification** | Removing edges to make a graph readable or to keep only strong relationships — by weight threshold, top percentage, or degree. |

## Part 1: What is a Network, and What is a Knowledge Graph?

### Data sources

| Term | Definition |
| --- | --- |
| **[Open Targets Platform](https://platform.opentargets.org/)** | Drug-target discovery platform aggregating gene–disease evidence from genetics, somatic mutations, drugs, pathways, expression and literature. Licensed [CC0 1.0](https://platform-docs.opentargets.org/licence) (public domain). Source of our genes, diseases, hierarchy and association scores. |
| **[MONDO Disease Ontology](https://mondo.monarchinitiative.org/)** | Unified disease ontology merging several older vocabularies, and the source of our MONDO → ICD-10 cross-references. Licensed [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). |
| **SSSOM** | Simple Standard for Sharing Ontological Mappings — the file format MONDO publishes its cross-references in (`mondo.sssom.tsv`). |
| **DisGeNET** | A widely cited gene–disease database. Commercially licensed since 2024, so it cannot be redistributed in a public teaching repository. |
| **Ensembl gene ID** | Stable gene identifier of the form `ENSG00000012048`. The identifier space shared by Open Targets and the Session 2 TCGA-BRCA matrix, which is what lets the two join. |
| **TCGA-BRCA** | The Cancer Genome Atlas breast cancer cohort. Source of the expression matrix in `coexpr_expression.csv.gz` (500 patients × 737 genes) and of the multi-omics data in Session 2. |

### Co-expression (Part 1 §7)

| Term | Definition |
| --- | --- |
| **Co-expression** | Two genes whose expression rises and falls together across patients. Evidence that they are *regulated* together — which may mean they work together, or merely that both track a third thing like proliferation. |
| **Correlation matrix** | Every gene correlated against every other. 737 genes give 271,216 pairs. The diagonal (each gene with itself) is zeroed, otherwise it dominates any ranking. |
| **Correlation threshold** | The cut-off above which a correlation becomes an edge. **There is no principled value.** Sweeping 0.3 → 0.8 takes the same data from 39,593 edges to 200, and every one of those networks is defensible. |
| **Proliferation confound** | The dominant structure in any tumour co-expression network: cell-cycle genes correlate with each other and with anything expressed in dividing cells. Responsible for `TOP2A`, `BUB1B` and `KNL1` appearing next to `BRCA1`. Not a discovery. |
| **Corroboration** | The same relationship supported by two sources whose errors are *independent* — e.g. `BRCA2`–`BRIP1` both curated as causal and co-expressed. Strong. Contrast with the tubulins, where both sources fail for related reasons. |
| **Candidate** | A gene co-expressed with a disease's causal genes that carries no curated edge to that disease. Where gene–disease discovery lives — and where the proliferation confound also lands, indistinguishably. |
| **Association score** | Open Targets' 0–1 summary of gene–disease evidence. **Not a probability** — a weighted aggregation across evidence types, so a high score can reflect genetics *or* drug treatment. |
| **Evidence count** | How many individual pieces of evidence support an association. |
| **Evidence type (datatype)** | The *kind* of claim behind an association: `genetic_association`, `somatic_mutation`, `known_drug`, `literature`, `rna_expression`, `animal_model`, `affected_pathway`, `genetic_literature`. Kept in `kg_evidence.csv`. |
| **Causal evidence** | `genetic_association` and `somatic_mutation` — the types asserting the gene helps *cause* the disease, as opposed to being a place to *treat* it. |
| **`known_drug` evidence** | An association arising because a drug hitting this target treats this disease. A statement about therapy, not aetiology — and the reason tubulins look like breast cancer genes. |

## Part 2: Building a Knowledge Graph with NetworkX

| Term | Definition |
| --- | --- |
| **NetworkX** | The Python graph library used throughout. `nx.Graph()` is undirected, `nx.DiGraph()` directed. |
| **Node attribute** | Arbitrary data stored on a node — here `type`, `name` and `extra`. |
| **Edge attribute** | Arbitrary data stored on an edge — here `type`, and `weight`/`evidence` on association edges only. |
| **Unweighted edge** | An edge with no `weight` attribute. `is_a` and `maps_to` are unweighted; treating a missing weight as zero silently deletes them when thresholding. |
| **Construction artefact** | A property of the graph caused by how it was built rather than by the underlying biology — e.g. diseases are hubs here only because the build kept the top 30 genes *per disease*. |
| **Annotation sparsity** | Missing edges that should exist, because nobody has studied the relationship or because the evidence was recorded against a different term. The default state of biomedical knowledge graphs. |
| **Vocabulary mismatch** | When a fact exists but is filed under an equivalent-but-different ontology term — e.g. basal-like breast carcinoma has zero genes here because the evidence sits on triple-negative breast carcinoma. |

## Part 3: Querying the Knowledge Graph

| Term | Definition |
| --- | --- |
| **[ICD-10](https://icd.who.int/browse10/2019/en)** | WHO International Classification of Diseases, 10th revision — the clinical coding system used in hospital records, death certificates and insurance claims. |
| **ICD10CM** | The US "Clinical Modification" of ICD-10, more finely subdivided than the WHO version. |
| **ICD10WHO** | The WHO's international version of ICD-10. |
| **[ICD-O-3](https://www.who.int/standards/classifications/other-classifications/international-classification-of-diseases-for-oncology)** | International Classification of Diseases *for Oncology* — the separate classification that codes tumour **morphology**. This is where subtype information lives; ICD-10 does not carry it. |
| **Granularity mismatch** | Two vocabularies describing the same thing at different resolutions. ICD-10 subdivides breast cancer *anatomically* (by quadrant), so it has no code for a molecular subtype — and no data source can supply one. |
| **Ontology** | A structured vocabulary with defined relationships between terms. The `is_a` edges make ours ontology-backed rather than flat. An ontology is a **network, not a tree** — a term can have several parents (triple-negative breast carcinoma has three, one per receptor). |
| **`is_a` / subsumption** | "X is a kind of Y". Directional and transitive, which is what makes climbing possible. |
| **Ancestor** | Any term reachable by following `is_a` edges upwards. |
| **Ontology traversal / climbing** | Following `is_a` edges to recover information missing on a specific term — raises ICD-10 coverage here from 19/90 to 82/90. |
| **Inherited mapping** | A code found on an ancestor rather than the node itself. Real, but a weaker claim: it describes the parent, and the mapping is many-to-one. |
| **Cross-reference (xref)** | A recorded correspondence between a term in one vocabulary and a term in another. |
| **PAM50** | A 50-gene expression signature classifying breast tumours into LumA, LumB, Basal, Her2 and Normal. The Session 2 prediction target — and all five labels are disease nodes in this graph. |
| **Versioned Ensembl ID** | An Ensembl ID carrying an annotation-release suffix (`ENSG00000012048.23`). TCGA uses these, Open Targets does not; joining without stripping the suffix matches **nothing**, silently. |
| **HBOC** | Hereditary Breast and Ovarian Cancer syndrome — driven by BRCA1/BRCA2/BRIP1, and recoverable from this graph purely by asking which genes breast and ovarian cancer share. |
| **Community detection** | Partitioning a graph into groups of nodes more connected to each other than to the rest. Louvain on our causal projection recovers cancers, autoimmune and metabolic diseases without being told they exist. |
| **Modularity** | What community detection optimises: how much denser the within-group connections are than chance would give. |
| **Promiscuous gene** | A gene associated with many unrelated diseases (TP53, tubulins). High degree, low specificity — a shared promiscuous gene rarely means shared biology. |
