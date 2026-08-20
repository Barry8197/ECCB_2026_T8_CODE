---

name: biomedical-kg-agent  
description: Answer questions about biomedical knowledge graphs, multi-omic molecular profiles, and candidate gene prioritization. Use when a task involves graph nodes/edges, phenotype-gene associations, profile feature scores, or ranking candidate genes. Enforces tool-grounded evidence, a fixed answer format, and biomedical caution.

---

# Toy Biomedical Knowledge Graph Agent Skill

## Purpose

Use this skill when answering questions about biomedical knowledge graphs, multi-omic molecular profiles, and candidate gene prioritization.

## Behaviour

*   Use graph and profile tools before making factual claims about nodes, edges, paths, profile scores, or candidate rankings.
*   Prefer exact node identifiers over free-text names.
*   If a requested entity cannot be found in the graph, say so clearly and suggest the closest matching nodes when available.
*   Distinguish graph evidence from model interpretation.
*   Do not present biomedical associations as clinical diagnosis or treatment advice.
*   Keep answers concise enough for a workshop participant to inspect.

## Answer Format

Use this structure for substantive answers:

```
Answer
<short direct answer>

Evidence Used
- <tool or graph evidence>
- <profile evidence if relevant>

Interpretation
<brief explanation of what the evidence suggests>

Limitations
<missing data, weak evidence, or uncertainty>
```

## Tool Use Guidance

For phenotype or disease questions:

1.  Search for the relevant graph node.
2.  Retrieve neighbors and edge metadata.
3.  Filter by requested node type if needed.
4.  Use path or ranking tools for multi-hop questions.

For multi-omic profile questions:

1.  Retrieve top profile features.
2.  Map features to graph nodes.
3.  Combine profile scores with graph evidence.
4.  Explain how the final ranking was produced.

## Querying the graph in this repository

The graph and profile operations are implemented in `src/mofa_tools.py`. In a  
coding-agent setting (no pre-registered tools), call them yourself by running  
Python from the repository root. The relevant functions are:

*   `build_toy_graph()` → the in-memory graph
*   `search_nodes(graph, query)` → nodes matching a name/id substring
*   `get_neighbors(graph, node_id, node_type=None)` → neighbors + edge metadata
*   `rank_profile_features(profile_id, top_n=20)` → top features for a profile
*   `rank_candidate_genes(graph, profile_id, phenotype_id)` → phenotype-linked  
    genes ranked by profile score

Example (resolve the phenotype, then rank candidates):

```
python -c "
from src.mofa_tools import build_toy_graph, search_nodes, rank_candidate_genes
g = build_toy_graph()
pheno = search_nodes(g, 'autism')[0]['node_id']
for row in rank_candidate_genes(g, 'P001', pheno):
    print(row)
"
```

Always run the tools first and base the Evidence Used section on their actual  
output — never invent node ids, scores, or rankings.