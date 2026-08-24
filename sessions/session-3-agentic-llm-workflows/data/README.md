## Toy Data Shape for Session 3

The LLM agent practical needs three linked objects:

1.  A knowledge graph with typed nodes and edges.
2.  Multi-omic profile scores from Session 2.
3.  A mapping from molecular features to graph nodes.

## Minimal Toy Schema

Graph nodes:

```
node_id,name,type
HP:0000729,Autism spectrum disorder,phenotype
GENE:SHANK3,SHANK3,gene
GENE:CHD8,CHD8,gene
```

Graph edges:

```
source,target,relation,evidence
GENE:SHANK3,HP:0000729,associated_with,toy literature edge
GENE:CHD8,HP:0000729,associated_with,toy literature edge
```

Profile scores:

```
profile_id,feature_id,score
P001,GENE:SHANK3,2.4
P001,GENE:SCN2A,1.7
```