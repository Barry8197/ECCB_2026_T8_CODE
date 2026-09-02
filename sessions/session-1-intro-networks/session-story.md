# The story of Session 1

*A one-page narrative of what we do and why. For the detail, see the notebooks.*

## What we want to understand

We want to understand what a biological network is, what a **knowledge graph** is,
and — the part people usually skip — how the two differ. A co-expression network is
*computed* from measurements: correlate every gene against every other gene and keep
what clears a threshold. A knowledge graph is *read out of recorded facts*: somebody
looked at the evidence and wrote down that this gene is implicated in that disease.
They look identical when drawn and they fail in opposite ways, so telling them apart
is the first thing worth learning.

So we build **both**, over the same genes, and spend the rest of the session on
what each one can and cannot say.

## How we do it

We build a real knowledge graph rather than a toy one. We download two public,
redistributable sources — [Open Targets](https://platform.opentargets.org/) (CC0)
for genes, diseases and gene–disease association scores, and the
[MONDO Disease Ontology](https://mondo.monarchinitiative.org/) (CC BY 4.0) for the
disease hierarchy and for cross-references out to clinical ICD-10 codes — and cut
them down to a breast-cancer-focused subset small enough to read by eye:
**881 nodes** (760 genes, 90 diseases, 31 ICD-10 codes) and **1,768 edges**.

Alongside it we ship the TCGA-BRCA transcriptomics matrix — the same one Session 2
works from — cut down to the 737 graph genes that are actually measured. From that
we compute a co-expression network: correlate all 737 genes, keep the pairs above
a threshold. Two networks, same genes, same identifiers, opposite provenance.

Around breast cancer and its subtypes we add **24 comparison diseases** — other
cancers, autoimmune and inflammatory conditions, neurological and cardiometabolic
disease — so that "which diseases resemble each other?" has a real answer rather
than one big breast-cancer blob.

We then build the graph in NetworkX, with **typed nodes** (gene, disease, ICD-10),
**typed edges** (`associated_with`, `is_a`, `maps_to`) and, crucially, **provenance**:
a separate table splitting every gene–disease edge by the *kind* of evidence
behind it.

## What we connect to what

Two joins do the work.

We connect **clinical vocabulary to molecular vocabulary**, by walking from an
ICD-10 code that a hospital would actually bill, through the MONDO disease
hierarchy, down to genes.

And we connect **the knowledge graph to the omics data** from Session 2. This is why
genes are keyed by Ensembl ID: it is the same identifier space as the TCGA-BRCA
transcriptomics matrix, so a gene list coming out of the multi-omics work can be
looked up in the graph directly.

## What we get out of it

**The same gene looks completely different in the two networks.** `BRCA1`
co-expresses with the proliferation programme — tubulin-adjacent mitotic genes
that have nothing to do with it — plus `PALB2` and `FANCD2`, which are genuine
repair partners the knowledge graph has no way to record. Meanwhile the curated
graph gives its hereditary breast–ovarian diseases. `TP53` is starker still: a
27-disease hub in the curated graph, and essentially invisible to co-expression
(strongest |r| = 0.18), because it is regulated by protein stability rather than
transcription. Neither absence means *unimportant*; they mean different things,
and knowing which is the skill.

**And the computed network has a knob the curated one does not.** Sweeping the
correlation threshold from 0.3 to 0.8 takes the same data from 39,593 edges to
200. Nothing in the data says where to stop. The knowledge graph has exactly
1,768 edges whatever we think — which removes that problem and introduces the
opposite one: we cannot loosen the criteria to see more.

**Breast and ovarian cancer share BRCA1, BRCA2 and BRIP1.** We never told the graph
these two diseases were related — we asked which genes they had in common and
hereditary breast-ovarian cancer syndrome fell out of the structure.

**But they also "share" tubulins.** TUBB, TUBA1B and the topoisomerases show up in
exactly the same list, not because of shared biology but because both cancers are
treated with taxanes and anthracyclines. The overall association score cannot tell
the two apart. Splitting by evidence type can: filter to genetic and somatic
evidence and every tubulin disappears, leaving the real hereditary genes. The
advantage of carrying provenance is that we can ask *what kind of claim is this*,
not just *how strong is it*.

**The diseases group themselves.** Running community detection on shared causal
genes alone recovers the clinical taxonomy — a cancer community, an
autoimmune/inflammatory community, a metabolic one — without ever being told those
categories exist. Structure nobody encoded turns out to be recoverable from
structure somebody did.

**And the gaps teach as much as the findings.** Only 19 of 90 diseases carry an
ICD-10 code, because ICD-10 subdivides breast cancer *anatomically* and has no code
for "triple-negative" — climbing the ontology recovers 82 of 90. Basal-like breast
carcinoma has zero genes, because the evidence was filed under a near-synonym.
Joining the omics matrix without stripping the Ensembl version suffix matches
*nothing at all*.

None of those are bugs. They are what working with real biomedical knowledge is
like, and recognising them is the skill.

## Putting the two networks back together

The session closes on the question its title actually asks: *can molecular data
strengthen the evidence for a gene–disease relationship?*

Take the 16 genes with causal evidence for breast cancer and ask what co-expresses
with them. Three things come back.

**Corroboration.** Ten of 120 causal gene pairs co-express — and the strongest are
`BRCA2`–`BRIP1`, `BARD1`–`BRIP1`, `BARD1`–`BRCA2`: the homologous-recombination
complex, recovered twice over. The graph knows because families were sequenced;
the matrix knows because they are transcribed together. Two independent routes,
one biology.

**A confound we had already removed.** The seven strongest partners carrying any
curated breast-cancer edge are `TOP2A`, `TYMS`, `TUBA1B`, `TUBA1C`, `TUBB`, `TOP1`
and `CDK6` — precisely the chemotherapy targets we filtered out an hour earlier.
They come back because they are proliferation genes, not because the filtering was
wrong. Agreement between two sources is only evidence when their errors are
independent, and here they are not.

**And 184 candidates nothing can sort.** Genes that co-express with the causal set
and carry no curated breast-cancer edge at all. `BLM`, `EXO1`, `FANCD2`, `RAD51`,
`MSH6` — real genome-stability genes worth following up — sit interleaved with
`BUB1B`, `HMMR`, `KNL1` and `RRM2`, which are proliferation artefacts. No column
in the table separates them; that judgement needs biology.

## Where it goes next

Session 3 builds LLM agents — tools, MCP and skills — and Session 4 turns them on
graphs like this one, queried from multi-omics profiles. Everything above — the
coverage gaps, the inherited codes, the technically-correct-but-misleading edges,
and that final unsortable list — is what those agents have to get right, and what
we need in order to check them.
