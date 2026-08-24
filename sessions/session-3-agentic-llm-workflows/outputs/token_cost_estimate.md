# Token & Cost Estimate — Notebooks 01 & 02

**Model:** `claude-haiku-4-5` — $1.00 / 1M input tokens, $5.00 / 1M output tokens
(cache read ~$0.10/1M, cache write ~$1.25/1M; no caching used, so both are 0).

Token counts are exact, taken from the Anthropic API `usage` on each response.
Each agent query makes one or more API calls (one per tool-call round-trip).

**Measured: 2026-08-17 (fresh clean pass).** Raw evidence:
[usage_0102_haiku_fresh.json](usage_0102_haiku_fresh.json). An earlier pass on
2026-07-27 produced the same call counts and token totals to within ~1%.

| Query | API calls | Input tok | Output tok | Cache tok | Cost |
|---|--:|--:|--:|--:|--:|
| 01-Q1 data shape | 2 | 2,399 | 291 | 0 | $0.00385 |
| 01-Q2 split | 2 | 2,350 | 194 | 0 | $0.00332 |
| 01-Q3 active factors | 2 | 2,426 | 327 | 0 | $0.00406 |
| 01-Q4 factor vs view | 2 | 2,376 | 207 | 0 | $0.00341 |
| 01-Q5 flagship drivers | 3 | 4,832 | 620 | 0 | $0.00793 |
| 01-Q6 predict + confusion | 2 | 2,373 | 206 | 0 | $0.00340 |
| 01-Q7 variance vs separation | 2 | 2,746 | 649 | 0 | $0.00599 |
| 01-Q8 out-of-sample | 2 | 2,649 | 675 | 0 | $0.00602 |
| 02-MCP flagship (Q5) | 3 | 5,353 | 570 | 0 | $0.00820 |
| **Total (01 + 02)** | **20** | **27,504** | **3,739** | **0** | **$0.0462** |

## Notes

- **Cache tokens = 0**: the notebooks don't use prompt caching, and each query is
  an independent conversation, so there is nothing to reuse. A shared cached
  system+tools prefix would lower input cost on repeated runs.
- Multi-step queries (Q5, Q7, Q8, and the MCP flagship) cost more because the
  agent makes an extra tool-call round-trip and writes a longer answer.
- Notebook 03 (Pi coding agent) is covered in
  [cost_estimation_notes.md](cost_estimation_notes.md) — its usage is parsed
  from Pi's session logs (evidence:
  [usage_03_pi_haiku_fresh.json](usage_03_pi_haiku_fresh.json)); measured
  $0.17–0.32 per pass (run length is nondeterministic).
