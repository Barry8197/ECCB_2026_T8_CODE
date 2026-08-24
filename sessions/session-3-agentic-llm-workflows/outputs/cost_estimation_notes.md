# Session 3 — API Cost Estimation Notes

## Summary

| Item | Value |
|---|---|
| **Model** | [`claude-haiku-4-5`](https://www.anthropic.com/claude/haiku) |
| **Per person, one clean pass** (measured 2026-08-17) | **$0.22** |
| **Per person, realistic** (room for 2–3× cell re-runs + 03 run-length variance) | **~$1.00** |
| **50-person audience** | **~$50 expected; set a $100 workspace spend cap** |

Details, evidence, and reproduction instructions below.

---

Per-participant cost estimate for running the Session 3 notebooks
(00 data → 01 tools → 02 MCP → 03 skills/Pi), plus how the numbers were
measured so a fresh session can reproduce or extend them.

## Model used for this estimate

**All numbers below are for `claude-haiku-4-5`** ($1.00/1M input, $5.00/1M
output, cache read $0.10/1M, cache write $1.25/1M). The notebooks set this
model in 01/02 (`MODEL = "claude-haiku-4-5"`) and 03
(`MODEL = "anthropic/claude-haiku-4-5"` for Pi).

Other model runs (e.g. Sonnet) to be added later for comparison — repeat the
measurement below with the model string swapped and pricing adjusted.

### Is Haiku good enough? (qualitative review of the agent outputs)

Yes, good enough for this workshop — with one caveat.

What Haiku did well, which is most of what the notebooks need:

- **Tool selection was flawless.** Every query in 01 picked the right tool(s)
  first try, no flailing — Q5 correctly chained `factor_subtype_association` →
  `top_features_for_factor`, Q7 pulled both rankings and compared them. That's
  the pedagogical core (watch an agent choose and chain tools), and it works.
- **Grounding held.** Numbers in the answers matched tool outputs; nothing
  invented. The system prompt's "cite which tool" discipline was followed.
- **03/Pi was genuinely solid** — it explored the repo, loaded the cached model
  (didn't try to re-fit), followed the skill's
  Answer/Evidence/Interpretation/Limitations format, generated the plots.
  25 tool calls of competent coding-agent behavior on Haiku is honestly
  impressive per dollar.

The caveat: **interpretation depth is where you feel the model tier.** The
biology narrative is correct but thin — e.g. it reports "Factor2 has η²=0.74
and separates subtypes" but won't spontaneously connect the top transcriptomic
drivers to known basal/luminal biology, or discuss why LumA↔LumB confusion is
biologically expected. Q7 (the "explain the mismatch" reasoning query) got the
right facts but a fairly mechanical explanation.

**Recommendation:** keep Haiku as the workshop default (fast, cheap, reliable
mechanics). Optionally pre-run the flagship Q5 / 03 interpretation once on a
bigger model (Sonnet) and include it as a "same query, bigger model"
comparison cell — a teaching moment paid for once, not ×50.

## Measured cost — clean passes

Two full clean passes were run, three weeks apart:

| Part | 2026-07-27 pass | 2026-08-17 pass (locked) |
|---|--:|--:|
| 00 data exploration (no LLM) | $0 | $0 |
| 01 — all 8 agent queries | $0.038 (17 calls) | $0.038 (17 calls) |
| 02 — MCP flagship query | $0.008 (3 calls) | $0.008 (3 calls) |
| 03 — Pi coding agent, 3 runs (Q9, Q10, no-skill) | $0.315 (48 calls) | $0.172 (27 calls) |
| **Total, one clean pass** | **$0.36** | **$0.22** |

2026-08-17 pass detail (Haiku 4.5):

| Part | Calls | Input | Output | CacheR | CacheW | Cost |
|---|--:|--:|--:|--:|--:|--:|
| 01 + 02 | 20 | 27,504 | 3,739 | 0 | 0 | $0.046 |
| 03 (Pi, 3 runs) | 27 | 14,376 | 15,426 | 221,826 | 46,416 | $0.172 |

Per-query breakdown for 01/02 is in `token_cost_estimate.md`.

Observations:

- **01/02 are near-deterministic**: token counts matched across the two passes
  to within ~1% (same tool choices, same call counts).
- **03 (Pi) varies run to run** — the coding agent's number of turns is
  nondeterministic (48 calls / $0.315 in July vs 27 calls / $0.172 in August
  for the same three prompts). Budget for the *upper* observed value.
- 03 dominates cost either way: large outputs + heavy prompt-cache traffic
  (Pi caches aggressively; cache reads are cheap but hundreds of thousands of
  them add up).
- 01/02 use no prompt caching (independent short conversations).
- **Working figure for one clean pass: $0.25–0.40.**

## Evidence / provenance for the numbers above

**2026-08-17 pass (the locked numbers) — raw evidence in this repo:**

- `outputs/usage_0102_haiku_fresh.json` — per-query token counts + cost for
  01/02, dumped directly from each response's `response_metadata["usage"]`
  during the run.
- `outputs/usage_03_pi_haiku_fresh.json` — per-run token counts + cost for
  the three Pi runs, summed from Pi's own session logs.
- Pi's raw session logs (the ultimate source for 03) are at:
  `~/.pi/agent/sessions/--Users-chaeeunlee-Documents-VSC_workspaces-ECCB2026_TEST-sessions-session-3-agentic-llm-workflows--/2026-08-17T*.jsonl`
  — sum `message.usage.cost.total` over assistant messages to recompute.

**2026-07-27 pass:**

- 03: same Pi log directory, the `2026-07-27T09-*.jsonl` files
  (probe $0.0027, Q9 $0.0776, Q10 $0.0872, no-skill $0.1496).
- 01/02: transcribed from live API usage during the run; the raw dump went to
  a temp scratchpad that was cleaned up, so no on-disk file for that pass —
  superseded by the 2026-08-17 in-repo evidence anyway (the two passes agree
  to ~1% on 01/02).

**The executed notebooks** (`notebooks/01_tools_v1.ipynb`, `02_MCP_v1.ipynb`,
`03_agent_skills_v1.ipynb`) contain the saved agent outputs (03's cells are
from the 2026-08-17 run) — evidence the runs happened, though cell outputs
don't include token counts.

## Per-participant estimate (50-person audience)

| Scenario | Per person | 50 people |
|---|--:|--:|
| One clean pass (upper of two measured passes) | $0.36 | $18 |
| Realistic (2–3× re-runs of agent cells, esp. 03) | $0.75–1.10 | $38–55 |
| Safe cap | $2.00 | $100 |

Main source of deviation: participants re-running cells.

**Recommendation:** dedicated Console *workspace* for the workshop with a
~$100 spend cap; one shared key (or a few keys to split blast radius);
revoke after the session.

## How the numbers were measured (reproduce in a fresh session)

### 01 / 02 (LangChain)

Each `AIMessage` returned by `ChatAnthropic` carries exact usage in
`ai.response_metadata["usage"]` — keys `input_tokens`, `output_tokens`,
`cache_read_input_tokens`, `cache_creation_input_tokens`. Accumulate across
every call of the agent loop (one API call per tool round-trip) and multiply
by pricing.

### 03 (Pi coding agent) — where to find usage

Pi does **not** print usage on stdout, but it logs everything per session:

```
~/.pi/agent/sessions/<cwd-slug>/<timestamp>_<id>.jsonl
```

- `<cwd-slug>` is the working directory Pi ran in, with `/` → `-`
  (e.g. `--Users-...-ECCB2026_TEST-sessions-session-3-agentic-llm-workflows--`).
- One file per `pi -p ...` invocation; match files to runs by mtime.
- Each line is a JSON record. Records with `type == "message"` and
  `message.role == "assistant"` carry a `message.usage` dict:

```json
"usage": {
  "input": 2155, "output": 166,
  "cacheRead": 0, "cacheWrite": 0,
  "totalTokens": 2321,
  "cost": {"input": 0.002155, "output": 0.00083, "total": 0.002985}
}
```

⚠️ Note the key names: `input` / `output` / `cacheRead` / `cacheWrite` —
**not** the Anthropic-style `input_tokens` etc. (a previous parse failed for
exactly this reason). Pi also self-computes `cost.total` at the model's
prices — summing `cost.total` over all assistant messages in the session
files gives the run cost directly.

### Optional Console cross-check

Create a project-named API key (e.g. `eccb-session3-estimate`), put it in
`.env`, run one clean pass, then read that key's usage line in the Anthropic
Console — an independent billing-side verification (~$0.40 on Haiku).
