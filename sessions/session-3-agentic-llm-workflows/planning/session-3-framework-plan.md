# Session 3 — Framework Plan & TODO

The three notebook tracks all answer the **same task** against the same  
`src/mofa_tools.py` backend:

> For profile `P001`, which genes connected to autism spectrum disorder  
> should we prioritize, and why? → `SHANK3, SCN2A, CHD8`

---

## Layer 1 — Tool calling (notebook 01)

How the LLM is given Python functions and decides when to call them.

**Agent frameworks** 

| Option | Notes | Link |
| --- | --- | --- |
| **LangChain / LangGraph** | De-facto standard; `@tool` + `bind_tools`. Most participants will have seen it. LangGraph adds explicit graph control flow. | https://www.langchain.com · https://langchain-ai.github.io/langgraph/ |
| LlamaIndex (Workflows) | Strong for RAG/retrieval-centric agents; heavier than we need here. | https://www.llamaindex.ai |
| Pydantic AI | Type-safe, clean API, lighter than LangChain; smaller ecosystem. | https://ai.pydantic.dev |
| Smolagents (Hugging Face) | Minimal, code-writing agents; very lightweight. | https://github.com/huggingface/smolagents |
| CrewAI | Role/crew-oriented multi-agent orchestration; more than we need for one task. | https://www.crewai.com |
| AG2 (AutoGen fork) | Conversational multi-agent framework. | https://github.com/ag2ai/ag2 |
| Microsoft Agent Framework | MS's unified agent stack (Semantic Kernel + AutoGen lineage). | https://github.com/microsoft/agent-framework |
| Haystack | deepset's pipeline/agent framework, retrieval-strong. | https://haystack.deepset.ai |
| Letta / MemGPT | Memory-centric stateful agents. | https://www.letta.com |
| DSPy | Programmatic prompt/program optimization rather than a tool-loop per se. | https://dspy.ai |
| Mastra (TS) | TypeScript agent framework — only if we go JS. | https://mastra.ai |
| Vercel AI SDK (TS) | TypeScript tool-calling SDK — only if we go JS. | https://sdk.vercel.ai |

**Vendor / lab SDKs**

| Option | Notes | Link |
| --- | --- | --- |
| Raw provider SDKs (no framework) | Call tool-use directly — maximum transparency, most boilerplate. | [Anthropic](https://docs.claude.com/en/api) · [OpenAI](https://platform.openai.com/docs) |
| Claude Agent SDK | Anthropic's agent loop; also does Skills (see Layer 3). | https://docs.claude.com/en/api/agent-sdk |
| OpenAI Agents SDK | Tidy agent loop; OpenAI-leaning, less neutral for teaching. | https://openai.github.io/openai-agents-python/ |
| Google ADK | Google's Agent Development Kit; Google-leaning. | https://google.github.io/adk-docs/ |

**➡️ Tentative choice: LangChain.** Minimal surprise for participants, `@tool`  
wraps our existing `mofa_tools.py` functions, and  
it keeps notebook 01 the simplest possible "LLM calls a Python function" demo.

---

## Layer 2 — MCP server + client (notebook 02)

How the same functions get exposed over the Model Context Protocol, then  
pulled back into the agent.

**Options  — server side**

| Option | Notes | Link |
| --- | --- | --- |
| **FastMCP** | Pythonic, decorator-per-tool; least plumbing. We already have a scaffold under `server/`. | https://github.com/PrefectHQ/fastmcp · https://gofastmcp.com |
| Official MCP Python SDK (low-level) | More control, more boilerplate; FastMCP is the high-level layer on top of it. | https://github.com/modelcontextprotocol/python-sdk |

**Options — client side**

| Option | Notes | Link |
| --- | --- | --- |
| **langchain-mcp-adapters** | Loads MCP tools back into the _same_ LangChain agent shape as notebook 01. | https://github.com/langchain-ai/langchain-mcp-adapters |
| Raw MCP client SDK | Direct protocol calls; loses the "agent code barely changes" payoff. | https://github.com/modelcontextprotocol/python-sdk |

**➡️ Tentative choice: FastMCP server + langchain-mcp-adapters client.**  
FastMCP exposes the same `mofa_tools.py` functions with a decorator each;  
langchain-mcp-adapters then loads them into the identical agent shape from  
notebook 01. 

---

## Layer 3 — Coding Agent Skills (notebook 03)

How reusable skill instructions (`skills/biomedical-kg-agent.SKILL.md`) get  
loaded and change model behaviour. The skill file itself is harness-agnostic. (Agent Skill Standard: https://agentskills.io)

There are two sub-families here: full **coding-agent harnesses/CLIs** (load  
the `SKILL.md` and run a real agent loop) and **library/DIY** approaches  
(inject the skill yourself). The skill file is the same either way.

**Coding-agent harnesses / CLIs**

| Option | Notes | Link |
| --- | --- | --- |
| **Pi** | Minimal harness (Mario Zechner); first-class Skills, very token-efficient, lightweight. Our tentative pick. | https://github.com/earendil-works/pi · https://pi.dev |
| Claude Code | Anthropic's CLI; native Agent Skills (`SKILL.md`) support, well documented. Heavier, opinionated. | https://github.com/anthropics/claude-code |
| OpenHands (ex-OpenDevin) | Full software-agent platform with a public Skills registry; powerful but heavy for a notebook session. | https://github.com/OpenHands/OpenHands |
| OpenCode | Open-source terminal agent, multi-provider TUI. Note: original `opencode-ai` repo archived Sep 2025; active at opencode.ai. | https://opencode.ai |
| OpenAI Codex CLI | OpenAI's terminal coding agent; OpenAI-leaning. | https://github.com/openai/codex |
| Gemini CLI | Google's terminal agent; Google-leaning. | https://github.com/google-gemini/gemini-cli |
| Goose | Block's extensible local agent (MCP-native extensions). | https://github.com/block/goose |
| Aider | Mature pair-programming CLI; git-centric, no formal "skills" concept. | https://github.com/Aider-AI/aider |

**Library / DIY (no separate CLI)**

| Option | Notes | Link |
| --- | --- | --- |
| Claude Agent SDK | Native Skills support as a Python library; more control than a CLI. | https://github.com/anthropics/claude-agent-sdk-python |
| LangChain + skill-as-system-prompt | Reuse Layer 1; inject skill text into the system message (already sketched in nb03). Lowest new-dependency cost, most notebook-native. | https://github.com/langchain-ai/langchain |

**➡️ Tentative choice: Pi.** Of the harnesses it's the lightest, Skills are a  
first-class concept, and it's not tied to one model provider. 

---

## Decision summary

| Layer | Choice | Status |
| --- | --- | --- |
| Tool calling (NB 01) | **LangChain** | ✅ tentative |
| MCP (NB 02) | **FastMCP** server + **langchain-mcp-adapters** client | ✅ tentative |
| Skills (NB 03) | **Pi** | ✅ tentative |

---

## General plan / build order

1.  **Finalise** `**mofa_tools.py**` as the single shared backend  
    (confirm tool set with Session 1 folks — current tools are dummies).
2.  **Notebook 01 — LangChain tool calling**  
    \- Wrap `search_nodes`, `get_neighbors`, `rank_profile_features`,  
    `rank_candidate_genes` with `@tool`.  
    \- Bind to a Claude chat model, run the agent loop on the shared task.  
    \- Keep the existing "call functions ourselves" cells as the no-LLM  
    fallback.
3.  **Notebook 02 — MCP**  
    \- Stand up FastMCP server in `server/` exposing the same functions.  
    \- Load via langchain-mcp-adapters into the same agent shape as 01.  
    \- Show the ranking matches notebook 01 (trace stub already there).
4.  **Notebook 03 — Skills (Pi)**  
    \- Load `biomedical-kg-agent.SKILL.md` into the Pi harness.  
    \- Demo without-skill vs with-skill behaviour (evidence discipline).
5.  **Fallback mode** — notebooks should work without API keys via direct  
    calls / replayed traces (already required in `notebooks/README.md`).
6.  **Slides** — finalise concept deck (tools → MCP → skills).

## Still-open decisions (not framework-related)

*   Confirm final tool set with Session 1 & 2 groups.
*   Model access options for participants (shared key vs bring-your-own vs  
    local vLLM instance).
*   Whether participants run MCP locally or only watch the instructor demo.
*   Confirm Session 1 graph + Session 2 MOFA outputs can be reused directly.