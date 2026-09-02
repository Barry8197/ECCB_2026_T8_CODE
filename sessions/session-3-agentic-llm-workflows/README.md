# ECCB 2026 Tutorial T08 Session3: Building Agentic LLM Workflows for Biomedical Knowledge Graphs

## Session Shape

Brief theoretical reviews on concepts -> hands-on practice session

## Repository Structure

```
ECCB2026/
  assets/       Internal assets
  data/         As of now just mock test data
  notebooks/    Jupyter notebooks for guided practice session
  planning/     Schedules, proposal notes, task lists
  server/       MCP server or tool-hosting demo code
  skills/       Example Skill.md-style reusable instructions
  slides/       Slide outlines and generated slide decks
  src/          Shared Python utilities
```

## MCP Server Setup

The repository includes a starter MCP server for Session 3:

*   [Toy MCP server](server/toy_kg_mcp_server.py)
*   [Underlying graph tools](src/mofa_tools.py)

## Local Setup

Install the starter Python dependencies with:

```
python3 -m pip install -r requirements.txt
```

The first scaffold includes a tiny toy graph in `src/mofa_tools.py` so the LLM-agent workflow can be prototyped before the final Session 1 and Session 2 outputs are ready.

## Session 3 Comparison Materials

The agentic workflow materials are split into three notebooks that solve largely the same biomedical task with different interfaces:

*   [01 · Direct tools](notebooks/01_tools.ipynb): locally implemented Python functions exposed to an LLM via LangChain `@tool` (or called directly in the notebook).
*   [02 · MCP-hosted tools](notebooks/02_MCP.ipynb): the same functions exposed through a portable FastMCP server and loaded with langchain-mcp-adapters.
*   [03 · Agent skills](notebooks/03_agent_skills.ipynb): reusable behavioral instructions (a `SKILL.md`) that guide how the agent uses tools and reports evidence, run with the Pi harness.

See [info.ipynb](info.ipynb) for more information on the three tracks.