# MCP Server

This folder contains the optional MCP server material for the Session 3 practical.

## Files

*   `toy_kg_mcp_server.py`: minimal FastMCP server exposing the toy graph/profile  
    tools (plus a `kg://summary` resource and a `prioritize_genes` prompt).
*   `../notebooks/02_MCP.ipynb`: walkthrough of connecting to this  
    server, both raw and via langchain-mcp-adapters.

## Intended Role

The server exposes the same graph/profile operations used in the notebook:

*   `search_nodes`
*   `get_neighbors`
*   `rank_profile_features`
*   `rank_candidate_genes`

The teaching point is that the biomedical tools can live outside the notebook and be reused by different LLM clients.

## Quick Start

From the project root (using the `eccb` conda env the notebooks use, or any  
env with `requirements.txt` installed):

```
python -m pip install -r requirements.txt
python server/toy_kg_mcp_server.py
```

In practice you rarely launch it by hand — `notebooks/02_MCP.ipynb`  
spawns it as a stdio subprocess.