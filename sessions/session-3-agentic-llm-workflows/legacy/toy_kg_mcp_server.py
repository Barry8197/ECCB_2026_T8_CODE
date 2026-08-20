"""Minimal MCP server for the ECCB 2026 toy knowledge graph practical.

Run with:

    python server/toy_kg_mcp_server.py

The server uses stdio transport, which is the usual local MCP pattern for
desktop clients and tutorial demos.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

from mcp.server.fastmcp import FastMCP
from mcp.types import ToolAnnotations

# Make `src` importable no matter what working directory the MCP client
# launches this script from (clients usually spawn it as a subprocess).
PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src import graph_tools  # noqa: E402  (import after sys.path setup)


mcp = FastMCP("eccb2026-toy-kg")
GRAPH = graph_tools.build_toy_graph()

# All of our tools only *read* the graph — they never mutate state. MCP lets us
# advertise that to clients via tool annotations (read-only / non-destructive /
# idempotent). A client can use these LLM-facing hints to decide, e.g., whether
# to auto-run a tool or ask the human first. These hints are exactly the kind of
# metadata a generic RPC framework would have no concept of.
READ_ONLY = ToolAnnotations(readOnlyHint=True, destructiveHint=False, idempotentHint=True)


# ---------------------------------------------------------------------------
# RESOURCE primitive: data a client loads *into context* (read, not "called").
# ---------------------------------------------------------------------------
@mcp.resource("kg://summary")
def graph_summary() -> dict[str, Any]:
    """Return a compact summary of the toy graph."""
    node_types: dict[str, int] = {}
    for _, attrs in GRAPH.nodes(data=True):
        node_type = str(attrs.get("type", "unknown"))
        node_types[node_type] = node_types.get(node_type, 0) + 1

    return {
        "name": "ECCB 2026 toy biomedical knowledge graph",
        "nodes": GRAPH.number_of_nodes(),
        "edges": GRAPH.number_of_edges(),
        "node_types": node_types,
        "available_profiles": sorted(graph_tools.TOY_PROFILES),
    }


# ---------------------------------------------------------------------------
# PROMPT primitive: a reusable, parameterized prompt template the server offers
# to clients (think slash-command). The server owns the wording; the client
# fills in the arguments. This is *not* a tool — the model does not call it; a
# user/client selects it.
# ---------------------------------------------------------------------------
@mcp.prompt(title="Prioritize candidate genes")
def prioritize_genes(profile_id: str, phenotype: str) -> str:
    """Reusable prompt template for the gene-prioritization task."""
    return (
        f"For profile {profile_id}, which genes connected to {phenotype} should "
        f"we prioritize, and why? Use the available tools to resolve the "
        f"phenotype id, rank candidate genes by profile score, and ground every "
        f"claim in tool evidence."
    )


# ---------------------------------------------------------------------------
# TOOL primitive: functions the model may choose to call.
# ---------------------------------------------------------------------------
@mcp.tool(annotations=READ_ONLY)
def search_nodes(query: str) -> list[dict[str, Any]]:
    """Search graph nodes by identifier or display name."""
    return graph_tools.search_nodes(GRAPH, query)


@mcp.tool(annotations=READ_ONLY)
def get_neighbors(
    node_id: str,
    node_type: str | None = None,
) -> list[dict[str, Any]]:
    """Return incoming and outgoing graph neighbors for a node."""
    return graph_tools.get_neighbors(GRAPH, node_id, node_type=node_type)


@mcp.tool(annotations=READ_ONLY)
def rank_profile_features(profile_id: str, top_n: int = 20) -> list[dict[str, Any]]:
    """Return the top molecular features for a toy multi-omic profile."""
    return graph_tools.rank_profile_features(profile_id, top_n=top_n)


@mcp.tool(annotations=READ_ONLY)
def rank_candidate_genes(
    profile_id: str,
    phenotype_id: str,
) -> list[dict[str, Any]]:
    """Rank phenotype-linked genes by their toy profile score."""
    return graph_tools.rank_candidate_genes(GRAPH, profile_id, phenotype_id)


if __name__ == "__main__":
    mcp.run()
