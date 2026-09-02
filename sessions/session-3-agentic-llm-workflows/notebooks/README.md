# Session 3: Agentic LLM workflow
Hands-on Practice Session

Session 3 main materials are split into three notebooks. Each one utilises the same queries suggested by Session 1 & 2, then changes the interface around the same data/tools/workflows.

For now we use the same mock task query for all three methods. We'll vary queries across three methods eventually to better demonstrate ideal use cases. 

## Example (mock) query

```
For profile P001, which genes connected to autism spectrum disorder should we prioritize, and why?
```

Expected toy-data answer:

```
SHANK3, SCN2A, CHD8
```

The ranking comes from profile scores among genes connected to phenotype node `HP:0000729`.

## Notebooks

`01_tools.ipynb`

*   Uses locally implemented Python functions directly.

`02_MCP.ipynb`

*   Shows how the same functions can be exposed through the MCP server.
*   We do need to figure out how to set up and port MCPs for participants.

`03_agent_skills.ipynb`

*   Shows how reusable skill instructions change the model's behaviour.
*   Local coding agent harness installation is required. (Pi for now)