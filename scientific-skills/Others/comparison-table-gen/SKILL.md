---
name: comparison-table-gen
description: Auto-generates comparison tables for concepts, drugs, or study results
  in Markdown format.
version: 1.0.0
category: Info
tags:
- comparison
- table
- markdown
- research
author: AIPOCH
license: MIT
status: Draft
risk_level: Medium
skill_type: Tool/Script
owner: AIPOCH
reviewer: ''
last_updated: '2026-02-06'
---

# Comparison Table Gen

Generates comparison tables for medical content.

## Features

- Side-by-side comparisons
- Markdown table output
- Drug comparison templates
- Study result comparisons

## Parameters

| Parameter | Type | Default | Required | Description |
|-----------|------|---------|----------|-------------|
| `--items`, `-i` | string | - | Yes | Items to compare (comma-separated) |
| `--attributes`, `-a` | string | - | Yes | Comparison attributes (comma-separated) |
| `--output`, `-o` | string | - | No | Output JSON file path |

## Usage

```bash
# Compare two drugs
python scripts/main.py --items "Drug A,Drug B" --attributes "Mechanism,Dose,Side Effects"

# Save to file
python scripts/main.py --items "Surgery,Chemo,Radiation" --attributes "Cost,Efficacy" --output comparison.json
```

## Input Format

- **items**: Comma-separated list of items to compare (e.g., "Drug A,Drug B")
- **attributes**: Comma-separated list of comparison attributes (e.g., "Mechanism,Dose")

## Output Format

```json
{
  "markdown_table": "string",
  "html_table": "string"
}
```

## Risk Assessment

| Risk Indicator | Assessment | Level |
|----------------|------------|-------|
| Code Execution | Python/R scripts executed locally | Medium |
| Network Access | No external API calls | Low |
| File System Access | Read input files, write output files | Medium |
| Instruction Tampering | Standard prompt guidelines | Low |
| Data Exposure | Output files saved to workspace | Low |

## Security Checklist

- [ ] No hardcoded credentials or API keys
- [ ] No unauthorized file system access (../)
- [ ] Output does not expose sensitive information
- [ ] Prompt injection protections in place
- [ ] Input file paths validated (no ../ traversal)
- [ ] Output directory restricted to workspace
- [ ] Script execution in sandboxed environment
- [ ] Error messages sanitized (no stack traces exposed)
- [ ] Dependencies audited
## Prerequisites

No additional Python packages required.

## Evaluation Criteria

### Success Metrics
- [ ] Successfully executes main functionality
- [ ] Output meets quality standards
- [ ] Handles edge cases gracefully
- [ ] Performance is acceptable

### Test Cases
1. **Basic Functionality**: Standard input → Expected output
2. **Edge Case**: Invalid input → Graceful error handling
3. **Performance**: Large dataset → Acceptable processing time

## Lifecycle Status

- **Current Stage**: Draft
- **Next Review Date**: 2026-03-06
- **Known Issues**: None
- **Planned Improvements**: 
  - Performance optimization
  - Additional feature support
