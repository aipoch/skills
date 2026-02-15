---
name: conference-poster-pitch
description: Generate elevator pitch for conference posters
version: 1.0.0
category: Present
tags: []
author: AIPOCH
license: MIT
status: Draft
risk_level: Medium
skill_type: Tool/Script
owner: AIPOCH
reviewer: ''
last_updated: '2026-02-06'
---

# Conference Poster Pitch

Generate elevator pitches for academic poster sessions.

## Parameters

| Parameter | Type | Default | Required | Description |
|-----------|------|---------|----------|-------------|
| `--poster-title`, `-t` | string | - | Yes | Poster title |
| `--duration`, `-d` | int | 60 | No | Pitch duration in seconds (30, 60, or 180) |

## Usage

```bash
# Generate 60-second pitch
python scripts/main.py --poster-title "CRISPR Therapy for Sickle Cell Disease" --duration 60

# Generate quick 30-second pitch
python scripts/main.py --poster-title "Novel Biomarkers in Cancer" --duration 30

# Generate detailed 3-minute pitch
python scripts/main.py --poster-title "AI in Drug Discovery" --duration 180
```

## Output

- 30s, 60s, and 3-minute pitch scripts
- Structured elevator pitch format
- Ready-to-practice delivery text

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
