---
name: date-calculator
description: Calculates gestational age and follow-up date windows.
version: 1.0.0
category: Utility
tags:
- dates
- pregnancy
- follow-up
- calculator
author: AIPOCH
license: MIT
status: Draft
risk_level: Medium
skill_type: Tool/Script
owner: AIPOCH
reviewer: ''
last_updated: '2026-02-06'
---

# Date Calculator

Calculates medical date windows.

## Features

- Gestational age
- Follow-up windows
- Visit scheduling
- Date adjustments

## Parameters

| Parameter | Type | Default | Required | Description |
|-----------|------|---------|----------|-------------|
| `--type`, `-t` | string | - | Yes | Calculation type (gestational or followup) |
| `--date`, `-d` | string | - | Yes | Date in YYYY-MM-DD format |
| `--weeks` | int | 4 | No | Number of weeks for follow-up |
| `--window-days` | int | 7 | No | Follow-up window size in days |
| `--output`, `-o` | string | - | No | Output JSON file path |

## Usage

```bash
# Calculate gestational age
python scripts/main.py --type gestational --date 2024-01-15

# Calculate 4-week follow-up window
python scripts/main.py --type followup --date 2024-03-01

# Calculate custom follow-up (6 weeks)
python scripts/main.py --type followup --date 2024-03-01 --weeks 6
```

## Output Format

**Gestational calculation:**
```json
{
  "lmp_date": "2024-01-15",
  "gestational_age": "12 weeks 3 days",
  "gestational_age_days": 87,
  "estimated_delivery_date": "2024-10-21",
  "calculation_date": "2024-04-12"
}
```

**Follow-up calculation:**
```json
{
  "start_date": "2024-03-01",
  "followup_weeks": 4,
  "window_start": "2024-03-29",
  "window_end": "2024-04-05",
  "window_range": "2024-03-29 to 2024-04-05"
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
