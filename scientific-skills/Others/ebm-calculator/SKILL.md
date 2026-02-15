---
name: ebm-calculator
description: Evidence-Based Medicine calculator for sensitivity, specificity, PPV,
  NPV, NNT, and likelihood ratios. Essential for clinical decision making and biostatistics
  education.
version: 1.0.0
category: Education
tags:
- ebm
- biostatistics
- calculator
- sensitivity
- specificity
- clinical-reasoning
author: AIPOCH
license: MIT
status: Draft
risk_level: Medium
skill_type: Tool/Script
owner: AIPOCH
reviewer: ''
last_updated: '2026-02-06'
---

# EBM Calculator

Evidence-Based Medicine diagnostic test calculator.

## Features

- Sensitivity / Specificity calculation
- PPV / NPV with prevalence adjustment
- Likelihood ratios (LR+ / LR-)
- Number Needed to Treat (NNT)
- Pre/post-test probability conversion

## Parameters

| Parameter | Type | Default | Required | Description |
|-----------|------|---------|----------|-------------|
| `--mode`, `-m` | string | diagnostic | No | Calculation mode (diagnostic, nnt, probability) |
| `--tp`, `--true-pos` | int | - | * | True positives (diagnostic mode) |
| `--fn`, `--false-neg` | int | - | * | False negatives (diagnostic mode) |
| `--tn`, `--true-neg` | int | - | * | True negatives (diagnostic mode) |
| `--fp`, `--false-pos` | int | - | * | False positives (diagnostic mode) |
| `--prevalence`, `-p` | float | - | No | Disease prevalence 0-1 (diagnostic mode) |
| `--control-rate` | float | - | ** | Control event rate 0-1 (nnt mode) |
| `--experimental-rate` | float | - | ** | Experimental event rate 0-1 (nnt mode) |
| `--pretest` | float | - | *** | Pre-test probability 0-1 (probability mode) |
| `--lr` | float | - | *** | Likelihood ratio (probability mode) |
| `--output`, `-o` | string | stdout | No | Output file path |

\* Required for diagnostic mode  
\** Required for nnt mode  
\*** Required for probability mode

## Output Format

```json
{
  "sensitivity": "float",
  "specificity": "float",
  "ppv": "float",
  "npv": "float",
  "lr_positive": "float",
  "lr_negative": "float",
  "interpretation": "string"
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
