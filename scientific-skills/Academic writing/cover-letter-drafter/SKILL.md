---
name: cover-letter-drafter
description: Generates professional cover letters for journal submissions and job
  applications in medical and academic contexts.
version: 1.0.0
category: Career
tags:
- cover-letter
- job-application
- journal-submission
- academic-writing
author: AIPOCH
license: MIT
status: Draft
risk_level: Medium
skill_type: Tool/Script
owner: AIPOCH
reviewer: ''
last_updated: '2026-02-06'
---

# Cover Letter Drafter

Creates tailored cover letters for academic and medical positions.

## Features

- Journal submission cover letters
- Job application cover letters
- Fellowship application letters
- Customizable templates

## Parameters

| Parameter | Type | Default | Required | Description |
|-----------|------|---------|----------|-------------|
| `--purpose` | string | job | No | Cover letter type (journal, job, fellowship) |
| `--recipient`, `-r` | string | - | Yes | Target journal or institution |
| `--key-points`, `-k` | string | - | Yes | Comma-separated key points to highlight |
| `--title` | string | - | No | Manuscript title (for journal submissions) |
| `--significance` | string | - | No | Significance statement (for journal submissions) |
| `--author`, `--applicant`, `-a` | string | Applicant | No | Author or applicant name |
| `--position` | string | - | No | Position title (for job applications) |
| `--fellowship` | string | - | No | Fellowship name (for fellowship applications) |
| `--output`, `-o` | string | - | No | Output JSON file path |

## Usage

```bash
# Journal submission cover letter
python scripts/main.py --purpose journal --recipient "Nature Medicine" \
  --key-points "Novel findings,Clinical relevance" \
  --title "Study X" --significance "major advance" --author "Dr. Smith"

# Job application cover letter
python scripts/main.py --purpose job --recipient "Harvard Medical School" \
  --key-points "10 years experience,Published 20 papers" \
  --position "Assistant Professor" --applicant "Dr. Jones"

# Fellowship application
python scripts/main.py --purpose fellowship --recipient "NIH" \
  --key-points "Research excellence,Leadership skills" \
  --fellowship "K99" --applicant "Dr. Lee"
```

## Output Format

```json
{
  "cover_letter": "string",
  "subject_line": "string",
  "word_count": "int"
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
