---
name: abstract-trimmer
description: "Compress academic abstracts to meet strict word limits while preserving 
  key information, scientific accuracy, and readability. Supports multiple compression 
  strategies for journal submissions, conference applications, and grant proposals."
version: 1.0.0
category: General
tags: ["abstract", "editing", "compression", "academic-writing", "word-count"]
author: AIPOCH
license: MIT
status: Draft
risk_level: Medium
skill_type: Tool/Script
owner: AIPOCH
reviewer: ''
last_updated: '2026-02-15'
---

# Abstract Trimmer

Precision editing tool that reduces abstract word count through intelligent compression techniques, maintaining scientific rigor while meeting strict journal and conference requirements.

## Features

- **Smart Compression**: Multiple strategies (aggressive, conservative, balanced)
- **Key Information Preservation**: Retains critical findings and statistics
- **Structural Integrity**: Maintains Background-Methods-Results-Conclusion flow
- **Quantitative Safety**: Protects numbers, P-values, and confidence intervals
- **Batch Processing**: Trim multiple abstracts efficiently
- **Quality Validation**: Post-trim readability and accuracy checks

## Usage

### Basic Usage

```bash
# Trim abstract from file
python scripts/main.py --input abstract.txt --target 250

# Trim abstract from command line
python scripts/main.py --text "Your abstract here..." --target 200

# Check word count only
python scripts/main.py --input abstract.txt --target 250 --check-only
```

### Parameters

| Parameter | Type | Default | Required | Description |
|-----------|------|---------|----------|-------------|
| `--input`, `-i` | str | None | No | Input file containing abstract |
| `--text`, `-t` | str | None | No | Abstract text (alternative to --input) |
| `--target`, `-T` | int | 250 | No | Target word count |
| `--strategy`, `-s` | str | balanced | No | Trimming strategy (conservative/balanced/aggressive) |
| `--output`, `-o` | str | None | No | Output file path |
| `--check-only`, `-c` | flag | False | No | Only check word count without trimming |
| `--format` | str | json | No | Output format (json/text) |

### Advanced Usage

```bash
# Aggressive trimming with text output
python scripts/main.py \
  --input abstract.txt \
  --target 200 \
  --strategy aggressive \
  --format text \
  --output trimmed.txt

# Batch check multiple abstracts
for file in *.txt; do
  python scripts/main.py --input "$file" --target 250 --check-only
done
```

## Trimming Strategies

| Strategy | Approach | Best For |
|----------|----------|----------|
| **Conservative** | Remove filler words, simplify sentences | Minor trims (10-20 words) |
| **Balanced** | Condense phrases, merge sentences | Moderate trims (20-50 words) |
| **Aggressive** | Remove secondary details, abbreviate | Major trims (50+ words) |

## Output Format

### JSON Output

```json
{
  "trimmed_abstract": "Compressed abstract text...",
  "original_words": 320,
  "final_words": 248,
  "reduction_percent": 22.5
}
```

### Text Output

```
Compressed abstract text...
```

## Technical Difficulty: **LOW**

⚠️ **AI自主验收状态**: 需人工检查

This skill requires:
- Python 3.7+ environment
- No external dependencies

## Dependencies

### Required Python Packages

```bash
pip install -r requirements.txt
```

### Requirements File

No external dependencies required (uses only Python standard library).

## Risk Assessment

| Risk Indicator | Assessment | Level |
|----------------|------------|-------|
| Code Execution | Python scripts executed locally | Low |
| Network Access | No network access | Low |
| File System Access | Read/write text files only | Low |
| Instruction Tampering | Standard prompt guidelines | Low |
| Data Exposure | No sensitive data exposure | Low |

## Security Checklist

- [x] No hardcoded credentials or API keys
- [x] No unauthorized file system access (../)
- [x] Output does not expose sensitive information
- [x] Prompt injection protections in place
- [x] Input file paths validated
- [x] Output directory restricted to workspace
- [x] Script execution in sandboxed environment
- [x] Error messages sanitized
- [x] Dependencies audited

## Prerequisites

```bash
# No dependencies required
python scripts/main.py --help
```

## Evaluation Criteria

### Success Metrics
- [ ] Successfully trims abstracts to target word count
- [ ] Preserves key scientific information
- [ ] Maintains grammatical correctness
- [ ] Handles edge cases gracefully

### Test Cases
1. **Basic Trimming**: Input abstract → Trimed to target word count
2. **Check Mode**: --check-only flag → Reports word count statistics
3. **File I/O**: Read from file, write to file → Correct file handling
4. **Different Strategies**: All three strategies work → Different compression levels

## Lifecycle Status

- **Current Stage**: Draft
- **Next Review Date**: 2026-03-15
- **Known Issues**: None
- **Planned Improvements**: 
  - Enhanced protection for quantitative data
  - Support for structured abstracts
  - Batch processing mode

## References

See `references/` for:
- Compression strategies documentation
- Protected elements guidelines
- Journal word limits by publisher

## Limitations

- **Language**: Optimized for English academic abstracts
- **Content Type**: Designed for structured abstracts (BMRC format)
- **No Rewriting**: Only removes/compresses; doesn't rephrase
- **Final Review**: Automated trimming requires human validation

---

**✂️ Remember: This tool helps meet word limits, but never sacrifice scientific accuracy. Always validate that trimmed abstracts maintain the integrity of your findings.**
