---
name: forest-plot-styler
description: Beautify meta-analysis forest plots with customizable odds ratio points
  and confidence intervals
version: 1.0.0
category: Visual
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

# Forest Plot Styler

ID: 157

Beautifies Meta-analysis or subgroup analysis forest plots, customizes Odds Ratio point sizes and confidence interval line styles.

---

## Features

- Reads Meta-analysis data (CSV/Excel format)
- Draws high-quality forest plots
- Customizes Odds Ratio point sizes, colors, and shapes
- Customizes confidence interval line styles (color, thickness, endpoint style)
- Supports subgroup analysis display
- Automatically calculates and displays pooled effect values
- Outputs to PNG, PDF, or SVG format

---

## Usage

```bash
python scripts/main.py --input <data.csv> [options]
```

### Parameters

| Parameter | Type | Default | Required | Description |
|-----------|------|---------|----------|-------------|
| `--input`, `-i` | string | - | Yes | Input data file (CSV or Excel) |
| `--output`, `-o` | string | forest_plot.png | No | Output file path |
| `--format`, `-f` | string | png | No | Output format (png/pdf/svg) |
| `--point-size` | int | 8 | No | OR point size |
| `--point-color` | string | #2E86AB | No | OR point color |
| `--ci-color` | string | #2E86AB | No | Confidence interval line color |
| `--ci-linewidth` | int | 2 | No | Confidence interval line thickness |
| `--ci-capwidth` | int | 5 | No | Confidence interval endpoint width |
| `--summary-color` | string | #A23B72 | No | Pooled effect point color |
| `--summary-shape` | string | diamond | No | Pooled effect point shape |
| `--subgroup` | string | - | No | Subgroup analysis column name |
| `--title`, `-t` | string | Forest Plot | No | Chart title |
| `--xlabel`, `-x` | string | Odds Ratio (95% CI) | No | X-axis label |
| `--reference-line` | float | 1.0 | No | Reference line position |
| `--width`, `-W` | int | 12 | No | Image width (inches) |
| `--height`, `-H` | int | auto | No | Image height (inches) |
| `--dpi` | int | 300 | No | Image resolution |
| `--font-size` | int | 10 | No | Font size |
| `--style`, `-s` | string | default | No | Preset style (default/minimal/dark) |

---

## Input Data Format

CSV/Excel files must contain the following columns:

| Column Name | Description | Type |
|------|------|------|
| `study` | Study name | Text |
| `or` | Odds Ratio value | Numeric |
| `ci_lower` | Confidence interval lower bound | Numeric |
| `ci_upper` | Confidence interval upper bound | Numeric |
| `weight` | Weight (optional, for point size) | Numeric |
| `subgroup` | Subgroup label (optional) | Text |

### Sample Data

```csv
study,or,ci_lower,ci_upper,weight,subgroup
Study A,0.85,0.65,1.12,15.2,Drug A
Study B,0.72,0.55,0.94,18.5,Drug A
Study C,1.15,0.88,1.50,12.3,Drug B
Study D,0.95,0.75,1.20,14.8,Drug B
```

---

## Examples

### Basic Usage

```bash
python scripts/main.py -i meta_data.csv
```

### Custom Style

```bash
python scripts/main.py -i meta_data.csv \
    --point-color="#E63946" \
    --ci-color="#457B9D" \
    --point-size=10 \
    --ci-linewidth=3 \
    -t "Meta-Analysis of Treatment Effects"
```

### Subgroup Analysis

```bash
python scripts/main.py -i meta_data.csv \
    --subgroup subgroup_column \
    --summary-color="#F4A261" \
    -o subgroup_forest.png
```

### Output PDF Vector Graphic

```bash
python scripts/main.py -i meta_data.csv \
    -f pdf \
    -o forest_plot.pdf
```

---

## Preset Styles

### default
- Blue color scheme
- Standard font size
- White background

### minimal
- Clean lines
- Grayscale color scheme
- No grid lines

### dark
- Dark background
- Bright data points
- Suitable for dark theme presentations

---

## Dependencies

- Python >= 3.8
- matplotlib >= 3.5.0
- pandas >= 1.3.0
- numpy >= 1.20.0
- openpyxl >= 3.0.0 (for reading Excel)

---

## Output Example

Generated forest plot contains:
- Left side: Study name list
- Middle: OR values and confidence intervals
- Right side: Weight percentage (if available)
- Bottom: Pooled effect value (diamond marker)
- Reference line (OR=1)

---

## Notes

1. Ensure input file encoding is UTF-8
2. OR values are automatically converted when log scale is suggested
3. Studies with confidence intervals crossing 1 are not statistically significant
4. Weight values are used to adjust point size, reflecting study contribution

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

```bash
# Python dependencies
pip install -r requirements.txt
```

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
