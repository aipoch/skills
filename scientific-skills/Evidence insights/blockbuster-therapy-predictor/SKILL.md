---
name: blockbuster-therapy-predictor
description: "Predict which early-stage biotechnology platforms (PROTAC, mRNA, gene 
  editing, etc.) have the highest potential to become blockbuster therapies.
  Analyzes clinical trial progression, patent landscape maturity, and venture 
  capital funding trends to generate investment and R&D prioritization scores.
  Trigger when: User asks about technology investment potential, platform 
  selection, or therapeutic modality comparison."
version: 1.0.0
category: Pharma
tags: ["investment", "prediction", "biotech", "clinical-trials", "patents"]
author: AIPOCH
license: MIT
status: Draft
risk_level: Medium
skill_type: Tool/Script
owner: AIPOCH
reviewer: ''
last_updated: '2026-02-15'
---

# Blockbuster Therapy Predictor

Comprehensive analytics tool for forecasting breakthrough therapeutic technologies by integrating multi-dimensional data sources including clinical development pipelines, intellectual property landscapes, and capital market indicators.

## Features

- **Multi-Source Data Integration**: Aggregates clinical trials, patents, and funding data
- **Predictive Scoring**: Calculates Blockbuster Index combining maturity, market potential, and momentum
- **Technology Landscape Mapping**: Tracks 10+ emerging therapeutic platforms
- **Investment Intelligence**: Provides data-driven R&D and investment recommendations
- **Trend Analysis**: Identifies acceleration patterns and inflection points

## Usage

### Basic Usage

```bash
# Run complete analysis with all technologies
python scripts/main.py

# Analyze specific technologies
python scripts/main.py --tech PROTAC,mRNA,CRISPR

# Output in JSON format
python scripts/main.py --output json
```

### Parameters

| Parameter | Type | Default | Required | Description |
|-----------|------|---------|----------|-------------|
| `--mode` | str | full | No | Analysis mode: full or quick |
| `--tech` | str | None | No | Comma-separated list of technologies to analyze |
| `--output` | str | console | No | Output format: console or json |
| `--threshold` | float | 0 | No | Minimum blockbuster index threshold (0-100) |
| `--save` | str | None | No | Save report to file path |

### Advanced Usage

```bash
# Analyze high-potential technologies only (index ≥70)
python scripts/main.py \
  --threshold 70 \
  --output json \
  --save high_potential_report.json

# Quick analysis of specific platforms
python scripts/main.py \
  --mode quick \
  --tech CAR-T,ADC,Bispecific \
  --output console
```

## Output

### Console Output

```
🏆 BLOCKBUSTER THERAPY PREDICTOR Report
Generated: 2026-02-15 10:30:00
Technologies analyzed: 10

📊 Technology Rankings
Rank  Technology       Blockbuster Index    Maturity    Market Potential    Momentum    Recommendation
🥇 1   mRNA             85.2                 78.5        92.1                88.0        Strongly Recommended
🥈 2   CAR-T            82.3                 85.2        78.5                75.0        Strongly Recommended
🥉 3   CRISPR           79.8                 72.3        88.2                68.0        Recommended
```

### JSON Output Structure

```json
{
  "generated_at": "2026-02-15T10:30:00",
  "total_routes": 10,
  "rankings": [
    {
      "rank": 1,
      "tech_name": "mRNA",
      "blockbuster_index": 85.2,
      "maturity_score": 78.5,
      "market_potential_score": 92.1,
      "momentum_score": 88.0,
      "recommendation": "Strongly Recommended",
      "key_drivers": ["Multiple Phase III trials", "Rapid patent growth"],
      "risk_factors": ["Regulatory uncertainties"],
      "timeline_prediction": "First product expected in 2-4 years"
    }
  ]
}
```

## Scoring Methodology

### Blockbuster Index Formula

```
Blockbuster Index = (Market Potential × 0.5) + (Maturity × 0.3) + (Momentum × 0.2)
```

### Component Scores

| Component | Weight | Factors |
|-----------|--------|---------|
| **Market Potential** | 50% | Market size, unmet need, competition |
| **Maturity** | 30% | Clinical stage, patent depth, funding stage |
| **Momentum** | 20% | Patent growth, funding activity, clinical progress |

### Investment Recommendation Thresholds

| Blockbuster Index | Recommendation | Action |
|-------------------|----------------|--------|
| ≥ 80 | **Strongly Recommended** | Prioritize R&D investment |
| 60-79 | **Recommended** | Active monitoring and early partnerships |
| 40-59 | **Watch** | Monitor milestones; reassess in 6-12 months |
| < 40 | **Cautious** | Minimal investment; consider divestment |

## Supported Technologies

| Technology | Category | Description |
|------------|----------|-------------|
| PROTAC | Protein Degradation | Proteolysis Targeting Chimera |
| mRNA | Nucleic Acid Drugs | Messenger RNA therapy platform |
| CRISPR | Gene Editing | CRISPR-Cas gene editing technology |
| CAR-T | Cell Therapy | Chimeric Antigen Receptor T-cell therapy |
| Bispecific | Antibody Drugs | Bispecific antibody technology |
| ADC | Antibody Drugs | Antibody-Drug Conjugate |
| RNAi | Nucleic Acid Drugs | RNA interference therapy |
| Gene Therapy | Gene Therapy | AAV vector gene therapy |
| Allogeneic | Cell Therapy | Universal/Allogeneic cell therapy |
| Cell Therapy | Cell Therapy | General cell therapy platform |

## Technical Difficulty: **MEDIUM**

⚠️ **AI自主验收状态**: 需人工检查

This skill requires:
- Python 3.8+ environment
- Basic understanding of biotech investment analysis
- Access to clinical trial, patent, and funding databases (optional)

## Dependencies

### Required Python Packages

```bash
pip install -r requirements.txt
```

### Requirements File

```
dataclasses
enum
```

## Risk Assessment

| Risk Indicator | Assessment | Level |
|----------------|------------|-------|
| Code Execution | Python scripts executed locally | Medium |
| Network Access | No external API calls in mock mode | Low |
| File System Access | Read/write report files only | Low |
| Instruction Tampering | Standard prompt guidelines | Low |
| Data Exposure | Output files saved to workspace | Low |

## Security Checklist

- [x] No hardcoded credentials or API keys
- [x] No unauthorized file system access (../)
- [x] Output does not expose sensitive information
- [x] Prompt injection protections in place
- [x] Input file paths validated (no ../ traversal)
- [x] Output directory restricted to workspace
- [x] Script execution in sandboxed environment
- [x] Error messages sanitized (no stack traces exposed)
- [x] Dependencies audited

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
1. **Basic Functionality**: Run without arguments → Expected output with all technologies
2. **Technology Filter**: Use --tech flag → Only specified technologies analyzed
3. **JSON Output**: Use --output json → Valid JSON format output
4. **Threshold Filter**: Use --threshold 70 → Only technologies with index ≥70 shown

## Lifecycle Status

- **Current Stage**: Draft
- **Next Review Date**: 2026-03-15
- **Known Issues**: None
- **Planned Improvements**: 
  - Integration with real-time data APIs
  - Additional technology platforms
  - Enhanced visualization capabilities

## References

See `references/` for:
- Historical blockbuster case studies
- Clinical trial data sources
- Patent analysis methodologies
- Investment scoring frameworks

## Limitations

- **Data Source**: Uses mock data for demonstration; real-time data integration required for production use
- **Prediction Accuracy**: Model provides indicative scores; not investment advice
- **Technology Coverage**: Limited to pre-configured technology platforms
- **Market Dynamics**: Cannot predict black swan events or regulatory changes
- **Regional Bias**: Data primarily focused on US/EU markets

---

**⚠️ DISCLAIMER: This tool provides quantitative analysis for decision support only. All investment and R&D decisions should incorporate qualitative domain expertise, regulatory consultation, and comprehensive due diligence. Past performance of historical blockbusters does not guarantee future success of emerging technologies.**
