---
name: adaptive-trial-simulator
description: "Design and simulate adaptive clinical trials with interim analyses, 
  sample size re-estimation, and early stopping rules. Evaluate Type I error control, 
  power, and expected sample size via Monte Carlo simulation before trial initiation."
version: 1.0.0
category: Clinical
tags: ["clinical-trials", "adaptive-design", "statistics", "simulation", "biostatistics"]
author: AIPOCH
license: MIT
status: Draft
risk_level: Medium
skill_type: Tool/Script
owner: AIPOCH
reviewer: ''
last_updated: '2026-02-15'
---

# Adaptive Trial Simulator

Statistical simulation platform for designing and validating adaptive clinical trial designs in silico. Enables optimization of interim analysis strategies, sample size adaptation, and early stopping rules while maintaining Type I error control.

## Features

- **Design Simulation**: Monte Carlo validation of adaptive designs
- **Sample Size Re-estimation**: Adapt sample size based on interim data
- **Early Stopping Rules**: Futility and efficacy boundary optimization
- **Type I Error Control**: Validate alpha spending strategies
- **Multi-Arm Designs**: Drop-the-loser and seamless Phase II/III
- **Power Optimization**: Identify designs with maximum power efficiency

## Usage

### Basic Usage

```bash
# Run standard group sequential design
python scripts/main.py

# Adaptive design with sample size re-estimation
python scripts/main.py --design adaptive_reestimate

# Optimize design parameters
python scripts/main.py --optimize
```

### Parameters

| Parameter | Type | Default | Required | Description |
|-----------|------|---------|----------|-------------|
| `--design` | str | group_sequential | No | Trial design type |
| `--n-simulations` | int | 10000 | No | Number of Monte Carlo simulations |
| `--sample-size` | int | 200 | No | Initial sample size per arm |
| `--effect-size` | float | 0.3 | No | Effect size (Cohen's d) |
| `--alpha` | float | 0.05 | No | Type I error rate |
| `--power` | float | 0.80 | No | Target statistical power |
| `--interim-looks` | int | 1 | No | Number of interim analyses |
| `--spending-function` | str | obrien_fleming | No | Alpha spending function |
| `--reestimate-method` | str | promising_zone | No | Sample size re-estimation method |
| `--output` | str | results.json | No | Output file path |
| `--visualize` | flag | False | No | Generate visualization charts |
| `--optimize` | flag | False | No | Search for optimal design parameters |

### Advanced Usage

```bash
# Full adaptive design with visualization
python scripts/main.py \
  --design adaptive_reestimate \
  --n-simulations 50000 \
  --sample-size 250 \
  --effect-size 0.35 \
  --interim-looks 2 \
  --spending-function obrien_fleming \
  --visualize \
  --output adaptive_results.json
```

## Design Types

| Design Type | Description | Use Case |
|-------------|-------------|----------|
| **Group Sequential** | Fixed interim looks with stopping boundaries | Standard adaptive trials |
| **Adaptive Re-estimate** | Sample size adjustment based on interim data | Uncertain effect size |
| **Drop the Loser** | Multi-arm trials dropping inferior arms | Phase II dose selection |

## Spending Functions

| Function | Characteristics | Early Boundary |
|----------|----------------|----------------|
| **O'Brien-Fleming** | Conservative early | High Z-scores early |
| **Pocock** | Aggressive early | Lower Z-scores throughout |
| **Power Family** | Moderate (ρ=3) | Balanced approach |

## Output Example

```json
{
  "design_config": {
    "design_type": "adaptive_reestimate",
    "sample_size_per_arm": 200,
    "effect_size": 0.3,
    "alpha": 0.05,
    "target_power": 0.8
  },
  "simulation_results": {
    "power": 0.8234,
    "type_i_error": 0.0481,
    "expected_sample_size": 385.2,
    "early_stop_rate": {
      "efficacy": 0.1523,
      "futility": 0.0841
    }
  }
}
```

## Technical Difficulty: **HIGH**

⚠️ **AI自主验收状态**: 需人工检查

This skill requires:
- Python 3.8+ environment
- NumPy, SciPy, and Matplotlib packages
- Understanding of clinical trial statistics

## Dependencies

```bash
pip install -r requirements.txt
```

### Requirements

```
numpy>=1.20.0
scipy>=1.7.0
matplotlib>=3.4.0
```

## Risk Assessment

| Risk Indicator | Assessment | Level |
|----------------|------------|-------|
| Code Execution | Python scripts with mathematical calculations | Medium |
| Network Access | No network access | Low |
| File System Access | Writes simulation results | Low |
| Instruction Tampering | Statistical parameters could affect results | Medium |
| Data Exposure | No sensitive data exposure | Low |

## Security Checklist

- [x] No hardcoded credentials or API keys
- [x] No unauthorized file system access
- [x] Output does not expose sensitive information
- [x] Input parameters validated
- [x] Error messages sanitized
- [x] Dependencies audited

## Prerequisites

```bash
pip install -r requirements.txt
python scripts/main.py --help
```

## Evaluation Criteria

### Success Metrics
- [ ] Simulations run without errors
- [ ] Type I error controlled at nominal level
- [ ] Power estimates are accurate
- [ ] Visualizations generated correctly

### Test Cases
1. **Basic Simulation**: Default parameters → Valid results
2. **Different Designs**: All design types → Appropriate behavior
3. **Optimization Mode**: --optimize flag → Finds optimal parameters
4. **Visualization**: --visualize flag → Charts generated

## Lifecycle Status

- **Current Stage**: Draft
- **Next Review Date**: 2026-03-15
- **Known Issues**: Type checking warnings with numpy arrays
- **Planned Improvements**: 
  - Bayesian adaptive designs
  - Multi-arm multi-stage (MAMS) support
  - Enhanced visualization options

## References

Available in `references/`:
- Adaptive design statistical theory
- Regulatory guidance documents
- Alpha spending function literature
- Sample size re-estimation methods

## Limitations

- **Statistical Complexity**: Requires biostatistics expertise
- **Simulation Time**: Large simulations may take hours
- **Simplified Models**: Does not capture all real-world complexities
- **Regulatory Consultation**: Results should be validated with regulators

---

**⚠️ DISCLAIMER: This tool provides simulation results for research and planning purposes only. All clinical trial designs should be reviewed by qualified biostatisticians and regulatory experts before implementation.**
