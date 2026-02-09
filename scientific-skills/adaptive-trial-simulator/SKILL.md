---
name: adaptive-trial-simulator
description: Design and simulate adaptive clinical trials with interim analyses, 
  sample size re-estimation, and early stopping rules. Evaluate Type I error control, 
  power, and expected sample size via Monte Carlo simulation before trial initiation.
allowed-tools: [Read, Write, Bash, Edit]
license: MIT
metadata:
    skill-author: AIPOCH
---

# Adaptive Trial Simulator

## Overview

Statistical simulation platform for designing and validating adaptive clinical trial designs in silico. Enables optimization of interim analysis strategies, sample size adaptation, and early stopping rules while maintaining Type I error control.

**Key Capabilities:**
- **Design Simulation**: Monte Carlo validation of adaptive designs
- **Sample Size Re-estimation**: Adapt sample size based on interim data
- **Early Stopping Rules**: Futility and efficacy boundary optimization
- **Type I Error Control**: Validate alpha spending strategies
- **Multi-Arm Designs**: Drop-the-loser and seamless Phase II/III
- **Power Optimization**: Identify designs with maximum power efficiency

## When to Use

**✅ Use this skill when:**
- Designing adaptive clinical trials with interim analyses
- Evaluating sample size re-estimation strategies
- Optimizing early stopping boundaries for efficacy/futility
- Comparing group sequential vs. adaptive designs
- Validating Type I error control via simulation
- Planning multi-arm multi-stage (MAMS) trials
- Training biostatisticians on adaptive methods
- Responding to regulatory questions about design validity

**❌ Do NOT use when:**
- Fixed-design trial is sufficient → Use `sample-size-power-calculator`
- Simple A/B testing → Use standard power analysis
- Trial already running (post-hoc adaptations) → Consult statistician
- Basket or umbrella trial design → Use `master-protocol-designer`
- Complex Bayesian adaptive designs → Use specialized Bayesian tools

**Integration:**
- **Upstream**: `sample-size-power-calculator` (initial sample size), `clinical-data-cleaner` (data prep)
- **Downstream**: `regulatory-consultation-prep` (FDA briefing docs), `grant-proposal-assistant` (methods section)

## Core Capabilities

### 1. Group Sequential Design Simulation

Simulate trials with planned interim analyses:

```python
from scripts.simulator import AdaptiveSimulator

sim = AdaptiveSimulator()

# Configure group sequential design
config = {
    "design_type": "group_sequential",
    "sample_size_per_arm": 300,
    "interim_looks": [0.33, 0.67],  # Information fractions
    "alpha_spending": "obrien_fleming",
    "effect_size": 0.3,
    "alpha": 0.05,
    "power": 0.80
}

# Run simulations
results = sim.simulate(
    config=config,
    n_simulations=10000,
    seed=42
)

print(f"Type I Error: {results.type_i_error:.4f}")
print(f"Power: {results.power:.4f}")
print(f"Expected Sample Size: {results.expected_n:.1f}")
```

**Boundary Types:**
| Spending Function | Characteristics | Best For |
|-------------------|----------------|----------|
| **O'Brien-Fleming** | Conservative early, lenient late | Confirmatory trials |
| **Pocock** | Aggressive early stopping | Early efficacy trials |
| **Power Family** | Tunable flexibility | Custom designs |
| **Lan-DeMets** | Flexible information rates | Variable enrollment |

### 2. Sample Size Re-estimation

Adapt sample size based on promising interim results:

```python
# Promising zone design
ssr_config = {
    "design_type": "adaptive_reestimate",
    "initial_n": 200,
    "max_n": 400,
    "promising_zone": [0.1, 0.5],  # Conditional power range
    "reestimate_method": "promising_zone"
}

results = sim.simulate(ssr_config, n_simulations=10000)

# Output shows:
# - Probability of SSR triggered
# - Expected sample size under various true effects
# - Type I error inflation (should be controlled)
```

**SSR Methods:**
- **Promising Zone**: Increase N only if conditional power in 10-50% range
- **Conditional Power**: Re-estimate to achieve target CP
- **Inverse Normal**: Combination test for weighted evidence
- **Weighted Statistics**: Adjust test statistic for adaptation

### 3. Drop-the-Loser Design

Simulate multi-arm adaptive selection:

```python
# Multi-arm multi-stage (MAMS)
mams_config = {
    "design_type": "drop_the_loser",
    "n_arms": 3,
    "selection_criteria": "best_response",
    "seamless_phase23": True,
    "control_arm": "shared"
}

results = sim.simulate_mams(mams_config)

# Evaluates:
# - Probability of selecting best arm
# - Family-wise error rate control
# - Expected sample size savings
```

**Applications:**
- Phase II dose selection (seamless II/III)
- Multiple treatment comparisons
- Adaptive enrichment designs
- Platform trial simulations

### 4. Operating Characteristics Analysis

Comprehensive evaluation of design performance:

```python
# Compare multiple designs under various scenarios
scenarios = [
    {"effect": 0.0, "label": "Null"},
    {"effect": 0.2, "label": "Small"},
    {"effect": 0.3, "label": "Expected"},
    {"effect": 0.5, "label": "Large"}
]

comparison = sim.compare_designs(
    designs=["fixed", "gs_obrien", "gs_pocock", "adaptive_ssr"],
    scenarios=scenarios,
    metrics=["power", "ess", "type_i_error", "early_stop_rate"]
)

comparison.generate_report("design_comparison.pdf")
```

## Common Patterns

### Pattern 1: Confirmatory Trial with Early Stopping

**Scenario**: Phase III trial where early efficacy claim is valuable.

```bash
# O'Brien-Fleming design for regulatory acceptance
python scripts/main.py \
  --design group_sequential \
  --sample-size 400 \
  --interim-looks 2 \
  --spending-function obrien_fleming \
  --alpha 0.025 \
  --one-sided \
  --n-simulations 50000 \
  --output phase3_design.json
```

**Design Features:**
- Two interim looks (33% and 67% information)
- Conservative O'Brien-Fleming boundaries
- 2.5% one-sided alpha (regulatory standard)
- High probability of early stop if large effect

**Expected Performance:**
- Type I Error: ~2.5% (controlled)
- Power: ~85% (slight loss vs. fixed design)
- Expected N: ~320 (savings if early stop)
- Early stop rate: ~40% under alternative

### Pattern 2: Promising Zone SSR

**Scenario**: Uncertain effect size; want to increase power if interim results promising.

```python
config = {
    "design_type": "adaptive_reestimate",
    "initial_n": 150,
    "max_n": 300,
    "promising_zone": [0.15, 0.85],  # Conditional power
    "min_effect_for_ssr": 0.15,  # Don't increase if CP < 15%
    "target_cp": 0.90  # Re-estimate to achieve 90% CP
}

results = sim.simulate(config, n_simulations=100000)
```

**Key Metrics:**
- SSR triggered: ~30% of trials
- Type I error inflation: <0.5% (acceptable)
- Expected N under H0: 152 (minimal inflation)
- Expected N under H1: 215 (power boost)

### Pattern 3: Seamless Phase II/III

**Scenario**: Dose selection followed by confirmatory testing without stopping enrollment.

```bash
# 3-arm dose selection design
python scripts/main.py \
  --design drop_the_loser \
  --n-arms 3 \
  --doses ["low", "medium", "high"] \
  --selection-stage-n 100 \
  --confirmatory-stage-n 200 \
  --control-arm shared \
  --selection-criteria best_response \
  --mam-output seamless_phase23_report.pdf
```

**Advantages:**
- No gap between Phase II and III
- Shared control arm reduces sample size
- Seamless transition if no safety concerns

**Considerations:**
- Family-wise error rate must be controlled
- Selection bias if response endpoint ≠ confirmatory endpoint
- Regulatory acceptance requires careful planning

### Pattern 4: Futility Analysis Design

**Scenario**: Want to stop early for lack of efficacy to save resources.

```python
# Conservative futility only
config = {
    "design_type": "group_sequential",
    "boundaries": {
        "efficacy": "none",  # No early efficacy claim
        "futility": "binding"  # Mandatory stop if crossed
    },
    "futility_boundary": "hsd",  # Hwang-Shih-DeCani spending
    "gamma": -4  # Conservative parameter
}

results = sim.simulate(config)
```

**Benefits:**
- Stop unpromising trials early
- Reallocate resources to other projects
- Ethical: avoid exposing patients to ineffective treatment

**Trade-offs:**
- Power loss if futility boundary too aggressive
- Must follow stopping rules (binding)
- Requires strong evidence of futility

## Complete Workflow Example

**From design concept to regulatory submission:**

```bash
# Step 1: Simulate fixed design baseline
python scripts/main.py \
  --design fixed \
  --sample-size 400 \
  --effect-size 0.3 \
  --n-simulations 10000 \
  --output fixed_baseline.json

# Step 2: Simulate adaptive alternatives
python scripts/main.py \
  --design adaptive_reestimate \
  --initial-n 300 \
  --max-n 500 \
  --promising-zone [0.1, 0.5] \
  --n-simulations 10000 \
  --output adaptive_design.json

# Step 3: Compare operating characteristics
python scripts/compare.py \
  --designs fixed_baseline.json adaptive_design.json \
  --metrics power,ess,type_i_error \
  --output comparison_report.pdf

# Step 4: Generate regulatory briefing document
python scripts/regulatory_report.py \
  --design adaptive_design.json \
  --template fda_adaptive_guidance \
  --output fda_briefing.pdf
```

**Python API:**

```python
from scripts.simulator import AdaptiveSimulator
from scripts.visualizer import DesignVisualizer
from scripts.regulatory import RegulatoryReport

# Initialize
sim = AdaptiveSimulator()
viz = DesignVisualizer()
reg = RegulatoryReport()

# Define candidate designs
designs = {
    "fixed": {"type": "fixed", "n": 400},
    "gs_of": {"type": "gs", "spending": "obrien_fleming", "interims": 2},
    "gs_pocock": {"type": "gs", "spending": "pocock", "interims": 2},
    "adaptive": {"type": "ssr", "initial_n": 300, "max_n": 500}
}

# Simulate all designs
results = {}
for name, config in designs.items():
    results[name] = sim.simulate(config, n_simulations=50000)

# Visualize comparison
viz.plot_power_curves(results, output="power_comparison.png")
viz.plot_ess_comparison(results, output="ess_comparison.png")

# Generate regulatory report
report = reg.generate(
    design=results["adaptive"],
    template="fda_adaptive_trial_guidance",
    sponsor="Your Company",
    indication="Indication X"
)
report.save("regulatory_submission.pdf")
```

## Quality Checklist

**Pre-Simulation:**
- [ ] Effect size assumptions justified (literature, pilot data)
- [ ] Alpha level appropriate (0.05 two-sided or 0.025 one-sided)
- [ ] Power target realistic (usually 80-90%)
- [ ] Interim timing feasible (accrual rate, data cleaning time)

**During Simulation:**
- [ ] Sufficient simulations (≥10,000 for Type I error, ≥5,000 for power)
- [ ] Multiple seeds tested (reproducibility check)
- [ ] Various effect sizes evaluated (null, small, expected, large)
- [ ] Assumption sensitivity analyses (variance, dropout rates)

**Post-Simulation:**
- [ ] Type I error controlled (≤ nominal alpha)
- [ ] Power meets target (≥ planned power)
- [ ] Expected sample size calculated under various scenarios
- [ ] Early stopping probabilities reasonable
- [ ] Boundary values clinically interpretable

**Before Implementation:**
- [ ] **CRITICAL**: Statistical analysis plan (SAP) drafted with adaptations
- [ ] **CRITICAL**: Independent Data Safety Monitoring Board (DSMB) charter prepared
- [ ] Regulatory pre-submission meeting conducted (if required)
- [ ] Simulation code and results documented and archived
- [ ] Independent statistician reviews design

## Common Pitfalls

**Statistical Issues:**
- ❌ **Type I error inflation** → Uncontrolled alpha due to naive adaptations
  - ✅ Always use proper alpha spending functions or combination tests
  
- ❌ **Over-optimistic assumptions** → Expected effect larger than realistic
  - ✅ Simulate under range of plausible effect sizes

- ❌ **Insufficient simulations** → Unstable operating characteristics
  - ✅ Use ≥10,000 simulations for Type I error estimation

**Design Issues:**
- ❌ **Interim analysis timing unrealistic** → Data not ready when planned
  - ✅ Account for data cleaning and database lock time

- ❌ **Sample size re-estimation too aggressive** → Large Type I error inflation
  - ✅ Limit maximum sample size increase (e.g., 2× initial N)

- ❌ **Ignoring overruns** → Patients enrolled between interim and decision
  - ✅ Account for accrual during decision-making period

**Implementation Issues:**
- ❌ **Unblinded interim analysis** → Bias in trial conduct
  - ✅ Use independent statistician for interim analyses

- ❌ **Changing primary endpoint post-hoc** → Invalidates design
  - ✅ Lock primary endpoint and adaptations in protocol/SAP

- ❌ **Not following stopping rules** → Operating characteristics void
  - ✅ Make stopping rules binding; document any deviations

## References

Available in `references/` directory:

- `group_sequential_theory.md` - O'Brien-Fleming, Pocock, spending functions
- `adaptive_design_methods.md` - SSR, combination tests, conditional power
- `regulatory_guidance.md` - FDA, EMA adaptive trial guidance documents
- `simulation_best_practices.md` - Monte Carlo methodology
- `mams_designs.md` - Multi-arm multi-stage methodology
- `software_comparison.md` - Comparison with East, ADDPLAN, rpact

## Scripts

Located in `scripts/` directory:

- `main.py` - CLI interface for simulations
- `simulator.py` - Core Monte Carlo simulation engine
- `boundaries.py` - Alpha spending function calculations
- `ssr_methods.py` - Sample size re-estimation algorithms
- `mams_simulator.py` - Multi-arm multi-stage simulations
- `visualizer.py` - Boundary plots, power curves, OC comparisons
- `regulatory_report.py` - FDA/EMA briefing document generation
- `validator.py` - Type I error validation and bias assessment

## Limitations

- **Frequentist Focus**: Designed for frequentist adaptive methods; Bayesian designs require different tools
- **Normal Endpoint**: Optimized for continuous endpoints; binary/survival may need adjustments
- **Two-Arm Focus**: Multi-arm simulations available but more complex
- **No Covariate Adjustment**: Standard analyses; adaptive covariate adjustment not implemented
- **Computational Intensity**: Large simulations (100K+) may take hours
- **Regulatory Acceptance**: Designs must be vetted with regulators; tool provides simulation, not regulatory approval

---

**📊 Remember: Adaptive designs offer flexibility but require careful planning and rigorous validation. Always conduct thorough simulations and engage regulators early in the design process.**
