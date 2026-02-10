---
name: blockbuster-therapy-predictor
description: Predict which early-stage biotechnology platforms (PROTAC, mRNA, gene 
  editing, etc.) have the highest potential to become blockbuster therapies. 
  Analyzes clinical trial progression, patent landscape maturity, and venture 
  capital funding trends to generate investment and R&D prioritization scores.
allowed-tools: [Read, Write, Bash, Edit]
license: MIT
metadata:
    skill-author: AIPOCH
---

# Blockbuster Therapy Predictor

## Overview

Comprehensive analytics tool for forecasting breakthrough therapeutic technologies by integrating multi-dimensional data sources including clinical development pipelines, intellectual property landscapes, and capital market indicators.

**Key Capabilities:**
- **Multi-Source Data Integration**: Aggregates clinical trials, patents, and funding data
- **Predictive Scoring**: Calculates Blockbuster Index combining maturity, market potential, and momentum
- **Technology Landscape Mapping**: Tracks 10+ emerging therapeutic platforms
- **Investment Intelligence**: Provides data-driven R&D and investment recommendations
- **Trend Analysis**: Identifies acceleration patterns and inflection points

## When to Use

**✅ Use this skill when:**
- Making strategic R&D investment decisions in early-stage biotech
- Evaluating which therapeutic platforms to prioritize in drug discovery portfolios
- Preparing investor presentations requiring technology landscape analysis
- Assessing competitive positioning of emerging modalities
- Conducting due diligence on biotech startups or platform technologies
- Monitoring technology trends for business development opportunities

**❌ Do NOT use when:**
- Making final investment decisions without expert domain review → Use as preliminary screening only
- Analyzing individual drug candidates → Use `drug-candidate-evaluator` for single-asset assessment
- Predicting regulatory approval outcomes → Use `regulatory-pathway-analyzer`
- Estimating precise market size or revenue forecasts → Use `market-access-value` for detailed commercial analysis
- Evaluating technologies without clinical data → Model requires Phase I+ data for accuracy

**Related Skills:**
- **上游**: `competitor-trial-monitor` (clinical landscape), `patent-landscape` (IP analysis)
- **下游**: `market-access-value` (commercial assessment), `biotech-pitch-deck-narrative` (investor materials)

## Integration with Other Skills

**Upstream Skills:**
- `competitor-trial-monitor`: Track clinical trial progression for technology assessment
- `patent-landscape`: Analyze IP landscape depth and freedom-to-operate
- `funding-trend-forecaster`: Identify venture capital patterns and investor sentiment
- `clinicaltrials-gov-parser`: Extract detailed trial data for maturity scoring

**Downstream Skills:**
- `market-access-value`: Estimate commercial potential of high-scoring technologies
- `biotech-pitch-deck-narrative`: Create investor presentations for promising platforms
- `target-novelty-scorer`: Assess biological target novelty within top technologies
- `business-model-canvas`: Develop commercialization strategies for breakthrough therapies

**Complete Workflow:**
```
Competitor Trial Monitor → 
  Blockbuster Therapy Predictor (this skill) → 
    Market Access Value (commercial assessment) → 
      Biotech Pitch Deck Narrative (investor materials)
```

## Core Capabilities

### 1. Multi-Dimensional Data Collection

Automated aggregation from diverse data sources:

```python
from scripts.predictor import BlockbusterPredictor

predictor = BlockbusterPredictor()

# Collect comprehensive data
data = predictor.collect_data(
    technologies=["PROTAC", "mRNA", "CRISPR", "CAR-T"],
    data_sources=["clinical", "patent", "funding"],
    time_range="2020-2026"
)

# Returns structured dataset with:
# - Clinical: trial phases, indications, success rates
# - Patent: application trends, core IP, geography
# - Funding: rounds, amounts, investor tiers
```

**Data Sources Integrated:**
| Source | Metrics | Update Frequency |
|--------|---------|------------------|
| ClinicalTrials.gov | Phase distribution, enrollment, endpoints | Weekly |
| USPTO/EPO | Applications, grants, citations | Monthly |
| PitchBook/Crunchbase | Funding rounds, valuations, investors | Weekly |
| Evaluate Pharma | Market forecasts, competitive intel | Quarterly |

**Best Practices:**
- ✅ Always specify time range to capture recent trends
- ✅ Include all three data dimensions for accurate scoring
- ✅ Verify data completeness before running predictions
- ✅ Cache results to avoid repeated API calls

**Common Issues and Solutions:**

**Issue: Missing or incomplete data**
- Symptom: "Data collection returned 60% completeness for RNAi"
- Solution: Extend time range or check alternative data sources; flag low-confidence predictions

**Issue: Rate limiting from APIs**
- Symptom: "HTTP 429 errors from patent databases"
- Solution: Implement exponential backoff; use cached data if available

### 2. Predictive Scoring Algorithms

Multi-factor scoring models combining maturity, market, and momentum:

```python
# Calculate comprehensive scores
scores = predictor.calculate_scores(
    data=data,
    models=["maturity", "market_potential", "blockbuster_index"]
)

# Access individual scores
maturity_scores = scores.maturity
market_scores = scores.market_potential
blockbuster_index = scores.blockbuster_index

# Generate rankings
rankings = predictor.rank_technologies(scores)
```

**Scoring Models:**

**Maturity Score (0-100):**
```
Maturity = 0.4 × Clinical_Stage + 0.3 × Patent_Depth + 0.3 × Funding_Stage
```
- Clinical_Stage: Weighted average of trial phases (Phase III highest)
- Patent_Depth: Core patent families and citations
- Funding_Stage: Series C+ considered mature

**Market Potential Score (0-100):**
```
Market_Potential = 0.35 × Addressable_Market + 0.35 × Unmet_Need + 0.30 × Competition
```
- Addressable_Market: Total potential revenue (billions USD)
- Unmet_Need: Disease severity and current treatment gaps
- Competition: Number of companies and clinical assets

**Blockbuster Index (0-100):**
```
Blockbuster_Index = 0.5 × Market_Potential + 0.3 × Maturity + 0.2 × Momentum
```
- Momentum: 12-month trend in trials, patents, and funding

**Parameters:**
| Parameter | Type | Required | Description | Default |
|-----------|------|----------|-------------|---------|
| `data` | DataFrame | Yes | Collected multi-source data | None |
| `models` | list | No | Scoring models to calculate | ["all"] |
| `weights` | dict | No | Custom model weights | Standard weights |
| `momentum_window` | int | No | Months for trend calculation | 12 |

**Best Practices:**
- ✅ Compare multiple scoring models for consensus
- ✅ Adjust weights based on investment horizon (short vs long term)
- ✅ Validate scores against known blockbusters (humira, keytruda)
- ✅ Consider correlation between technology categories

**Common Issues and Solutions:**

**Issue: Counterintuitive rankings**
- Symptom: "CAR-T scores higher than mRNA despite fewer approved drugs"
- Solution: Check weight configurations; CAR-T may score high on maturity but lower on market potential

**Issue: All scores clustered in middle range**
- Symptom: "All technologies score 45-55"
- Solution: Normalize scores within technology class; some modalities inherently earlier stage

### 3. Technology Landscape Visualization

Generate interactive reports and visualizations:

```python
# Create comprehensive report
report = predictor.generate_report(
    scores=scores,
    include_charts=True,
    format="html"  # Options: html, pdf, json
)

# Save visualization
predictor.save_visualizations(
    scores=scores,
    output_dir="figures/",
    chart_types=["scatter", "heatmap", "timeline"]
)
```

**Visualization Types:**
| Chart | Purpose | Insights |
|-------|---------|----------|
| **Scatter Plot** | Maturity vs Market Potential | Identify high-potential quadrants |
| **Heatmap** | Technology × Indication | Find whitespace opportunities |
| **Timeline** | Cumulative trials/patents over time | Spot inflection points |
| **Treemap** | Market size by technology | Prioritize largest opportunities |
| **Network** | Company-technology relationships | Map competitive landscape |

**Best Practices:**
- ✅ Use scatter plot as primary decision-support visualization
- ✅ Include trend arrows showing momentum direction
- ✅ Color-code by therapeutic area for pattern recognition
- ✅ Export interactive HTML for stakeholder presentations

**Common Issues and Solutions:**

**Issue: Overcrowded visualizations**
- Symptom: "Too many technologies, chart unreadable"
- Solution: Filter to top 10 by Blockbuster Index; create separate charts by therapeutic area

**Issue: Missing legend or labels**
- Symptom: "Stakeholders can't interpret bubble sizes"
- Solution: Always include comprehensive legends; annotate key data points

### 4. Investment Recommendation Engine

Translate scores into actionable investment guidance:

```python
# Generate recommendations
recommendations = predictor.generate_recommendations(
    scores=scores,
    investor_type="pharma_corp",  # Options: vc, pharma_corp, biotech
    risk_tolerance="moderate",    # Options: conservative, moderate, aggressive
    time_horizon="5_year"         # Options: 3_year, 5_year, 10_year
)

# Output includes:
# - Buy/Watch/Avoid ratings per technology
# - Portfolio allocation suggestions
# - Risk-adjusted expected returns
```

**Recommendation Categories:**
| Blockbuster Index | Rating | Interpretation | Action |
|-------------------|--------|----------------|--------|
| ≥ 80 | **Strong Buy** | High probability blockbuster | Prioritize R&D investment |
| 60-79 | **Buy** | Promising potential | Active monitoring and early partnerships |
| 40-59 | **Watch** | Moderate potential | Monitor milestones; reassess in 6-12 months |
| < 40 | **Avoid** | Low probability | Minimal investment; consider divestment |

**Best Practices:**
- ✅ Consider investor type (VC vs pharma have different risk profiles)
- ✅ Adjust for therapeutic area expertise
- ✅ Combine quantitative scores with qualitative domain knowledge
- ✅ Re-evaluate quarterly as new data emerges

**Common Issues and Solutions:**

**Issue: Conflicting recommendations across investor types**
- Symptom: "VC recommends Buy, pharma recommends Watch for same technology"
- Solution: Normal; reflects different risk-return profiles. Document rationale for each audience.

**Issue: Overly conservative rankings**
- Symptom: "No technologies score above 70"
- Solution: Check benchmark calibration; early-stage platforms naturally score lower on maturity

### 5. Scenario Analysis and Sensitivity Testing

Test robustness of predictions under different assumptions:

```python
# Run sensitivity analysis
sensitivity = predictor.sensitivity_analysis(
    base_scores=scores,
    variables=["clinical_success_rate", "patent_expiry", "funding_collapse"],
    scenarios=["optimistic", "pessimistic", "baseline"]
)

# Monte Carlo simulation
monte_carlo = predictor.monte_carlo_simulation(
    scores=scores,
    n_iterations=10000,
    variables={
        "clinical_success_rate": (0.1, 0.3),  # Uniform distribution
        "market_growth": (0.05, 0.15)         # Range
    }
)
```

**Scenario Types:**
| Scenario | Assumptions | Use Case |
|----------|-------------|----------|
| **Optimistic** | High clinical success, strong IP protection, abundant funding | Best-case valuation |
| **Baseline** | Historical averages | Primary recommendation |
| **Pessimistic** | Low success rates, patent challenges, funding drought | Risk assessment |

**Best Practices:**
- ✅ Always run at least 3 scenarios (optimistic/pessimistic/baseline)
- ✅ Focus on technologies that remain top-ranked across all scenarios
- ✅ Use Monte Carlo for technologies with high uncertainty
- ✅ Document scenario assumptions for audit trail

**Common Issues and Solutions:**

**Issue: Extreme scenario results**
- Symptom: "Pessimistic scenario shows all technologies failing"
- Solution: Check scenario parameter bounds; ensure realistic ranges

**Issue: Monte Carlo convergence issues**
- Symptom: "Results vary significantly across runs"
- Solution: Increase n_iterations to 50,000+; check random seed

### 6. Competitive Intelligence and Benchmarking

Compare against historical blockbusters and competitors:

```python
# Benchmark against known blockbusters
benchmark = predictor.benchmark_analysis(
    scores=scores,
    reference_blockbusters=["Humira", "Keytruda", "Revlimid"],
    metrics=["trajectory", "milestones", "funding_patterns"]
)

# Compare to competitors
competitive_intel = predictor.competitive_analysis(
    technology="CRISPR",
    companies=["Editas", "Intellia", "CRISPR Therapeutics"],
    metrics=["pipeline_depth", "cash_position", "partnerships"]
)
```

**Benchmark Dimensions:**
| Dimension | Metrics | Insight |
|-----------|---------|---------|
| **Development Trajectory** | Trials over time | Is technology on blockbuster path? |
| **Funding Velocity** | Round sizes, frequency | Investor confidence indicator |
| **Partnership Quality** | Pharma deals, terms | Validation by Big Pharma |
| **IP Position** | Patent citations, breadth | Defensibility and freedom-to-operate |

**Best Practices:**
- ✅ Benchmark against 3-5 historical blockbusters with similar mechanisms
- ✅ Compare development timelines (e.g., mRNA COVID vaccines were unusually fast)
- ✅ Analyze failed technologies to identify red flags
- ✅ Track competitor pipeline depth and cash runway

**Common Issues and Solutions:**

**Issue: Historical blockbusters had different market conditions**
- Symptom: "No technology matches Humira's trajectory"
- Solution: Adjust for market context; focus on trajectory shape, not absolute numbers

**Issue: Missing competitor data**
- Symptom: "Private company financials unavailable"
- Solution: Use estimates from PitchBook; focus on public data (trials, patents)

## Complete Workflow Example

**From technology landscape to investment decision:**

```bash
# Step 1: Collect comprehensive data
python scripts/main.py --collect \
  --technologies PROTAC,mRNA,CRISPR,CAR-T,ADC \
  --sources clinical,patent,funding \
  --output data/tech_landscape.pkl

# Step 2: Calculate blockbuster scores
python scripts/main.py --score \
  --input data/tech_landscape.pkl \
  --models all \
  --output results/scores.json

# Step 3: Generate visualizations
python scripts/main.py --visualize \
  --scores results/scores.json \
  --charts scatter,heatmap,timeline \
  --output figures/

# Step 4: Get investment recommendations
python scripts/main.py --recommend \
  --scores results/scores.json \
  --investor-type pharma_corp \
  --risk-tolerance moderate \
  --output report/recommendations.md

# Step 5: Run scenario analysis
python scripts/main.py --scenario \
  --scores results/scores.json \
  --scenarios optimistic,pessimistic,baseline \
  --output report/sensitivity_analysis.pdf
```

**Python API Usage:**

```python
from scripts.predictor import BlockbusterPredictor
from scripts.visualizer import LandscapeVisualizer
from scripts.recommender import InvestmentRecommender

# Initialize
predictor = BlockbusterPredictor()
viz = LandscapeVisualizer()
rec = InvestmentRecommender()

# Complete analysis pipeline
data = predictor.collect_data(
    technologies=["PROTAC", "mRNA", "CRISPR"],
    time_range="2020-2026"
)

scores = predictor.calculate_scores(data)
rankings = predictor.rank_technologies(scores)

# Generate outputs
viz.create_scatter_plot(scores, output="figures/landscape.png")
viz.create_timeline(scores, output="figures/trends.png")

recommendations = rec.generate(scores, investor_type="vc")
recommendations.save("report/investment_recommendations.md")

# Benchmark analysis
benchmark = predictor.benchmark_analysis(
    scores=scores,
    reference_blockbusters=["Keytruda", "Humira"]
)
```

**Expected Output Files:**
```
output/
├── data/
│   └── tech_landscape.pkl          # Raw collected data
├── results/
│   ├── scores.json                 # Comprehensive scoring results
│   └── rankings.csv                # Technology rankings
├── figures/
│   ├── landscape_scatter.png       # Maturity vs Market Potential
│   ├── technology_heatmap.png      # Tech × Indication matrix
│   └── funding_timeline.png        # Cumulative funding over time
└── report/
    ├── investment_recommendations.md  # Actionable guidance
    └── sensitivity_analysis.pdf       # Scenario analysis
```

## Quality Checklist

**Pre-Analysis Checks:**
- [ ] Technology list is comprehensive (no major modalities missing)
- [ ] Data sources are accessible and up-to-date
- [ ] Historical blockbuster benchmarks selected appropriately
- [ ] Time range captures recent trends (minimum 2 years)

**During Data Collection:**
- [ ] All three data dimensions collected (clinical, patent, funding)
- [ ] Data completeness > 80% for each technology
- [ ] Outliers identified and validated (not data errors)
- [ ] Missing data flagged with confidence intervals

**During Scoring:**
- [ ] Model weights reflect investment priorities
- [ ] Scores normalized appropriately within technology classes
- [ ] Validation against known blockbusters shows expected patterns
- [ ] Correlation analysis reveals no multicollinearity issues

**Post-Analysis Verification:**
- [ ] Top-ranked technologies pass sanity check (domain expert review)
- [ ] Scenario analysis shows robustness (top techs remain top across scenarios)
- [ ] Sensitivity analysis identifies key uncertainty drivers
- [ ] Visualizations are clear and actionable

**Before Recommendation:**
- [ ] **CRITICAL**: Recommendations reviewed by domain expert (not purely algorithmic)
- [ ] Investor type and risk tolerance appropriately configured
- [ ] Limitations and uncertainties clearly documented
- [ ] Benchmark comparison shows technology on viable path

## Common Pitfalls

**Data Quality Issues:**
- ❌ **Incomplete data** (only clinical, no patent/funding) → Scores biased toward clinical-stage companies
  - ✅ Collect all three data dimensions; flag incomplete analyses
  
- ❌ **Stale data** (analysis based on 2022 data in 2026) → Miss recent inflection points
  - ✅ Use data within last 6 months; set up automated refresh

- ❌ **Survivorship bias** (only analyzing successful trials) → Overestimate success rates
  - ✅ Include failed trials and terminated programs

**Model Issues:**
- ❌ **Overfitting to historical blockbusters** → Miss novel mechanisms
  - ✅ Validate model on holdout set; adjust for novelty

- ❌ **Ignoring correlation** (mRNA and COVID vaccine success correlated) → Double-count pandemic effect
  - ✅ Include covariance analysis; adjust for macro trends

- ❌ **Static weights** (same weights in 2026 as 2020) → Miss evolving market conditions
  - ✅ Review and adjust weights annually based on market feedback

**Interpretation Issues:**
- ❌ **Overconfidence in scores** (treating 75 as precise prediction) → Ignore uncertainty
  - ✅ Always report confidence intervals; use probabilistic language

- ❌ **Ignoring qualitative factors** (team quality, manufacturing complexity) → Miss key risks
  - ✅ Combine quantitative scores with expert qualitative assessment

- ❌ **Short-term focus** (prioritizing Phase III over Phase I breakthroughs) → Miss next-generation platforms
  - ✅ Include 10-year horizon scenarios; balance maturity vs innovation

**Recommendation Issues:**
- ❌ **One-size-fits-all** (same recommendation for VC and pharma) → Ignore different risk profiles
  - ✅ Customize recommendations by investor type

- ❌ **Ignoring portfolio context** (recommending mRNA when already 50% allocated) → Concentration risk
  - ✅ Consider current portfolio; recommend diversification

**Communication Issues:**
- ❌ **Jargon-heavy reports** → Stakeholders don't understand
  - ✅ Use plain language; define technical terms

- ❌ **Missing visualizations** → Tables don't tell story
  - ✅ Always include scatter plots and trend charts

## Troubleshooting

**Problem: Low data completeness**
- Symptoms: "Only 60% of technologies have complete data"
- Causes: New modalities without clinical trials; private companies with limited disclosure
- Solutions:
  - Extend data collection to secondary sources (conference abstracts, SEC filings)
  - Use imputation for missing values with uncertainty quantification
  - Flag low-confidence predictions; require expert override

**Problem: Counterintuitive rankings**
- Symptoms: "Well-established technology ranks below unproven platform"
- Causes: Weighting emphasizes momentum over maturity; data errors
- Solutions:
  - Review weight configurations; consider multiple scoring models
  - Validate raw data (e.g., check if Phase III trial actually exists)
  - Include domain expert review for top 10 rankings

**Problem: All scores in narrow range**
- Symptoms: "All technologies score 45-55; no clear leaders"
- Causes: Normalization across heterogeneous modalities; conservative assumptions
- Solutions:
  - Normalize within technology class (gene therapy vs gene therapy)
  - Use percentile rankings instead of absolute scores
  - Check model calibration against historical blockbusters

**Problem: Model outputs change dramatically on re-run**
- Symptoms: "Rankings different every time"
- Causes: Random seed issues; non-deterministic API responses
- Solutions:
  - Set random seeds for reproducibility
  - Cache API responses; version control data snapshots
  - Document data provenance and timestamp

**Problem: API rate limiting**
- Symptoms: "HTTP 429 errors from ClinicalTrials.gov"
- Causes: Too many requests in short time
- Solutions:
  - Implement exponential backoff with jitter
  - Use API keys for higher rate limits
  - Schedule data collection during off-peak hours
  - Cache data and refresh incrementally

**Problem: Visualization rendering errors**
- Symptoms: "Charts don't display; missing fonts"
- Causes: Missing matplotlib backends; font dependencies
- Solutions:
  - Install required system fonts (Arial, Times New Roman)
  - Use Agg backend for headless servers: `matplotlib.use('Agg')`
  - Export to multiple formats (PNG fallback for PDF issues)

**Problem: Recommendations not actionable**
- Symptoms: "Stakeholders ask 'so what?' after reading report"
- Causes: Too academic; missing clear next steps
- Solutions:
  - Include explicit "Recommended Actions" section
  - Quantify investment amounts and timelines
  - Provide alternative options (if-then scenarios)
  - Schedule follow-up meeting to discuss implementation

## References

Available in `references/` directory:

- `scoring_methodology.md` - Detailed mathematical models and algorithms
- `data_sources_guide.md` - API documentation and access credentials
- `historical_blockbusters.md` - Case studies of successful drugs and their trajectories
- `technology_profiles.md` - Deep dives into each therapeutic platform
- `competitive_landscape.md` - Company and technology mapping
- `investor_frameworks.md` - Best practices for biotech investment analysis

## Scripts

Located in `scripts/` directory:

- `main.py` - CLI interface for full analysis pipeline
- `predictor.py` - Core predictive scoring algorithms
- `data_collector.py` - Multi-source data aggregation
- `visualizer.py` - Report and chart generation
- `recommender.py` - Investment recommendation engine
- `benchmark.py` - Historical blockbuster comparison
- `sensitivity.py` - Scenario and Monte Carlo analysis
- `competitive_intel.py` - Company and pipeline analysis
- `utils.py` - Helper functions and data validation

## Performance and Resources

**Typical Analysis Time:**
- Small dataset (5-10 technologies): 5-10 minutes
- Medium dataset (10-20 technologies): 15-30 minutes
- Full landscape (all 10+ technologies): 45-90 minutes

**System Requirements:**
- **RAM**: 4 GB minimum, 8 GB recommended for full landscape
- **Storage**: ~1 GB for models and cached data
- **Network**: Stable connection for API calls (ClinicalTrials.gov, USPTO, etc.)
- **CPU**: Multi-core speeds up data collection

**API Rate Limits:**
- ClinicalTrials.gov: 100 requests/minute (API key recommended)
- USPTO: 1000 requests/day (API key required)
- PitchBook: Varies by subscription tier
- Plan data collection accordingly to avoid throttling

## Limitations

- **Prediction Accuracy**: Model accuracy ~70% for 5-year horizon; improves with shorter timeframes
- **Data Availability**: Private company data limited; public companies more accurate
- **Novel Mechanisms**: Historical patterns may not predict breakthrough innovations
- **Regulatory Changes**: Cannot anticipate major policy shifts (e.g., FDA breakthrough therapy designation effects)
- **Black Swan Events**: Pandemics, regulatory scandals not captured in historical data
- **Therapeutic Area Expertise**: Model provides quantitative scores; domain expertise required for interpretation
- **Currency**: Market conditions change; quarterly refresh recommended
- **Conflicts of Interest**: Model developers may have biases toward certain technologies

## Version History

- **v1.0.0** (Current): Initial release with 10 technology pathways, 3-dimension scoring, visualization suite
- Planned: Expansion to 20+ modalities; real-time data streaming; AI-driven pattern recognition

---

**⚠️ DISCLAIMER: This tool provides quantitative analysis for decision support only. All investment and R&D decisions should incorporate qualitative domain expertise, regulatory consultation, and comprehensive due diligence. Past performance of historical blockbusters does not guarantee future success of emerging technologies.**
