---
name: crispr-screen-analyzer
description: Process CRISPR screening data to identify essential genes and hit candidates. Performs quality control, statistical analysis (RRA), and hit calling for pooled CRISPR screens including viability screens and drug resistance/sensitivity studies.
allowed-tools: [Read, Write, Bash, Edit, Grep]
license: MIT
metadata:
  skill-author: AIPOCH
---

# CRISPR Screen Analyzer

Analyze pooled CRISPR screening data to identify essential genes, drug resistance/sensitivity candidates, and screen quality metrics. Supports Robust Rank Aggregation (RRA) analysis, quality control assessment, and hit identification for functional genomics studies.

**Key Capabilities:**
- **Quality Control Assessment**: Calculate Gini index, read depth, and dropout metrics to evaluate screen quality
- **Log Fold Change Calculation**: Compute sgRNA-level fold changes between treatment and control conditions
- **Statistical Analysis**: Perform Robust Rank Aggregation (RRA) to identify significantly enriched or depleted sgRNAs
- **Hit Identification**: Apply FDR and fold change thresholds to identify candidate genes
- **Multi-Sample Support**: Process multiple replicates and treatment conditions simultaneously

---

## When to Use

**✅ Use this skill when:**
- Analyzing **genome-wide viability screens** to identify essential genes required for cell survival
- Performing **drug resistance screens** to find genes whose knockout confers resistance
- Conducting **drug sensitivity screens** to identify synthetic lethal interactions
- Performing **quality control assessment** of CRISPR screen data before downstream analysis
- **Comparing multiple treatment conditions** (e.g., drug vs DMSO, hypoxia vs normoxia)
- **Validating screen quality** before publication or further experimental validation
- **Generating hit lists** for secondary screens or validation experiments

**❌ Do NOT use when:**
- Analyzing **single-cell CRISPR data** (Perturb-seq, CROP-seq) → Use specialized single-cell analysis tools
- Working with **arrayed CRISPR screens** (well-by-well format) → Use standard differential expression analysis
- Performing **CRISPR activation (CRISPRa) or interference (CRISPRi)** screens → May need adjusted normalization
- Requiring **Bayesian or MAGeCK** statistical analysis → This tool uses RRA; use MAGeCK for alternative algorithms
- Analyzing **small custom libraries** (<1000 sgRNAs) → Statistical power may be insufficient
- **Time-course CRISPR screens** → Requires specialized trajectory analysis methods

**Related Skills:**
- **上游 (Upstream)**: `crispr-grna-designer`, `fastqc-report-interpreter`
- **下游 (Downstream)**: `go-kegg-enrichment`, `pathway-visualization`, `hit-validation-planner`

---

## Integration with Other Skills

**Upstream Skills:**
- `crispr-grna-designer`: Design sgRNA libraries before screening; validate library composition
- `fastqc-report-interpreter`: Assess sequencing quality before CRISPR screen analysis
- `alignment-quality-checker`: Verify sgRNA alignment rates and mapping quality

**Downstream Skills:**
- `go-kegg-enrichment`: Perform pathway enrichment on identified hit genes
- `pathway-visualization`: Visualize hits in pathway contexts
- `hit-validation-planner`: Design follow-up experiments for candidate genes
- `gene-essentiality-predictor`: Compare screen results with known essential gene databases

**Complete Workflow:**
```
Library Design (crispr-grna-designer) → Transduction → Sequencing → fastqc-report-interpreter → crispr-screen-analyzer → go-kegg-enrichment → Hit Validation
```

---

## Core Capabilities

### 1. Quality Control Metrics Calculation

Assess CRISPR screen quality using established metrics including Gini index, read depth, and sgRNA dropout rates.

```python
from scripts.main import CRISPRScreenAnalyzer

# Initialize analyzer with count matrix and sample annotations
analyzer = CRISPRScreenAnalyzer(
    counts_file="sgrna_counts.txt",
    samplesheet="samples.csv"
)

# Calculate QC metrics
qc_results = analyzer.qc_metrics()

# Review key metrics
print("Quality Control Metrics:")
print(f"Total reads per sample:")
for sample, reads in qc_results['total_reads'].items():
    print(f"  {sample}: {reads:,} reads")

print(f"\nGini index (library representation):")
for sample, gini in qc_results['gini_index'].items():
    status = "✅ Good" if gini < 0.3 else "⚠️  Check" if gini < 0.4 else "❌ Poor"
    print(f"  {sample}: {gini:.3f} {status}")

print(f"\nZero-count sgRNAs (potential dropout):")
for sample, zeros in qc_results['zero_count_sgrnas'].items():
    pct = (zeros / len(analyzer.counts)) * 100
    print(f"  {sample}: {zeros} ({pct:.1f}%)")
```

**QC Metrics Explained:**

| Metric | Target Range | Interpretation |
|--------|--------------|----------------|
| **Gini Index** | <0.3 | Measures library evenness; lower = more uniform |
| **Total Reads** | >10M per sample | Sufficient depth for statistical power |
| **Zero-count sgRNAs** | <5% | Acceptable dropout; higher indicates library loss |
| **Read Distribution** | Log-normal | Should follow expected distribution |

**Best Practices:**
- ✅ **Check Gini index first**: Values >0.4 indicate potential library bias or bottleneck
- ✅ **Compare replicates**: QC metrics should be consistent across replicates
- ✅ **Assess time points**: Later time points typically show higher dropout
- ✅ **Validate early**: Poor QC may require screen repetition

**Common Issues and Solutions:**

**Issue: High Gini index (>0.4)**
- Symptom: Uneven sgRNA representation suggesting library bottleneck
- Solution: Check MOI (multiplicity of infection); verify puromycin selection; consider repeating screen

**Issue: Excessive zero-count sgRNAs (>10%)**
- Symptom: Many sgRNAs not detected in final samples
- Causes: Low sequencing depth, library degradation, or strong selection
- Solution: Increase sequencing depth; verify library quality at transduction

### 2. Log Fold Change Calculation

Calculate log2 fold changes between treatment and control conditions to identify enriched or depleted sgRNAs.

```python
from scripts.main import CRISPRScreenAnalyzer

analyzer = CRISPRScreenAnalyzer("counts.txt", "samples.csv")

# Define sample groups
control_samples = ["Control_1", "Control_2", "Control_3"]
treatment_samples = ["Drug_1", "Drug_2", "Drug_3"]

# Calculate log fold changes
lfc = analyzer.calculate_lfc(control_samples, treatment_samples)

# Analyze distribution
print("Log Fold Change Statistics:")
print(f"  Mean: {lfc.mean():.3f}")
print(f"  Std:  {lfc.std():.3f}")
print(f"  Max:  {lfc.max():.3f}")
print(f"  Min:  {lfc.min():.3f}")

# Identify extreme changes
strong_depletion = lfc[lfc < -2]  # Strong negative selection
strong_enrichment = lfc[lfc > 2]   # Strong positive selection

print(f"\nStrongly depleted sgRNAs: {len(strong_depletion)}")
print(f"Strongly enriched sgRNAs: {len(strong_enrichment)}")
```

**LFC Calculation:**
```
lfc = log2((treatment_mean + 1) / (control_mean + 1))
```

**Interpretation:**

| LFC Range | Interpretation | Biological Meaning |
|-----------|---------------|-------------------|
| **LFC < -2** | Strong depletion | Essential gene or drug sensitivity |
| **LFC -2 to -1** | Moderate depletion | Moderate effect |
| **LFC -1 to 1** | No change | No significant effect |
| **LFC 1 to 2** | Moderate enrichment | Moderate resistance |
| **LFC > 2** | Strong enrichment | Resistance gene or suppressor |

**Best Practices:**
- ✅ **Use pseudocount of 1** to avoid log(0) issues
- ✅ **Average replicates** to reduce technical variance
- ✅ **Visualize distribution** to identify batch effects or outliers
- ✅ **Check positive controls** (known essential genes should have negative LFC)

**Common Issues and Solutions:**

**Issue: Skewed LFC distribution**
- Symptom: Mean LFC significantly different from 0
- Causes: Library size differences, batch effects, or strong selection
- Solution: Apply TMM or DESeq2 normalization; check for batch effects

**Issue: Extreme outliers**
- Symptom: Few sgRNAs with very large LFC values
- Solution: Winsorize extreme values; verify these are not technical artifacts

### 3. Robust Rank Aggregation (RRA) Statistical Analysis

Perform statistical analysis to identify significantly enriched or depleted sgRNAs using z-score and FDR correction.

```python
from scripts.main import CRISPRScreenAnalyzer

analyzer = CRISPRScreenAnalyzer("counts.txt", "samples.csv")

# Calculate LFC first
lfc = analyzer.calculate_lfc(
    control_samples=["Ctrl_1", "Ctrl_2"],
    treatment_samples=["Treat_1", "Treat_2"]
)

# Perform RRA analysis
results = analyzer.rra_analysis(lfc, fdr_threshold=0.05)

# Review top hits
print("Top 10 Most Significant sgRNAs:")
top_hits = results.nsmallest(10, 'fdr')
print(top_hits[['sgrna', 'lfc', 'pvalue', 'fdr']].to_string(index=False))

# Summary statistics
print(f"\nTotal sgRNAs tested: {len(results)}")
print(f"Significant at FDR < 0.05: {sum(results['fdr'] < 0.05)}")
print(f"Significant depletions: {sum((results['fdr'] < 0.05) & (results['lfc'] < 0))}")
print(f"Significant enrichments: {sum((results['fdr'] < 0.05) & (results['lfc'] > 0))}")
```

**RRA Analysis Steps:**

1. **Z-score calculation**: `z = (lfc - mean) / std`
2. **P-value calculation**: Two-tailed normal test
3. **FDR correction**: Benjamini-Hochberg procedure

**Statistical Output:**

| Column | Description | Usage |
|--------|-------------|-------|
| `sgrna` | sgRNA identifier | Mapping to genes |
| `lfc` | Log fold change | Effect size |
| `pvalue` | Raw p-value | Statistical significance |
| `fdr` | Adjusted p-value (FDR) | Multiple testing correction |

**Best Practices:**
- ✅ **Use FDR < 0.05** as standard significance threshold
- ✅ **Consider FDR < 0.01** for high-confidence hits
- ✅ **Combine p-value and LFC** for hit prioritization
- ✅ **Validate top hits** experimentally before publication

**Common Issues and Solutions:**

**Issue: No significant hits despite visible effects**
- Symptom: Biological effects present but no FDR-significant results
- Causes: High variance, insufficient replicates, or weak effects
- Solution: Increase replicate number; use more permissive FDR threshold; use gene-level aggregation

**Issue: Too many significant hits**
- Symptom: Hundreds or thousands of FDR-significant sgRNAs
- Causes: Low variance, strong selection, or batch effects
- Solution: Apply more stringent FDR threshold; add LFC cutoff; filter by effect size

### 4. Hit Identification with Thresholds

Apply statistical and biological thresholds to identify candidate genes for follow-up validation.

```python
from scripts.main import CRISPRScreenAnalyzer

analyzer = CRISPRScreenAnalyzer("counts.txt", "samples.csv")
lfc = analyzer.calculate_lfc(["Ctrl_1", "Ctrl_2"], ["Treat_1", "Treat_2"])
results = analyzer.rra_analysis(lfc)

# Identify hits with multiple thresholds
threshold_configs = [
    {"fdr": 0.05, "lfc": 1.0, "name": "Standard"},
    {"fdr": 0.01, "lfc": 1.5, "name": "Stringent"},
    {"fdr": 0.1, "lfc": 0.5, "name": "Permissive"}
]

for config in threshold_configs:
    hits = analyzer.identify_hits(
        results, 
        fdr_threshold=config['fdr'],
        lfc_threshold=config['lfc']
    )
    
    depletions = hits[hits['lfc'] < 0]
    enrichments = hits[hits['lfc'] > 0]
    
    print(f"\n{config['name']} (FDR<{config['fdr']}, |LFC|>{config['lfc']}):")
    print(f"  Total hits: {len(hits)}")
    print(f"  Depletions: {len(depletions)}")
    print(f"  Enrichments: {len(enrichments)}")

# Save hits for downstream analysis
standard_hits = analyzer.identify_hits(results, fdr_threshold=0.05, lfc_threshold=1.0)
standard_hits.to_csv("hits_standard.csv", index=False)
```

**Hit Classification:**

| Category | Criteria | Biological Interpretation |
|----------|----------|---------------------------|
| **Essential** | FDR<0.05, LFC<-1 | Required for cell viability |
| **Drug Sensitive** | FDR<0.05, LFC<-1 | Synthetic lethal with treatment |
| **Drug Resistant** | FDR<0.05, LFC>1 | Confers resistance to treatment |
| **Suppressor** | FDR<0.05, LFC>1 | Suppresses phenotype of interest |

**Best Practices:**
- ✅ **Use consistent thresholds** across related screens for comparability
- ✅ **Require multiple sgRNAs** per gene for confidence (≥2 recommended)
- ✅ **Validate with orthogonal methods** (siRNA, rescue experiments)
- ✅ **Compare with known essential genes** as positive controls

**Common Issues and Solutions:**

**Issue: Single sgRNA hits**
- Symptom: Only one sgRNA per gene significant
- Solution: Require ≥2 significant sgRNAs per gene; check for off-target effects

**Issue: Off-target effects dominating**
- Symptom: Known essential genes not identified; unexpected hits prominent
- Solution: Use second-generation libraries with improved specificity; validate with rescue

### 5. Gene-Level Aggregation

Aggregate sgRNA-level results to gene-level statistics for biological interpretation.

```python
import pandas as pd
from scripts.main import CRISPRScreenAnalyzer

analyzer = CRISPRScreenAnalyzer("counts.txt", "samples.csv")
lfc = analyzer.calculate_lfc(["Ctrl_1", "Ctrl_2"], ["Treat_1", "Treat_2"])
results = analyzer.rra_analysis(lfc)

# Add gene annotations (example mapping)
sgrna_to_gene = pd.read_csv("library_annotation.csv")  # sgRNA, Gene columns
results_with_gene = results.merge(sgrna_to_gene, on='sgrna')

# Aggregate to gene level
gene_results = results_with_gene.groupby('Gene').agg({
    'lfc': 'mean',           # Average LFC across sgRNAs
    'pvalue': 'min',         # Best p-value
    'fdr': 'min',            # Best FDR
    'sgrna': 'count'         # Number of sgRNAs
}).rename(columns={'sgrna': 'sgrna_count'})

# Filter genes with multiple sgRNAs
gene_results = gene_results[gene_results['sgrna_count'] >= 2]

# Identify gene-level hits
gene_hits = gene_results[
    (gene_results['fdr'] < 0.05) & 
    (abs(gene_results['lfc']) > 1.0)
]

print(f"Gene-level hits: {len(gene_hits)}")
print("\nTop 10 hits:")
print(gene_hits.nsmallest(10, 'fdr')[['lfc', 'pvalue', 'fdr', 'sgrna_count']])
```

**Gene Aggregation Methods:**

| Method | Description | Best For |
|--------|-------------|----------|
| **Mean LFC** | Average across sgRNAs | General hit calling |
| **Best FDR** | Most significant sgRNA | Conservative approach |
| **Second-best** | Second most significant | Reduces outlier effects |
| **STARS/RRA** | Rank-based aggregation | Standard CRISPR analysis |

**Best Practices:**
- ✅ **Require ≥3 sgRNAs per gene** for reliable gene-level calling
- ✅ **Use mean LFC** for primary analysis; best FDR for validation
- ✅ **Check sgRNA concordance** - all should show same direction
- ✅ **Remove genes with conflicting sgRNAs** from hit list

**Common Issues and Solutions:**

**Issue: Discordant sgRNAs for same gene**
- Symptom: Some sgRNAs positive, others negative for same gene
- Causes: Off-target effects, library errors, or complex biology
- Solution: Exclude genes with discordant sgRNAs; investigate specific cases

### 6. Multi-Condition Comparison

Compare CRISPR screen results across multiple treatment conditions or time points.

```python
from scripts.main import CRISPRScreenAnalyzer

analyzer = CRISPRScreenAnalyzer("counts.txt", "samples.csv")

# Define multiple comparisons
comparisons = {
    "Drug_A": {
        "control": ["DMSO_1", "DMSO_2"],
        "treatment": ["DrugA_1", "DrugA_2"]
    },
    "Drug_B": {
        "control": ["DMSO_1", "DMSO_2"], 
        "treatment": ["DrugB_1", "DrugB_2"]
    },
    "Combination": {
        "control": ["DMSO_1", "DMSO_2"],
        "treatment": ["Combo_1", "Combo_2"]
    }
}

# Analyze all conditions
all_results = {}
for comp_name, samples in comparisons.items():
    lfc = analyzer.calculate_lfc(samples['control'], samples['treatment'])
    results = analyzer.rra_analysis(lfc)
    hits = analyzer.identify_hits(results)
    
    all_results[comp_name] = {
        'lfc': lfc,
        'results': results,
        'hits': hits
    }
    
    print(f"{comp_name}: {len(hits)} hits")

# Find common hits across conditions
common_hits = set(all_results['Drug_A']['hits'].index)
for comp in ['Drug_B', 'Combination']:
    common_hits &= set(all_results[comp]['hits'].index)

print(f"\nCommon hits across all conditions: {len(common_hits)}")

# Compare LFC correlations between conditions
import matplotlib.pyplot as plt

lfc_drugA = all_results['Drug_A']['lfc']
lfc_drugB = all_results['Drug_B']['lfc']

correlation = lfc_drugA.corr(lfc_drugB)
print(f"\nCorrelation between Drug A and Drug B: {correlation:.3f}")
```

**Multi-Condition Analysis:**

| Comparison Type | Question Addressed | Interpretation |
|----------------|-------------------|----------------|
| **Drug vs Control** | What genes mediate drug response? | Resistance/sensitivity mechanisms |
| **Condition A vs B** | Differential genetic dependencies | Context-specific essentiality |
| **Time-course** | How does genetic dependency change? | Temporal dynamics |
| **Cell line comparison** | Cell-type specific dependencies | Lineage-specific vulnerabilities |

**Best Practices:**
- ✅ **Use same control** across multiple treatments for comparability
- ✅ **Check correlation** between replicates and conditions
- ✅ **Look for condition-specific hits** for mechanism insights
- ✅ **Validate common hits** as robust findings

**Common Issues and Solutions:**

**Issue: High variability between replicates**
- Symptom: Low correlation between replicates of same condition
- Solution: Increase replicate number; check for technical batch effects

---

## Complete Workflow Example

**From count matrix to hit identification:**

```bash
# Step 1: Run QC assessment
python scripts/main.py --counts sgrna_counts.txt --samples samples.csv --output qc_results

# Step 2: Perform differential analysis
python scripts/main.py \
  --counts sgrna_counts.txt \
  --samples samples.csv \
  --control "Ctrl_1,Ctrl_2,Ctrl_3" \
  --treatment "Drug_1,Drug_2,Drug_3" \
  --output drug_screen \
  --fdr 0.05

# Step 3: Review results
cat drug_screen_sgrna_results.csv | head -20
```

**Python API Usage:**

```python
from scripts.main import CRISPRScreenAnalyzer
import pandas as pd

def analyze_crispr_screen(
    counts_file: str,
    samplesheet: str,
    control_samples: list,
    treatment_samples: list,
    output_prefix: str,
    fdr_threshold: float = 0.05,
    lfc_threshold: float = 1.0
) -> dict:
    """
    Complete CRISPR screen analysis workflow.
    """
    # Initialize analyzer
    analyzer = CRISPRScreenAnalyzer(counts_file, samplesheet)
    
    print(f"Loaded {analyzer.counts.shape[0]} sgRNAs x {analyzer.counts.shape[1]} samples")
    
    # Quality control
    print("\n1. Quality Control Assessment...")
    qc = analyzer.qc_metrics()
    
    # Check QC status
    qc_pass = all(gini < 0.4 for gini in qc['gini_index'].values())
    if not qc_pass:
        print("⚠️  Warning: High Gini index detected - check library representation")
    
    # Calculate fold changes
    print("\n2. Calculating log fold changes...")
    lfc = analyzer.calculate_lfc(control_samples, treatment_samples)
    
    # Statistical analysis
    print("\n3. Running RRA analysis...")
    results = analyzer.rra_analysis(lfc, fdr_threshold)
    
    # Identify hits
    print("\n4. Identifying significant hits...")
    hits = analyzer.identify_hits(results, fdr_threshold, lfc_threshold)
    
    # Categorize hits
    depletions = hits[hits['lfc'] < 0]
    enrichments = hits[hits['lfc'] > 0]
    
    # Save results
    results.to_csv(f"{output_prefix}_sgrna_results.csv", index=False)
    hits.to_csv(f"{output_prefix}_hits.csv", index=False)
    
    # Compile summary
    summary = {
        'total_sgrnas': len(results),
        'significant_hits': len(hits),
        'depletions': len(depletions),
        'enrichments': len(enrichments),
        'qc_metrics': qc,
        'output_files': {
            'full_results': f"{output_prefix}_sgrna_results.csv",
            'hits': f"{output_prefix}_hits.csv"
        }
    }
    
    # Print summary
    print(f"\n{'='*60}")
    print("ANALYSIS SUMMARY")
    print(f"{'='*60}")
    print(f"Total sgRNAs: {summary['total_sgrnas']}")
    print(f"Significant hits (FDR<{fdr_threshold}, |LFC|>{lfc_threshold}): {summary['significant_hits']}")
    print(f"  - Depletions: {summary['depletions']}")
    print(f"  - Enrichments: {summary['enrichments']}")
    print(f"\nResults saved:")
    print(f"  - {summary['output_files']['full_results']}")
    print(f"  - {summary['output_files']['hits']}")
    print(f"{'='*60}")
    
    return summary

# Execute workflow
results = analyze_crispr_screen(
    counts_file="sgrna_counts.txt",
    samplesheet="samples.csv",
    control_samples=["Ctrl_1", "Ctrl_2", "Ctrl_3"],
    treatment_samples=["Drug_1", "Drug_2", "Drug_3"],
    output_prefix="drug_resistance_screen",
    fdr_threshold=0.05,
    lfc_threshold=1.0
)
```

**Expected Output Files:**

```
analysis_results/
├── drug_resistance_screen_sgrna_results.csv  # All sgRNA statistics
├── drug_resistance_screen_hits.csv          # Significant hits only
└── qc_report.txt                            # Quality control summary
```

---

## Common Patterns

### Pattern 1: Viability Screen (Essential Gene Identification)

**Scenario**: Identify genes essential for cell survival by comparing T0 (transduction) vs T14 (14 days post-transduction).

```json
{
  "screen_type": "viability",
  "comparison": "T14_vs_T0",
  "expected_depletions": "Essential genes (ribosomal, splicing, etc.)",
  "expected_enrichments": "None (unless suppressors of toxicity)",
  "positive_controls": ["RPL30", "RPS19", "PCNA"],
  "negative_controls": ["LacZ", "NTC"],
  "analysis_parameters": {
    "fdr_threshold": 0.05,
    "lfc_threshold": 1.0,
    "gene_aggregation": "mean"
  }
}
```

**Workflow:**
1. Collect cells at T0 (immediately after transduction)
2. Maintain parallel culture for 14 days (T14)
3. Harvest T14 cells when control cells reach confluence
4. Sequence both T0 and T14 samples
5. Analyze depletion of sgRNAs at T14 relative to T0
6. Identify genes with significantly depleted sgRNAs (essential genes)
7. Validate top hits with individual sgRNA validation

**Output Example:**
```
Essential Gene Screen Results:
  Total sgRNAs tested: 65,383
  Significantly depleted: 3,847 sgRNAs (FDR<0.05, LFC<-1)
  
Top Essential Genes:
  RPL30: mean LFC = -4.2, 5/5 sgRNAs significant
  RPS19: mean LFC = -3.8, 4/5 sgRNAs significant
  PCNA:  mean LFC = -3.5, 5/5 sgRNAs significant
  
QC Metrics:
  Gini index: 0.25 (excellent library representation)
  Read depth: 25M per sample (sufficient)
```

### Pattern 2: Drug Resistance Screen

**Scenario**: Identify genes whose knockout confers resistance to a cytotoxic drug (e.g., vemurafenib in BRAF-mutant melanoma).

```json
{
  "screen_type": "drug_resistance",
  "treatment": "vemurafenib (2 μM)",
  "control": "DMSO",
  "duration": "14 days",
  "expected_depletions": "Drug sensitizers, synthetic lethal",
  "expected_enrichments": "Drug resistance genes",
  "known_resistance_genes": ["NRAS", "MAP2K1", "MEK1"],
  "analysis_parameters": {
    "fdr_threshold": 0.05,
    "lfc_threshold": 1.0,
    "focus": "enrichments"
  }
}
```

**Workflow:**
1. Transduce cells with genome-wide sgRNA library
2. Split into drug-treated and DMSO control groups
3. Treat with drug at appropriate concentration (IC70-IC90)
4. Maintain for 2-3 weeks until control cells die
5. Harvest resistant colonies from drug-treated group
6. Compare sgRNA representation: Drug vs DMSO
7. Identify enriched sgRNAs (resistance genes)
8. Validate resistance with individual sgRNAs and drug dose-response

**Output Example:**
```
Drug Resistance Screen Results (Vemurafenib):
  Significant enrichments: 156 sgRNAs (FDR<0.05, LFC>1)
  
Top Resistance Genes:
  NRAS:   mean LFC = +2.8, 4/5 sgRNAs enriched
  MAP2K1: mean LFC = +2.5, 5/5 sgRNAs enriched
  MED12:  mean LFC = +2.1, 3/5 sgRNAs enriched
  
Validation recommended:
  - Test individual sgRNAs in dose-response assay
  - Confirm resistance phenotype with cell viability assay
  - Check for known resistance mechanisms
```

### Pattern 3: Drug Sensitivity/Synthetic Lethality Screen

**Scenario**: Identify genes that, when knocked out, sensitize cells to drug treatment (synthetic lethal interactions).

```json
{
  "screen_type": "drug_sensitivity",
  "treatment": "PARP inhibitor (olaparib)",
  "control": "DMSO",
  "cell_line": "BRCA1-mutant ovarian cancer",
  "expected_depletions": "DNA repair genes (synthetic lethal)",
  "expected_enrichments": "Drug resistance mechanisms",
  "known_synthetic_lethal": ["PARP1", "BRCA2", "PALB2"],
  "analysis_parameters": {
    "fdr_threshold": 0.05,
    "lfc_threshold": 1.0,
    "focus": "depletions"
  }
}
```

**Workflow:**
1. Transduce cells with sgRNA library
2. Treat with sub-lethal drug concentration (IC30)
3. Maintain for 2 weeks under drug selection
4. Compare sgRNA representation: Drug-treated vs control
5. Identify depleted sgRNAs (synthetic lethal/sensitizer genes)
6. Validate with individual sgRNAs and combination assays
7. Compare with genetic dependency maps (DepMap)

**Output Example:**
```
Synthetic Lethality Screen (Olaparib in BRCA1-mutant):
  Significant depletions: 234 sgRNAs (FDR<0.05, LFC<-1)
  
Top Synthetic Lethal Hits:
  BRCA2:   mean LFC = -3.2, 5/5 sgRNAs depleted
  PALB2:   mean LFC = -2.8, 4/5 sgRNAs depleted
  RAD51C:  mean LFC = -2.5, 5/5 sgRNAs depleted
  
Biological Interpretation:
  - Strong enrichment of homologous recombination genes
  - Consistent with known synthetic lethal interactions
  - Potential combination therapy targets identified
```

### Pattern 4: Comparative Screen (Cell Line vs Cell Line)

**Scenario**: Compare genetic dependencies between two cell lines to identify lineage-specific vulnerabilities.

```json
{
  "screen_type": "comparative",
  "comparison": "Melanoma_vs_Lung_cancer",
  "cell_lines": ["A375", "SKMEL28", "A549", "H1299"],
  "analysis_type": "differential_essentiality",
  "expected_lineage_specific": {
    "melanoma": ["MITF", "SOX10", "TYR"],
    "lung": ["NKX2-1", "TP63"]
  },
  "analysis_parameters": {
    "fdr_threshold": 0.05,
    "lfc_threshold": 1.0,
    "replicate_requirement": 2
  }
}
```

**Workflow:**
1. Perform viability screens in multiple cell lines in parallel
2. Normalize each screen independently
3. Compare gene-level essentiality scores across lines
4. Identify genes essential in one lineage but not another
5. Validate lineage-specific dependencies
6. Explore therapeutic relevance (tumor-type specific targets)

**Output Example:**
```
Comparative Screen: Melanoma vs Lung Cancer
  Melanoma-specific essential: 127 genes
  Lung-specific essential: 203 genes
  Common essential: 1,847 genes
  
Top Melanoma-Specific Dependencies:
  MITF:   LFC diff = -4.5 (essential in melanoma, not lung)
  SOX10:  LFC diff = -3.8
  TYR:    LFC diff = -3.2
  
Top Lung-Specific Dependencies:
  NKX2-1: LFC diff = -3.9
  TP63:   LFC diff = -3.1
  
Therapeutic Implications:
  - Lineage-specific targets identified
  - Potential for tumor-type selective therapy
```

---

## Quality Checklist

**Pre-Analysis Checks:**
- [ ] **CRITICAL**: Verify library composition matches expected sgRNA list
- [ ] Check sequencing depth (>10M reads per sample recommended)
- [ ] Confirm sample annotations match count matrix columns
- [ ] Verify control and treatment sample assignments are correct
- [ ] Check for batch effects (different sequencing runs, library preps)
- [ ] Review positive control performance (known essential genes)
- [ ] Confirm negative controls show no significant effects
- [ ] Validate replicate consistency (correlation >0.7 expected)

**During Analysis:**
- [ ] Calculate and review QC metrics (Gini, read depth, dropout)
- [ ] **CRITICAL**: Check Gini index <0.4 for library quality
- [ ] Examine LFC distribution for normality and outliers
- [ ] Verify positive controls are significantly depleted (viability screens)
- [ ] Check for batch effects using PCA or correlation heatmaps
- [ ] Apply appropriate statistical thresholds (FDR < 0.05 standard)
- [ ] Require multiple sgRNAs per gene for hit calling (≥2 recommended)
- [ ] Compare hit lists with published data for similar screens

**Post-Analysis Verification:**
- [ ] **CRITICAL**: Validate top hits show concordance across sgRNAs
- [ ] Check known positive controls are recovered
- [ ] Assess negative control performance (should not be significant)
- [ ] Compare replicate correlation for hits vs non-hits
- [ ] Review hit gene functions for biological plausibility
- [ ] Check for potential off-target effects (seed sequence analysis)
- [ ] Verify hit numbers are reasonable (10s-100s, not 1000s)
- [ ] Generate visualization (MA plots, volcano plots, heatmaps)

**Before Validation or Publication:**
- [ ] **CRITICAL**: Validate top 5-10 hits with individual sgRNAs
- [ ] Perform rescue experiments to confirm on-target effects
- [ ] Compare with orthogonal datasets (DepMap, published screens)
- [ ] Check for cell line-specific vs pan-essential classification
- [ ] Assess therapeutic relevance of identified hits
- [ ] Plan secondary screens if primary screen quality issues found
- [ ] Document all parameters and thresholds used
- [ ] Prepare data for public deposition (if applicable)

---

## Common Pitfalls

**Experimental Design Issues:**
- ❌ **Insufficient sequencing depth** → Poor statistical power, missed hits
  - ✅ Minimum 10M reads per sample; 20M+ for complex libraries
  
- ❌ **Library bottleneck** → Gini index >0.4, skewed representation
  - ✅ Maintain MOI <0.3; use sufficient cell numbers (500-1000x library coverage)
  
- ❌ **Inadequate replicates** → High variance, irreproducible results
  - ✅ Use ≥3 biological replicates per condition
  
- ❌ **Wrong time point** → Too early (no selection) or too late (extensive dropout)
  - ✅ Optimize time point based on doubling time and selection pressure

**Analysis Issues:**
- ❌ **Ignoring QC metrics** → Analyzing poor quality data
  - ✅ Always review Gini index, read depth, and dropout before analysis
  
- ❌ **Incorrect sample assignment** → Control/treatment mix-up
  - ✅ Double-check sample annotation file; validate with positive controls
  
- ❌ **Single sgRNA hits** → Potential off-target effects
  - ✅ Require ≥2 significant sgRNAs per gene; check concordance
  
- ❌ **Over-reliance on p-values** → Many false positives with large library
  - ✅ Use FDR correction; add LFC threshold; validate experimentally

**Interpretation Issues:**
- ❌ **Ignoring cell number effects** → Different growth rates confound results
  - ✅ Normalize for cell doublings; use appropriate controls
  
- ❌ **Off-target effects dominating** → False positive hits
  - ✅ Use improved libraries (e.g., Brunello, Brie); validate with rescue
  
- ❌ **Pan-essential vs selective** → Misclassifying broadly essential genes
  - ✅ Compare with DepMap data; use differential analysis for specificity
  
- ❌ **Not validating hits** → Publishing false positives
  - ✅ Validate top hits with individual sgRNAs; perform rescue experiments

**Technical Issues:**
- ❌ **Batch effects** → Confounding by library prep or sequencing batch
  - ✅ Randomize samples across batches; include batch in statistical model
  
- ❌ **Contamination** → Cross-sample contamination affects quantification
  - ✅ Use unique molecular identifiers (UMIs); check for index hopping
  
- ❌ **Reference genome mismatch** → sgRNAs not mapping correctly
  - ✅ Use same genome version as library design; check sgRNA sequences
  
- ❌ **Incomplete annotation** → sgRNAs missing gene mapping
  - ✅ Verify library annotation file is complete and current

---

## Troubleshooting

**Problem: No significant hits despite strong biological effect**
- Symptoms: Clear phenotype but no FDR-significant sgRNAs
- Causes:
  - High variance between replicates
  - Insufficient sequencing depth
  - Weak effect sizes
  - Stringent statistical thresholds
- Solutions:
  - Increase replicate number
  - Increase sequencing depth
  - Use more permissive FDR threshold (0.1)
  - Consider gene-level aggregation

**Problem: Too many significant hits (1000s)**
- Symptoms: Excessive number of hits, many likely false positives
- Causes:
  - Low variance (overdispersion underestimated)
  - Strong selection pressure
  - Library quality issues
  - Noisy data
- Solutions:
  - Use more stringent FDR threshold (0.01)
  - Increase LFC threshold (1.5 or 2.0)
  - Filter by sgRNA concordance
  - Review QC metrics and repeat if poor quality

**Problem: High Gini index (>0.4)**
- Symptoms: Library representation highly skewed
- Causes:
  - Library bottleneck at transduction
  - Insufficient cell numbers
  - High MOI leading to multiple integrations
- Solutions:
  - Use lower MOI (<0.3)
  - Increase cell numbers (500-1000x library size)
  - Improve transduction efficiency
  - Consider repeating screen

**Problem: Known essential genes not identified**
- Symptoms: Positive controls (RPL30, RPS19) not significantly depleted
- Causes:
  - Insufficient selection time
  - Library quality issues
  - Analysis errors
- Solutions:
  - Extend time point for viability screens
  - Check library composition and representation
  - Verify analysis parameters (control vs treatment assignment)

**Problem: Discordant sgRNAs for same gene**
- Symptoms: Only 1-2 of 5 sgRNAs significant for hit genes
- Causes:
  - Off-target effects
  - Variable sgRNA efficiency
  - Library design issues
- Solutions:
  - Require ≥3 significant sgRNAs for gene-level hits
  - Check sgRNA sequences for off-target potential
  - Use improved second-generation libraries
  - Validate with independent sgRNAs

**Problem: Batch effects between replicates**
- Symptoms: Low correlation between replicates of same condition
- Causes:
  - Different library prep batches
  - Different sequencing runs
  - Technical variation
- Solutions:
  - Include batch as covariate in analysis
  - Use ComBat or similar batch correction
  - Re-sequence inconsistent replicates
  - Randomize samples across batches in future

**Problem: Negative controls showing significant effects**
- Symptoms: Non-targeting controls (NTC) or safe-targeting sgRNAs in hit list
- Causes:
  - Technical artifacts
  - Random chance with large library
  - Library design issues
- Solutions:
  - Review NTC performance; should not be systematically enriched/depleted
  - If systematic, investigate technical issues
  - Use NTC distribution to set empirical thresholds

---

## References

Available in `references/` directory:

- (No reference files currently available for this skill)

**External Resources:**
- AddGene CRISPR Libraries: https://www.addgene.org/crispr/libraries/
- DepMap Portal: https://depmap.org/portal/
- MAGeCK Documentation: https://sourceforge.net/p/mageck/wiki/Home/
- BAGEL Algorithm: https://github.com/hart-lab/bagel
- CRISPR Screen Analysis Best Practices: https://pubmed.ncbi.nlm.nih.gov/29651053/

---

## Scripts

Located in `scripts/` directory:

- `main.py` - CRISPR screen analysis engine with QC, RRA, and hit identification

---

## Common CRISPR Screen Types

| Screen Type | Comparison | Expected Hits | Typical Duration |
|-------------|-----------|---------------|------------------|
| **Viability** | T14 vs T0 | Essential genes depleted | 10-14 days |
| **Drug Resistance** | Drug vs DMSO | Resistance genes enriched | 14-21 days |
| **Drug Sensitivity** | Drug vs DMSO | Sensitizers depleted | 14-21 days |
| **Comparative** | Cell A vs Cell B | Lineage-specific dependencies | 10-14 days |
| **Sensitizer** | Drug A+B vs Drug A | Combination targets | 10-14 days |

## Parameters

| Parameter | Type | Default | Required | Description |
|-----------|------|---------|----------|-------------|
| `--counts`, `-c` | string | - | Yes | sgRNA count matrix file |
| `--samples`, `-s` | string | - | Yes | Sample annotation file |
| `--control` | string | - | No | Control samples (comma-separated) |
| `--treatment`, `-t` | string | - | No | Treatment samples (comma-separated) |
| `--output`, `-o` | string | - | No | Output directory |
| `--fdr` | float | 0.05 | No | FDR threshold |

## Usage

### Basic Usage

```bash
# Analyze CRISPR screen data
python scripts/main.py --counts sgrna_counts.txt --samples samplesheet.csv

# With specific control and treatment
python scripts/main.py --counts counts.txt --samples samples.csv --control "Ctrl1,Ctrl2" --treatment "Treat1,Treat2"

# Custom FDR threshold
python scripts/main.py --counts counts.txt --samples samples.csv --fdr 0.01 --output ./results
```

## Risk Assessment

| Risk Indicator | Assessment | Level |
|----------------|------------|-------|
| Code Execution | Python script executed locally | Low |
| Network Access | No external API calls | Low |
| File System Access | Read count files, write results | Low |
| Data Exposure | Processes genomic screening data | Medium |
| PHI Risk | May contain cell line genetic info | Low |

## Security Checklist

- [x] No hardcoded credentials or API keys
- [x] No unauthorized file system access
- [x] Input validation for file paths
- [x] Output directory restricted
- [x] Error messages sanitized
- [x] Script execution in sandboxed environment

## Prerequisites

```bash
# Python 3.7+
numpy
pandas
scipy
```

## Evaluation Criteria

### Success Metrics
- [x] Successfully loads sgRNA count matrices
- [x] Calculates QC metrics (Gini index, zero counts)
- [x] Performs RRA analysis
- [x] Identifies significant hits with FDR control

### Test Cases
1. **Basic Analysis**: Count matrix + samplesheet → QC metrics + hit list
2. **RRA Analysis**: Control vs Treatment → Ranked gene list with p-values
3. **QC Metrics**: Count data → Gini scores, zero sgRNA counts

## Lifecycle Status

- **Current Stage**: Active
- **Next Review Date**: 2026-03-09
- **Known Issues**: None
- **Planned Improvements**:
  - Add MAGeCK integration
  - Support for multiple analysis methods
  - Enhanced visualization

---

**Last Updated**: 2026-02-09  
**Skill ID**: 183  
**Version**: 2.0 (K-Dense Standard)
