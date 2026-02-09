---
name: cnv-caller-plotter
description: Detect copy number variations from whole genome sequencing data and generate publication-quality genome-wide CNV plots. Supports CNV calling, segmentation, and visualization for cancer genomics and rare disease analysis.
allowed-tools: [Read, Write, Bash, Edit, Grep]
license: MIT
metadata:
  skill-author: AIPOCH
---

# CNV Caller & Plotter

Detect copy number variations (CNVs) from whole genome sequencing (WGS) data and generate genome-wide visualization plots for cancer genomics, rare disease analysis, and population genetics studies. Provides CNV calling, segmentation analysis, and publication-ready visualization.

**Key Capabilities:**
- **CNV Detection from WGS**: Identify copy number gains and losses from aligned sequencing data
- **Genomic Segmentation**: Divide genome into bins/windows for copy number estimation
- **Flexible Input Support**: Process BAM, VCF, and other standard genomics formats
- **Publication-Quality Plots**: Generate genome-wide CNV profiles in PNG, PDF, or SVG formats
- **Standard Output Formats**: Export CNV calls in BED format for downstream analysis

---

## When to Use

**✅ Use this skill when:**
- Analyzing **cancer genomes** to identify somatic copy number alterations (SCNAs)
- Studying **rare diseases** with suspected copy number variation etiology
- Performing **population genetics** studies comparing CNV frequencies across groups
- Generating **genome-wide CNV visualizations** for publications or reports
- Creating **BED format CNV calls** for integration with other analysis pipelines
- Performing **comparative CNV analysis** between tumor and normal samples
- Validating **CNV calls from SNP arrays** with sequencing data

**❌ Do NOT use when:**
- Working with **targeted sequencing panels** (exome/targeted capture) → Use specialized tools like CNVkit or ExomeDepth
- Detecting **structural variations** involving translocations or inversions → Use `structural-variant-caller`
- Analyzing **single-cell RNA-seq data** → Use single-cell specific CNV tools (e.g., inferCNV)
- Detecting **small indels** (<50bp) → Use `variant-caller` for small variant detection
- Requiring **clinical-grade CNV detection** for diagnostic purposes → Use validated clinical pipelines with proper QC
- Working with **low-coverage data** (<10x) → Results may be unreliable; consider SNP array-based methods

**Related Skills:**
- **上游 (Upstream)**: `fastqc-report-interpreter`, `alignment-quality-checker`, `variant-caller`
- **下游 (Downstream)**: `circos-plot-generator`, `go-kegg-enrichment`, `heatmap-beautifier`

---

## Integration with Other Skills

**Upstream Skills:**
- `fastqc-report-interpreter`: Assess sequencing quality before CNV calling; low quality data may produce unreliable CNVs
- `alignment-quality-checker`: Verify BAM file quality and coverage uniformity; uneven coverage causes CNV artifacts
- `variant-caller`: Generate SNV/indel calls for combined CNV-SNV analysis in cancer samples

**Downstream Skills:**
- `circos-plot-generator`: Create circular genome plots integrating CNVs with other genomic features
- `go-kegg-enrichment`: Perform pathway enrichment on genes within CNV regions
- `heatmap-beautifier`: Visualize CNV profiles across multiple samples

**Complete Workflow:**
```
Raw WGS Data → fastqc-report-interpreter → alignment-quality-checker → cnv-caller-plotter → circos-plot-generator → Publication Figures
```

---

## Core Capabilities

### 1. Copy Number Variation Detection

Identify genomic regions with copy number gains (amplifications) or losses (deletions) from WGS data by analyzing read depth patterns.

```python
from scripts.main import CNVCaller

# Initialize CNV caller with bin size
caller = CNVCaller(bin_size=1000)

# Call CNVs from BAM file
cnv_calls = caller.call_cnvs(
    input_file="sample.bam",
    reference="hg38.fa"
)

# Review detected CNVs
for cnv in cnv_calls:
    print(f"{cnv['chrom']}:{cnv['start']}-{cnv['end']}")
    print(f"  Copy Number: {cnv['cn']}")
    if cnv['cn'] > 2:
        print(f"  Type: Amplification (gain)")
    elif cnv['cn'] < 2:
        print(f"  Type: Deletion (loss)")
```

**Parameters:**

| Parameter | Type | Required | Description | Default |
|-----------|------|----------|-------------|---------|
| `input_file` | str | Yes | Path to input BAM or VCF file | None |
| `reference` | str | Yes | Path to reference genome FASTA | None |
| `bin_size` | int | No | Size of genomic bins for segmentation (bp) | 1000 |

**CNV Calling Strategy:**

| Approach | Best For | Sensitivity | Specificity |
|----------|----------|-------------|-------------|
| **Read Depth Analysis** | Large CNVs (>10kb) | High | Medium |
| **Paired-end Mapping** | Medium CNVs (1-10kb) | Medium | High |
| **Split-read Analysis** | Small CNVs (<1kb) | Medium | High |
| **Combined Approach** | Comprehensive detection | High | High |

**Best Practices:**
- ✅ **Use appropriate bin size**: 1000bp for WGS, smaller for targeted analysis
- ✅ **Ensure sufficient coverage**: Minimum 15-20x for reliable CNV detection
- ✅ **Match reference genome**: Use same reference as alignment (hg19 vs hg38)
- ✅ **Check coverage uniformity**: GC bias can cause false positive CNVs

**Common Issues and Solutions:**

**Issue: False positive CNVs in repetitive regions**
- Symptom: Many CNV calls in centromeres, telomeres, or segmental duplications
- Solution: Filter CNVs overlapping known problematic regions; use mappability filters

**Issue: Low sensitivity for small CNVs**
- Symptom: Missing CNVs <5kb despite adequate coverage
- Solution: Reduce bin size; use split-read or paired-end signals in addition to depth

### 2. Genomic Segmentation and Binning

Divide the genome into windows/bins for copy number estimation, enabling systematic analysis of the entire genome.

```python
from scripts.main import CNVCaller

# Different bin sizes for different applications
bin_configs = {
    "high_resolution": 100,    # For small CNV detection
    "standard": 1000,          # Default for WGS
    "low_resolution": 10000    # For large-scale alterations
}

for config_name, bin_size in bin_configs.items():
    caller = CNVCaller(bin_size=bin_size)
    print(f"\n{config_name} (bin_size={bin_size}bp):")
    
    # Calculate approximate number of bins for human genome
    genome_size = 3_000_000_000  # 3 Gb
    num_bins = genome_size // bin_size
    print(f"  Estimated bins: ~{num_bins:,}")
    print(f"  Resolution: {bin_size}bp")
```

**Bin Size Selection Guide:**

| Bin Size | Resolution | Use Case | Coverage Required |
|----------|------------|----------|-------------------|
| **100 bp** | High | Small CNVs (<5kb) | >30x |
| **1000 bp** | Standard | General WGS analysis | >15x |
| **10000 bp** | Low | Large chromosomal alterations | >5x |
| **Variable** | Adaptive | Mixed resolution | >20x |

**Best Practices:**
- ✅ **Match bin size to expected CNV size**: Use smaller bins for detecting small CNVs
- ✅ **Consider coverage depth**: Higher coverage enables smaller bins
- ✅ **Exclude unmappable regions**: Filter bins with zero or very low mappability
- ✅ **Normalize for GC content**: GC-rich regions have different coverage patterns

**Common Issues and Solutions:**

**Issue: Noisy segmentation due to small bins**
- Symptom: Erratic copy number estimates with high variance
- Solution: Increase bin size; apply smoothing algorithms; use larger bins for baseline

**Issue: Missing large CNVs with large bins**
- Symptom: Large deletions/amplifications not called when spanning multiple bins
- Solution: Use statistical segmentation (CBS, PSCBS) to join adjacent altered bins

### 3. Genome-Wide Visualization

Generate publication-quality plots showing copy number profiles across all chromosomes for visual interpretation and presentation.

```python
from scripts.main import CNVCaller

caller = CNVCaller(bin_size=1000)

# Example CNV calls for plotting
cnv_calls = [
    {"chrom": "chr1", "start": 1000000, "end": 2000000, "cn": 3},   # Gain
    {"chrom": "chr7", "start": 50000000, "end": 55000000, "cn": 1}, # Loss
    {"chrom": "chr17", "start": 35000000, "end": 36000000, "cn": 4} # High-level amplification
]

# Generate plots in different formats
output_dir = "./cnv_results"

for fmt in ["png", "pdf", "svg"]:
    plot_file = caller.plot_genome_wide(
        cnv_calls=cnv_calls,
        output_path=output_dir,
        fmt=fmt
    )
    print(f"Generated: {plot_file}")

# Plot features:
# - Genome-wide view with all chromosomes
# - Copy number on Y-axis (0-6 typical range)
# - Chromosomal position on X-axis
# - Color coding: red=loss, blue=gain, black=neutral
```

**Output Formats:**

| Format | Extension | Best For | File Size |
|--------|-----------|----------|-----------|
| **PNG** | .png | Web, presentations, quick viewing | Medium |
| **PDF** | .pdf | Publications, high-quality printing | Large |
| **SVG** | .svg | Vector editing, scalable graphics | Small |

**Best Practices:**
- ✅ **Use PDF for publications**: Vector format maintains quality at any zoom
- ✅ **Include baseline (CN=2)**: Reference line helps interpret gains/losses
- ✅ **Color-blind friendly palette**: Use distinct colors for gains vs losses
- ✅ **Annotate key regions**: Mark known cancer genes or regions of interest

**Common Issues and Solutions:**

**Issue: Plot too crowded with many CNVs**
- Symptom: Overlapping points make plot unreadable
- Solution: Use segmentation to merge adjacent calls; adjust point size/alpha

**Issue: ChrY not displayed for female samples**
- Symptom: Missing chromosome in plot for female subjects
- Solution: Dynamically detect sex from coverage; adjust plot accordingly

### 4. BED Format Export

Export CNV calls in standard BED format for compatibility with genome browsers and downstream analysis tools.

```python
from scripts.main import CNVCaller

caller = CNVCaller()

# Example CNV calls
cnv_calls = [
    {"chrom": "chr1", "start": 1000000, "end": 2000000, "cn": 3},
    {"chrom": "chr7", "start": 50000000, "end": 55000000, "cn": 1},
]

# Export to BED format
bed_file = caller.save_bed(cnv_calls, "./output")

# BED format structure:
# chrom  start    end       name      score  strand
# chr1   1000000  2000000   CN=3      .      .
# chr7   50000000 55000000  CN=1      .      .

print(f"BED file saved: {bed_file}")

# Read and display BED content
with open(bed_file, 'r') as f:
    print("\nBED file content:")
    for line in f:
        print(line.strip())
```

**BED Format Specification:**

| Column | Field | Description | Example |
|--------|-------|-------------|---------|
| 1 | chrom | Chromosome name | chr1, chrX |
| 2 | start | Start position (0-based) | 1000000 |
| 3 | end | End position (1-based) | 2000000 |
| 4 | name | CNV annotation | CN=3 |
| 5 | score | Optional quality score | . |
| 6 | strand | Strand info (usually .) | . |

**Best Practices:**
- ✅ **Use 0-based coordinates**: Standard BED format uses 0-based start, 1-based end
- ✅ **Include copy number in name**: Makes CNV status immediately visible
- ✅ **Sort by chromosome and position**: Required for many tools (bedtools, IGV)
- ✅ **Validate format**: Check with `bedtools` or genome browser before distribution

**Common Issues and Solutions:**

**Issue: BED file rejected by genome browser**
- Symptom: IGV or UCSC Genome Browser shows error loading BED
- Solution: Ensure proper chromosome naming (chr1 vs 1); sort file; check for tabs vs spaces

**Issue: Coordinate system confusion**
- Symptom: CNVs appear shifted by 1bp in different tools
- Solution: BED is 0-based, GFF/VCF are 1-based; convert if necessary

### 5. Tumor-Normal Comparison

Compare CNV profiles between tumor and matched normal samples to identify somatic copy number alterations (SCNAs).

```python
from scripts.main import CNVCaller

caller = CNVCaller(bin_size=1000)

# Call CNVs in tumor and normal samples
tumor_cnvs = caller.call_cnvs("tumor.bam", "hg38.fa")
normal_cnvs = caller.call_cnvs("normal.bam", "hg38.fa")

# Identify somatic CNVs (present in tumor, not in normal)
def find_somatic_cnvs(tumor_calls, normal_calls):
    """Identify CNVs present in tumor but not normal."""
    somatic_cnvs = []
    
    for t_cnv in tumor_calls:
        is_somatic = True
        
        # Check if similar CNV exists in normal
        for n_cnv in normal_calls:
            if (t_cnv['chrom'] == n_cnv['chrom'] and
                abs(t_cnv['start'] - n_cnv['start']) < 10000 and
                abs(t_cnv['end'] - n_cnv['end']) < 10000 and
                t_cnv['cn'] == n_cnv['cn']):
                is_somatic = False
                break
        
        if is_somatic:
            somatic_cnvs.append(t_cnv)
    
    return somatic_cnvs

somatic_cnvs = find_somatic_cnvs(tumor_cnvs, normal_cnvs)

print(f"Total tumor CNVs: {len(tumor_cnvs)}")
print(f"Somatic CNVs: {len(somatic_cnvs)}")

# Categorize somatic alterations
amplifications = [c for c in somatic_cnvs if c['cn'] > 2]
deletions = [c for c in somatic_cnvs if c['cn'] < 2]

print(f"  Amplifications: {len(amplifications)}")
print(f"  Deletions: {len(deletions)}")
```

**Somatic vs Germline Classification:**

| Category | Tumor CN | Normal CN | Interpretation |
|----------|----------|-----------|----------------|
| **Somatic Amplification** | >2 | 2 | Tumor-specific gain |
| **Somatic Deletion** | <2 | 2 | Tumor-specific loss |
| **Germline CNV** | ≠2 | ≠2 | Inherited CNV |
| **LOH** | 1 | 2 | Loss of heterozygosity |

**Best Practices:**
- ✅ **Use matched normal when available**: Essential for distinguishing somatic vs germline
- ✅ **Consider tumor purity**: Low purity samples have attenuated CNV signals
- ✅ **Validate key findings**: Use orthogonal methods (FISH, qPCR) for important CNVs
- ✅ **Account for clonality**: Subclonal CNVs may be present at lower frequencies

**Common Issues and Solutions:**

**Issue: Normal sample contamination in tumor**
- Symptom: CNV signals weaker than expected; fractional copy numbers
- Solution: Estimate tumor purity; use purity-corrected CNV calling

**Issue: Germline CNVs misclassified as somatic**
- Symptom: Many "somatic" CNVs that look like common polymorphisms
- Solution: Filter against population CNV databases (DGV, gnomAD-SV)

### 6. Quality Control and Filtering

Apply quality filters to remove artifactual CNV calls and improve result reliability.

```python
from scripts.main import CNVCaller

caller = CNVCaller()

# Example raw CNV calls with QC metrics
cnv_calls = [
    {
        "chrom": "chr1", "start": 1000000, "end": 2000000, "cn": 3,
        "quality_score": 50, "supporting_reads": 150
    },
    {
        "chrom": "chr7", "start": 50000000, "end": 50001000, "cn": 0,
        "quality_score": 10, "supporting_reads": 5  # Likely artifact
    },
]

# Apply quality filters
def filter_cnvs(cnv_list, min_quality=20, min_size=1000, min_support=20):
    """Filter CNVs based on quality metrics."""
    filtered = []
    
    for cnv in cnv_list:
        size = cnv['end'] - cnv['start']
        quality = cnv.get('quality_score', 0)
        support = cnv.get('supporting_reads', 0)
        
        # Apply filters
        if quality < min_quality:
            continue
        if size < min_size:
            continue
        if support < min_support:
            continue
            
        filtered.append(cnv)
    
    return filtered

# Filter with different stringencies
for min_q in [10, 20, 30]:
    filtered = filter_cnvs(cnv_calls, min_quality=min_q)
    print(f"Quality >= {min_q}: {len(filtered)} CNVs retained")

# Additional filters to consider:
# - Exclude segmental duplications
# - Exclude centromeres and telomeres
# - Minimum number of supporting bins
# - Concordance with paired-end or split-read signals
```

**Quality Metrics:**

| Metric | Threshold | Purpose |
|--------|-----------|---------|
| **Quality Score** | >20 | Overall confidence in CNV call |
| **Size** | >1kb | Remove small artifactual calls |
| **Supporting Reads** | >20 | Sufficient evidence depth |
| **Log2 Ratio** | |0.3| | Significant deviation from diploid |
| **Mappability** | >0.8 | Reliable unique mapping |

**Best Practices:**
- ✅ **Apply size filters**: Remove CNVs <1kb (often artifacts)
- ✅ **Filter repetitive regions**: Exclude known problematic regions
- ✅ **Use multiple evidence types**: Combine depth, paired-end, and split-read signals
- ✅ **Validate high-impact CNVs**: Use orthogonal methods for therapeutic targets

**Common Issues and Solutions:**

**Issue: Too many low-quality CNV calls**
- Symptom: Hundreds or thousands of CNVs called
- Solution: Increase quality thresholds; apply population frequency filters

**Issue: True CNVs filtered out**
- Symptom: Known cancer driver CNVs missing from results
- Solution: Use gene-specific filters; manually review regions of interest

---

## Complete Workflow Example

**From WGS data to CNV visualization:**

```bash
# Step 1: Call CNVs from tumor sample
python scripts/main.py \
  --input tumor_sample.bam \
  --reference hg38.fa \
  --output tumor_cnv/ \
  --bin-size 1000 \
  --plot-format pdf

# Step 2: Call CNVs from matched normal
python scripts/main.py \
  --input normal_sample.bam \
  --reference hg38.fa \
  --output normal_cnv/ \
  --bin-size 1000

# Step 3: Compare and identify somatic CNVs
# (Use Python API for comparison logic)

# Step 4: Generate final plots
python scripts/main.py \
  --input tumor_sample.bam \
  --reference hg38.fa \
  --output final_results/ \
  --plot-format pdf
```

**Python API Usage:**

```python
from scripts.main import CNVCaller
from pathlib import Path

def analyze_cancer_genome(
    tumor_bam: str,
    normal_bam: str,
    reference: str,
    output_dir: str
) -> dict:
    """
    Complete cancer genome CNV analysis workflow.
    """
    caller = CNVCaller(bin_size=1000)
    
    # Create output directory
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    # Call CNVs in both samples
    print("Calling CNVs in tumor sample...")
    tumor_cnvs = caller.call_cnvs(tumor_bam, reference)
    
    print("Calling CNVs in normal sample...")
    normal_cnvs = caller.call_cnvs(normal_bam, reference)
    
    # Identify somatic alterations
    somatic_cnvs = identify_somatic(tumor_cnvs, normal_cnvs)
    
    # Generate outputs
    tumor_bed = caller.save_bed(tumor_cnvs, output_dir)
    somatic_bed = caller.save_bed(somatic_cnvs, f"{output_dir}/somatic")
    plot_file = caller.plot_genome_wide(tumor_cnvs, output_dir, "pdf")
    
    # Calculate statistics
    stats = {
        "total_tumor_cnvs": len(tumor_cnvs),
        "somatic_cnvs": len(somatic_cnvs),
        "amplifications": len([c for c in somatic_cnvs if c['cn'] > 2]),
        "deletions": len([c for c in somatic_cnvs if c['cn'] < 2]),
        "output_files": {
            "tumor_bed": tumor_bed,
            "somatic_bed": somatic_bed,
            "genome_plot": plot_file
        }
    }
    
    return stats

# Execute workflow
results = analyze_cancer_genome(
    tumor_bam="tumor.bam",
    normal_bam="normal.bam",
    reference="hg38.fa",
    output_dir="./cnv_analysis"
)

print(f"\nAnalysis complete!")
print(f"Total tumor CNVs: {results['total_tumor_cnvs']}")
print(f"Somatic CNVs: {results['somatic_cnvs']}")
print(f"  Amplifications: {results['amplifications']}")
print(f"  Deletions: {results['deletions']}")
```

**Expected Output Files:**

```
cnv_analysis/
├── cnv_calls.bed              # All CNV calls in BED format
├── somatic/
│   └── cnv_calls.bed         # Somatic CNVs only
├── cnv_plot.pdf              # Genome-wide visualization
└── analysis_summary.json     # Statistics and metadata
```

---

## Common Patterns

### Pattern 1: Cancer Genome Analysis (Tumor-Normal Pair)

**Scenario**: Identify somatic copy number alterations in a cancer sample compared to matched normal tissue.

```json
{
  "analysis_type": "cancer_genome",
  "samples": {
    "tumor": "tumor_wgs.bam",
    "normal": "blood_normal.bam"
  },
  "reference": "hg38.fa",
  "parameters": {
    "bin_size": 1000,
    "min_cnv_size": 10000,
    "plot_format": "pdf"
  },
  "expected_outputs": [
    "Somatic CNV calls (BED format)",
    "Genome-wide CNV profile plot",
    "CNV statistics and summary"
  ]
}
```

**Workflow:**
1. Process both tumor and normal BAM files
2. Call CNVs in each sample independently
3. Compare to identify somatic alterations
4. Filter germline polymorphisms against population databases
5. Annotate cancer genes within CNV regions
6. Generate publication-quality visualization
7. Validate key driver alterations with orthogonal methods

**Output Example:**
```
Somatic CNV Summary:
  Total alterations: 47
  Amplifications: 12 (including MYC, EGFR)
  Deletions: 35 (including TP53, PTEN)
  
High-impact alterations:
  chr8:128000000-129000000 CN=8 (MYC amplification)
  chr17:7000000-8000000 CN=0 (TP53 deletion)
```

### Pattern 2: Rare Disease CNV Detection

**Scenario**: Detect pathogenic CNVs in a patient with suspected genomic disorder.

```json
{
  "analysis_type": "rare_disease",
  "sample": "patient.bam",
  "reference": "hg38.fa",
  "parameters": {
    "bin_size": 500,
    "min_cnv_size": 1000,
    "max_frequency": 0.01
  },
  "annotation": [
    "OMIM genes",
    "ClinVar pathogenic variants",
    "Decipher syndromes"
  ]
}
```

**Workflow:**
1. Call CNVs with high sensitivity settings
2. Filter against common population CNVs (DGV, gnomAD)
3. Prioritize rare CNVs (<1% frequency)
4. Annotate with disease-associated genes
5. Assess inheritance pattern (if parental data available)
6. Cross-reference with phenotype/HPO terms
7. Generate clinical report with prioritized findings

**Output Example:**
```
Rare CNV Findings:
  chr22:19000000-21000000 CN=1 (22q11.2 deletion syndrome)
    Size: 2.0 Mb
    Genes: TBX1, COMT, etc.
    Frequency: <0.1% in population
    Phenotype match: Cardiac, thymic, facial anomalies
    Classification: Pathogenic
```

### Pattern 3: Population CNV Analysis

**Scenario**: Compare CNV profiles across multiple samples to identify recurrent alterations.

```json
{
  "analysis_type": "population",
  "samples": [
    "sample1.bam", "sample2.bam", "sample3.bam",
    ...
  ],
  "cohorts": {
    "cases": 50,
    "controls": 50
  },
  "parameters": {
    "bin_size": 1000,
    "plot_format": "png"
  },
  "analysis": [
    "Recurrent CNV detection",
    "Burden analysis",
    "Association testing"
  ]
}
```

**Workflow:**
1. Call CNVs in all samples with consistent parameters
2. Merge and harmonize CNV calls across samples
3. Identify recurrent CNV regions
4. Perform burden analysis (total CNV load)
5. Test association with phenotype/status
6. Correct for multiple testing
7. Visualize CNV landscape across cohort

**Output Example:**
```
Population CNV Analysis:
  Samples analyzed: 100
  Total CNVs detected: 2,847
  
Recurrent alterations:
  chr1:1000000-2000000: 23% frequency
  chr16:15000000-16000000: 18% frequency
  
Case vs Control association:
  Significant enrichment: 3 CNV regions
  Most significant: chr8:128000000-129000000 (p=0.001)
```

### Pattern 4: Cell Line Characterization

**Scenario**: Characterize CNV profile of a cancer cell line for research or quality control.

```json
{
  "analysis_type": "cell_line",
  "sample": "mcf7_cell_line.bam",
  "reference": "hg38.fa",
  "parameters": {
    "bin_size": 1000,
    "plot_format": "pdf"
  },
  "comparison": {
    "reference_profile": "mcf7_ccle_cnvs.bed",
    "expected_alterations": ["chr8_MYC_amp", "chr20_ZNF217_amp"]
  }
}
```

**Workflow:**
1. Generate high-quality CNV profile from WGS
2. Compare to reference profiles (CCLE, COSMIC)
3. Verify expected cancer driver alterations
4. Identify subclonal populations
5. Assess genome stability metrics
6. Generate QC report for cell line authentication
7. Document for reproducibility

**Output Example:**
```
Cell Line: MCF-7
Identity confirmed: Yes (99.2% match to reference)

Expected alterations detected:
  chr8:128000000-129000000: CN=8 (MYC) ✓
  chr20:50000000-52000000: CN=6 (ZNF217) ✓

Additional alterations:
  chr17:35000000-37000000: CN=3 (ERBB2) ✓
  
Ploidy: 2.8 (aneuploid)
Genome instability score: High
```

---

## Quality Checklist

**Pre-analysis Checks:**
- [ ] **CRITICAL**: Verify input BAM file is properly aligned and indexed
- [ ] Confirm reference genome version matches alignment (hg19 vs hg38)
- [ ] Check sequencing coverage is sufficient (>15x for WGS, >30x for high resolution)
- [ ] Assess coverage uniformity (low uniformity causes CNV artifacts)
- [ ] Review FASTQC reports for quality issues
- [ ] Ensure matched normal sample is available for cancer analysis
- [ ] Verify sample identity (check sex chromosomes match metadata)
- [ ] Confirm no sample swaps or contamination

**During Analysis:**
- [ ] Select appropriate bin size for expected CNV size and coverage
- [ ] Apply GC content normalization if necessary
- [ ] Check for batch effects if analyzing multiple samples
- [ ] Monitor for high false positive rates in repetitive regions
- [ ] Validate sex chromosome calls against known sex
- [ ] Assess mitochondrial CNVs as quality control metric
- [ ] Review coverage plots for technical artifacts
- [ ] Check concordance with SNP array data if available

**Post-analysis Verification:**
- [ ] **CRITICAL**: Filter CNVs in known problematic regions (centromeres, telomeres)
- [ ] Remove common germline CNVs using population databases (DGV, gnomAD)
- [ ] Validate cancer driver alterations in known genes
- [ ] Check for CNV calls that disrupt single exons (often artifacts)
- [ ] Review very large CNVs (>50Mb) for technical artifacts
- [ ] Assess CNV burden against population norms
- [ ] Verify BED file format compliance
- [ ] Generate and review genome-wide plots

**Before Clinical or Publication Use:**
- [ ] **CRITICAL**: Have results reviewed by experienced analyst
- [ ] Validate pathogenic CNVs with orthogonal methods (FISH, qPCR, MLPA)
- [ ] Cross-reference with clinical databases (ClinVar, OMIM, Decipher)
- [ ] Document all parameters and filters applied
- [ ] Assess reproducibility by re-running with different parameters
- [ ] Check for batch effects in multi-sample analyses
- [ ] Confirm CNV coordinates with latest genome build
- [ ] Archive raw data and analysis scripts for reproducibility

---

## Common Pitfalls

**Input Data Issues:**
- ❌ **Using low coverage data** → Noisy CNV calls with many false positives
  - ✅ Minimum 15-20x coverage for reliable WGS CNV calling
  
- ❌ **Mismatched reference genomes** → CNVs called in wrong coordinates
  - ✅ Verify BAM uses same reference as CNV caller (hg19 vs hg38)
  
- ❌ **Not using matched normal for tumors** → Cannot distinguish somatic vs germline
  - ✅ Always use matched normal when available; use population controls otherwise
  
- ❌ **Poor coverage uniformity** → GC bias causes false CNVs
  - ✅ Check coverage plots; apply GC correction algorithms

**Analysis Parameter Issues:**
- ❌ **Bin size too large** → Miss small CNVs (<10kb)
  - ✅ Use 100-500bp bins for high-resolution analysis; 1000bp for standard WGS
  
- ❌ **Bin size too small** → Excessive noise in low coverage regions
  - ✅ Balance resolution with coverage; use adaptive binning if available
  
- ❌ **Inadequate quality filtering** → Too many false positive CNVs
  - ✅ Apply minimum quality scores; filter by size and read support
  
- ❌ **Not filtering common CNVs** → Report common polymorphisms as pathogenic
  - ✅ Filter against DGV, gnomAD, and other population databases

**Interpretation Issues:**
- ❌ **Ignoring tumor purity** → Misinterpret subclonal CNVs
  - ✅ Estimate tumor purity; adjust CNV calling thresholds accordingly
  
- ❌ **Not validating key findings** → Report false positive driver alterations
  - ✅ Validate cancer-relevant CNVs with orthogonal methods
  
- ❌ **Over-interpreting small CNVs** → Single-exon deletions are often artifacts
  - ✅ Focus on larger CNVs (>10kb) unless supported by multiple evidence types
  
- ❌ **Ignoring parental data** → Cannot determine inheritance in rare disease
  - ✅ Include parental samples for de novo vs inherited classification

**Output and Reporting Issues:**
- ❌ **Unclear coordinate system** → Confusion between 0-based and 1-based
  - ✅ Clearly document coordinate system used; BED is 0-based, VCF is 1-based
  
- ❌ **Missing quality metrics** → Cannot assess confidence in CNV calls
  - ✅ Include quality scores, supporting reads, and log2 ratios
  
- ❌ **Not archiving raw data** → Results cannot be reproduced
  - ✅ Save BAM files, parameter settings, and analysis scripts
  
- ❌ **Inadequate documentation** → Others cannot interpret results
  - ✅ Document all filters, thresholds, and databases used

---

## Troubleshooting

**Problem: No CNVs detected**
- Symptoms: Empty or nearly empty CNV call set
- Causes:
  - Coverage too low (<10x)
  - Bin size too large for small CNVs
  - Quality thresholds too stringent
  - Sample is actually diploid with no CNVs
- Solutions:
  - Verify coverage depth from BAM file
  - Reduce bin size for higher resolution
  - Relax quality filters temporarily
  - Check coverage uniformity across genome

**Problem: Too many CNV calls (hundreds or thousands)**
- Symptoms: Excessive number of CNV calls, many small or low-quality
- Causes:
  - Low coverage or high noise
  - Bin size too small
  - No quality filtering applied
  - Sample from highly polymorphic population
- Solutions:
  - Apply minimum quality score filter (Q>20)
  - Filter by minimum size (>1kb)
  - Remove calls in segmental duplications
  - Filter against population CNV databases

**Problem: False positives in repetitive regions**
- Symptoms: CNVs concentrated in centromeres, telomeres, or SDs
- Causes:
  - Low mappability in repetitive regions
  - Uneven coverage due to alignment issues
  - Reference genome gaps
- Solutions:
  - Filter CNVs overlapping known problematic regions
  - Use mappability filters (require mappability >0.8)
  - Exclude centromeres and telomeres from analysis
  - Use high-mappability reads only

**Problem: CNV signals too weak in tumor samples**
- Symptoms: Known cancer alterations not detected or weak signal
- Causes:
  - Low tumor purity (<20%)
  - Normal cell contamination
  - Subclonal alterations at low frequency
- Solutions:
  - Estimate tumor purity from VAF distribution
  - Use purity-corrected CNV calling
  - Lower thresholds for detection
  - Consider single-cell sequencing for subclonal analysis

**Problem: Sex chromosomes have unexpected copy numbers**
- Symptoms: XX sample showing CN=1 for X, or XY showing CN=2
- Causes:
  - Sex chromosome aneuploidy (e.g., Klinefelter, Turner syndromes)
  - Mislabeled sample sex
  - Pseudoautosomal region miscalls
- Solutions:
  - Verify sample sex from coverage ratios (X/Y)
  - Check clinical records for known sex chromosome abnormalities
  - Exclude pseudoautosomal regions from analysis
  - Analyze autosomes and sex chromosomes separately

**Problem: Batch effects in multi-sample analysis**
- Symptoms: CNV patterns correlate with sequencing batch rather than biology
- Causes:
  - Different sequencing platforms or chemistries
  - Coverage differences between batches
  - Different alignment parameters
- Solutions:
  - Normalize coverage across batches
  - Use same alignment and processing pipeline for all samples
  - Include batch as covariate in association testing
  - Perform batch correction algorithms

**Problem: Cannot install or run tool**
- Symptoms: Import errors, missing dependencies, execution failures
- Causes:
  - Missing Python packages (pysam, numpy, matplotlib)
  - Incompatible Python version
  - Missing reference genome index files
- Solutions:
  - Install required packages: `pip install pysam numpy matplotlib pandas`
  - Use Python 3.8 or higher
  - Create reference genome index: `samtools faidx reference.fa`
  - Check BAM file index exists: `sample.bam.bai`

---

## References

Available in `references/` directory:

- (No reference files currently available for this skill)

**External Resources:**
- Database of Genomic Variants (DGV): http://dgv.tcag.ca
- gnomAD Structural Variants: https://gnomad.broadinstitute.org
- ClinVar: https://www.ncbi.nlm.nih.gov/clinvar
- DECIPHER: https://www.deciphergenomics.org
- COSMIC: https://cancer.sanger.ac.uk

---

## Scripts

Located in `scripts/` directory:

- `main.py` - Main CNV calling and plotting engine

---

## CNV Detection Methods Comparison

| Method | Input | Sensitivity | Resolution | Best For |
|--------|-------|-------------|------------|----------|
| **Read Depth (this tool)** | BAM | Medium | 1-10 kb | Large CNVs, WGS |
| **Paired-end Mapping** | BAM | Medium | 100bp-10kb | Deletions, insertions |
| **Split-read Analysis** | BAM | High | 1bp-1kb | Breakpoint detection |
| **SNP Array** | CEL/IDAT | High | 5-25kb | Cost-effective screening |
| **Optical Mapping** | Bionano | High | 500bp+ | Very large SVs |

---

**Last Updated**: 2026-02-09  
**Skill ID**: 162  
**Version**: 2.0 (K-Dense Standard)
