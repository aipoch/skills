---
name: adverse-event-narrative
description: Generate CIOMS-compliant adverse event narratives for Individual Case 
  Safety Reports (ICSR). Creates structured pharmacovigilance documents following 
  CIOMS I and ICH E2B standards from case data for regulatory submission to 
  health authorities.
allowed-tools: [Read, Write, Bash, Edit]
license: MIT
metadata:
    skill-author: AIPOCH
---

# Adverse Event Narrative Generator

## Overview

Regulatory-grade narrative generation tool that transforms adverse event case data into CIOMS-compliant ICSR narratives suitable for submission to FDA, EMA, and other health authorities.

**Key Capabilities:**
- **CIOMS I Compliance**: Standardized narrative structure per international guidelines
- **ICH E2B Integration**: Electronic submission format compatibility
- **Temporal Analysis**: Timeline reconstruction and causality assessment
- **Medical Accuracy**: Clinical terminology and MedDRA coding
- **Multi-Case Processing**: Batch narrative generation for periodic reporting
- **Quality Validation**: Automated checks for completeness and consistency

## When to Use

**✅ Use this skill when:**
- Drafting ICSR narratives for regulatory submissions
- Converting safety case data to standardized text
- Preparing adverse event reports for health authorities
- Generating case summaries for signal detection
- Creating pharmacovigilance documentation for clinical trials
- Standardizing narrative format across safety teams
- Training new drug safety associates on narrative writing

**❌ Do NOT use when:**
- Case requires medical judgment or causality assessment → Use qualified safety physician
- Narrative for litigation or legal proceedings → Use legal documentation standards
- Patient-facing communications → Use `lay-summary-gen`
- Aggregate safety summaries → Use `safety-summary-reports`
- Coding MedDRA terms from verbatim → Use `meddra-coder`

**Integration:**
- **Upstream**: `meddra-coder` (MedDRA term coding), `clinical-data-cleaner` (data preparation)
- **Downstream**: `safety-summary-reports` (aggregate analysis), `regulatory-submission-prep` (FDA/EMA filing)

## Core Capabilities

### 1. CIOMS I Narrative Structure

Generate standardized sections per CIOMS guidelines:

```python
from scripts.narrative_generator import NarrativeGenerator

generator = NarrativeGenerator()

# Generate complete narrative
narrative = generator.generate(
    case_data=case_json,
    format="cioms_i",  # or "ich_e2b", "fda_medwatch"
    include_meddra=True
)

narrative.save("ICSR_2024_001_narrative.txt")
```

**Standard Sections:**
1. **Patient Demographics** - Age, sex, weight, relevant characteristics
2. **Medical History** - Significant pre-existing conditions
3. **Concomitant Medications** - Other drugs at time of event
4. **Suspect Drug(s)** - Medication(s) in question with dosing
5. **Adverse Event** - Detailed reaction description with MedDRA terms
6. **Diagnostic Results** - Lab values, imaging, procedures
7. **Treatment** - Medical management of the event
8. **Dechallenge/Rechallenge** - Effect of drug withdrawal/reintroduction
9. **Outcome** - Final patient status and sequelae
10. **Causality Assessment** - Reporter's relationship evaluation

### 2. Temporal Relationship Analysis

Reconstruct timeline and assess temporal plausibility:

```python
# Analyze temporal relationships
timeline = generator.analyze_timeline(
    drug_start="2024-01-15",
    drug_stop="2024-02-01",
    ae_onset="2024-01-28",
    dechallenge_date="2024-02-01",
    rechallenge_date=None
)

# Output shows temporal assessment
# "AE onset 13 days after drug initiation, positive dechallenge within 24h"
```

**Assessments Generated:**
- Time to onset (latency period)
- Dechallenge response (positive/negative/unknown)
- Rechallenge response (if applicable)
- Temporal plausibility (consistent with known drug profile)

### 3. Causality Evaluation Support

Structure causality assessment per WHO-UMC criteria:

```python
# Generate causality section
causality = generator.assess_causality(
    case_data=case,
    criteria="who_umc",  # or "naranjo", "cochrane"
    include_rationale=True
)

# Output structured assessment with points for each criterion
```

**WHO-UMC Categories:**
- **Certain** - Event reproduced on rechallenge
- **Probable/Likely** - Reasonable time, positive dechallenge, alternative causes unlikely
- **Possible** - Compatible time, but alternative causes possible
- **Unlikely** - Incompatible time or alternative cause probable
- **Conditional/Unclassified** - Insufficient information
- **Unassessable/Unclassifiable** - Data contradictory or incomplete

### 4. Multi-Format Output

Generate narratives for different regulatory contexts:

```python
# FDA MedWatch Form 3500A
fda_narrative = generator.generate(
    case_data=case,
    format="fda_medwatch",
    max_length=2000  # Character limit
)

# EMA E2B(R3) electronic format
ema_narrative = generator.generate(
    case_data=case,
    format="ich_e2b",
    version="R3"
)

# CIOMS I paper format
cioms_narrative = generator.generate(
    case_data=case,
    format="cioms_i"
)
```

## Common Patterns

### Pattern 1: Serious Adverse Event (Hospitalization)

**Scenario**: Patient hospitalized for severe drug reaction.

```json
{
  "case_id": "2024-SAE-001",
  "patient_age": "58 years",
  "patient_sex": "Female",
  "suspect_drugs": [{
    "drug_name": "Metformin",
    "dose": "1000 mg BID",
    "dates": "2024-01-15 to 2024-02-01"
  }],
  "adverse_events": [{
    "meddra_pt": "Lactic acidosis",
    "seriousness": "Hospitalization",
    "onset": "2024-01-28"
  }],
  "outcome": "Recovered with sequelae"
}
```

**Narrative Emphasis:**
- Hospitalization details (admission/discharge dates)
- Severity markers (ICU stay, intubation)
- Lactate levels and trend
- Renal function status
- Complete recovery timeline

### Pattern 2: Fatal Outcome Case

**Scenario**: Death suspected to be drug-related.

```python
# Fatal case handling
narrative = generator.generate(
    case_data=fatal_case,
    format="cioms_i",
    include_autopsy=True,
    cause_of_death_analysis=True
)
```

**Critical Elements:**
- Complete medical history relevant to death
- Concomitant medications contributing
- Autopsy findings (if performed)
- Cause of death per death certificate
- Alternative causes ruled out
- Reporter's opinion on contribution to death

### Pattern 3: Rechallenge Case

**Scenario**: Positive rechallenge confirms drug causation.

**Key Documentation:**
- First exposure dates and reaction
- Dechallenge response
- Rechallenge dates and circumstances
- Recurrence of same reaction
- Any differences in severity
- Conclusion on causality

**Narrative Structure:**
```
First Exposure:
- Drug X initiated [date]
- AE occurred [date], [description]
- Drug discontinued [date]
- Dechallenge: [positive/negative]

Rechallenge:
- Drug X reintroduced [date]
- Same AE recurred [date]
- Drug discontinued [date]
- Outcome: [status]

Causality: Certain (positive rechallenge)
```

### Pattern 4: Multi-Drug Reaction

**Scenario**: Multiple suspect drugs, need to identify most likely culprit.

**Analysis Approach:**
- Temporal sequence of drug initiations
- Time to onset for each drug
- Known adverse reaction profiles
- Dechallenge attempts (if any)
- Rechallenge data (if any)
- Reporter's primary suspect

**Narrative Organization:**
1. List all suspect drugs with indication and dates
2. Describe temporal relationship for each
3. Discuss dechallenge/rechallenge for each
4. Present reporter's causality assessment
5. Include alternative explanations

## Complete Workflow Example

**From case data to regulatory submission:**

```bash
# Step 1: Generate narrative
python scripts/main.py \
  --input case_data.json \
  --format cioms_i \
  --output narrative.txt

# Step 2: Validate completeness
python scripts/validate.py \
  --narrative narrative.txt \
  --check cioms_completeness \
  --output validation_report.txt

# Step 3: Generate E2B format for electronic submission
python scripts/main.py \
  --input case_data.json \
  --format ich_e2b \
  --output e2b_narrative.xml

# Step 4: Medical review markup
python scripts/review.py \
  --narrative narrative.txt \
  --output review_version.txt
```

**Python API:**

```python
from scripts.narrative_generator import NarrativeGenerator
from scripts.validator import NarrativeValidator

# Initialize
generator = NarrativeGenerator()
validator = NarrativeValidator()

# Load case data
import json
with open("case_001.json", "r") as f:
    case = json.load(f)

# Generate narrative
narrative = generator.generate(
    case_data=case,
    format="cioms_i",
    include_meddra=True,
    language="en"
)

# Validate
validation = validator.check(
    narrative=narrative,
    criteria=["cioms_completeness", "temporal_logic", "meddra_accuracy"]
)

if validation.passed:
    with open("final_narrative.txt", "w") as f:
        f.write(narrative.text)
    print("✓ Narrative validated and saved")
else:
    print(f"⚠ Issues: {validation.issues}")
```

## Quality Checklist

**Pre-Generation:**
- [ ] Case ID unique and formatted per SOP
- [ ] Patient age/sex complete
- [ ] Suspect drug(s) clearly identified
- [ ] Adverse event(s) coded with MedDRA PT
- [ ] Dates consistent (no future dates)
- [ ] Reporter information included

**Narrative Content:**
- [ ] All CIOMS I sections present
- [ ] Temporal sequence clear and logical
- [ ] Dechallenge/rechallenge described (if applicable)
- [ ] Lab values with reference ranges
- [ ] Concomitant medications listed
- [ ] Medical history relevant to event
- [ ] Outcome clearly stated
- [ ] Causality assessment justified

**Post-Generation:**
- [ ] MedDRA terms accurate and current
- [ ] No contradictory information
- [ ] Language objective and factual
- [ ] No speculation or opinion (except causality section)
- [ ] Patient identifiers removed or de-identified
- [ ] **CRITICAL**: Medical review completed
- [ ] **CRITICAL**: Causality assessment by qualified physician

## Common Pitfalls

**Completeness Issues:**
- ❌ **Missing dechallenge information** → Cannot assess causality
  - ✅ Always document effect after drug discontinuation

- ❌ **Vague temporal information** → "Recently started" vs. specific dates
  - ✅ Use exact dates when available

- ❌ **Incomplete concomitant medication list** → Alternative causes missed
  - ✅ Include all medications within relevant timeframe

**Medical Accuracy Issues:**
- ❌ **Incorrect MedDRA coding** → Wrong medical concept
  - ✅ Use current MedDRA version; verify with medical reviewer

- ❌ **Confusing correlation with causation** → Temporal = causal
  - ✅ Clearly state "temporally associated" vs. "causally related"

- ❌ **Omitting alternative diagnoses** → Biased toward drug causation
  - ✅ Include all differential diagnoses considered

**Regulatory Issues:**
- ❌ **Opinion in narrative body** → "Clearly caused by drug"
  - ✅ Reserve opinion for causality section; narrative should be factual

- ❌ **Patient identifiers** → HIPAA/privacy violation
  - ✅ De-identify per regulatory requirements

- ❌ **Abbreviations not defined** → Assumes reader knowledge
  - ✅ Spell out on first use in each narrative

## References

Available in `references/` directory:

- `cioms_i_guidelines.pdf` - CIOMS I international reporting standards
- `ich_e2b_specifications.md` - ICH E2B(R3) electronic format details
- `meddra_coding_guide.md` - MedDRA terminology and coding principles
- `who_umc_causality.md` - WHO causality assessment criteria
- `fda_medwatch_guide.md` - FDA Form 3500A instructions
- `gvp_module_vi.md` - EU Good Pharmacovigilance Practices
- `narrative_templates.md` - Example narratives by case type

## Scripts

Located in `scripts/` directory:

- `main.py` - CLI interface for narrative generation
- `narrative_generator.py` - Core narrative composition engine
- `temporal_analyzer.py` - Timeline reconstruction and analysis
- `causality_assessor.py` - Causality evaluation support
- `meddra_integrator.py` - Medical terminology and coding
- `validator.py` - Completeness and quality checks
- `format_converter.py` - Convert between CIOMS, E2B, MedWatch formats
- `batch_processor.py` - Multi-case narrative generation

## Limitations

- **Medical Review Required**: Generates draft only; requires physician review before submission
- **Causality Assessment**: Structures reporter's assessment; does not perform independent causality evaluation
- **MedDRA Version**: Uses installed MedDRA version; may not have latest terms
- **Language**: Optimized for English; other languages may need translation
- **Literature Integration**: Does not automatically search literature for similar cases
- **Signal Detection**: Individual case narratives only; aggregate analysis requires other tools
- **Legal Proceedings**: Not suitable for litigation support or expert witness reports

---

**⚠️ CRITICAL: This tool generates draft narratives for efficiency. All adverse event narratives require review by qualified drug safety physicians before regulatory submission. Causality assessment must be performed by healthcare professionals with access to complete medical records.**
