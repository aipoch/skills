---
name: acronym-unpacker
description: Disambiguate medical acronyms and abbreviations with context-aware full 
  form lookup. Resolves ambiguous abbreviations (e.g., "PID" = Pelvic Inflammatory 
  Disease vs. Prolapsed Intervertebral Disc) based on clinical specialty, document 
  context, and usage patterns.
allowed-tools: [Read, Write, Bash, Edit]
license: MIT
metadata:
    skill-author: AIPOCH
---

# Acronym Unpacker

## Overview

Intelligent medical abbreviation disambiguation tool that resolves ambiguous acronyms using clinical context, specialty-specific knowledge, and document-level semantic analysis.

**Key Capabilities:**
- **Context-Aware Disambiguation**: Uses clinical specialty to rank expansions
- **Semantic Analysis**: Analyzes surrounding text for contextual clues
- **Frequency-Based Ranking**: Prioritizes common usage patterns
- **Multi-Specialty Support**: Covers medicine, nursing, pharmacy, and research
- **Batch Processing**: Expand acronyms in entire documents
- **Learning System**: Improves accuracy with usage feedback

## When to Use

**✅ Use this skill when:**
- Reading clinical notes filled with unfamiliar abbreviations
- Writing patient education materials requiring full terms
- Preparing manuscripts for non-specialist audiences
- Reviewing multi-disciplinary case reports
- Onboarding new staff to specialty-specific terminology
- Translating medical documents for lay audiences
- Auditing documentation for unclear abbreviations

**❌ Do NOT use when:**
- Acronyms are already defined in document → Use find/replace
- Standardizing abbreviations across organization → Use `abbreviation-standardizer`
- Creating abbreviation lists for glossaries → Use `glossary-generator`
- Non-medical technical documents → Use general NLP tools
- Legal or regulatory documents → Use specialized legal tools

**Integration:**
- **Upstream**: `medical-translation` (preprocess foreign documents)
- **Downstream**: `lay-summary-gen` (patient-friendly text), `medical-email-polisher` (professional communication)

## Core Capabilities

### 1. Context-Aware Expansion

Resolve ambiguous acronyms using clinical context:

```python
from scripts.acronym_resolver import AcronymResolver

resolver = AcronymResolver()

# Single acronym with context
result = resolver.resolve(
    acronym="PID",
    context="gynecology",  # or document text
    source_specialty="obstetrics"
)

# Returns ranked options:
# 1. Pelvic Inflammatory Disease (confidence: 0.92)
# 2. Prolapsed Intervertebral Disc (confidence: 0.05)
# 3. Photoionization Detector (confidence: 0.03)
```

**Context Signals:**
| Signal | Examples | Impact |
|--------|----------|--------|
| **Specialty** | Cardiology vs. Orthopedics | Primary filter |
| **Surrounding Terms** | "uterus", "cervix" near "PID" | Boosts gynecological meaning |
| **Document Type** | Discharge summary vs. Research paper | Different abbreviation patterns |
| **Frequency** | Common usage in corpus | Prioritizes standard meanings |

### 2. Document-Level Analysis

Process entire documents to resolve acronyms consistently:

```python
# Batch process document
document = resolver.process_document(
    text=clinical_note,
    strategy="contextual",  # or "first_definition", "interactive"
    min_confidence=0.70
)

# Output shows all expansions with locations
for acronym, expansions in document.acronyms.items():
    print(f"{acronym}: {expansions[0].full_form} (confidence: {expansions[0].confidence})")
```

**Strategies:**
- **First Definition**: Uses first explicit definition in document
- **Contextual**: Analyzes surrounding text for clues
- **Interactive**: Prompts user for ambiguous cases
- **Conservative**: Only high-confidence expansions

### 3. Specialty-Specific Databases

Access domain-specific abbreviation knowledge:

```python
# Load specialty database
resolver.load_specialty("cardiology")

# Now prioritizes cardiology meanings
result = resolver.resolve("MI")  # Returns Myocardial Infarction first

# Available specialties
databases = resolver.list_specialties()
# ['cardiology', 'oncology', 'neurology', 'pediatrics', ...]
```

**Supported Specialties:**
- Cardiology
- Oncology/Hematology
- Neurology/Neurosurgery
- Obstetrics/Gynecology
- Pediatrics
- Emergency Medicine
- Radiology
- Pathology
- Pharmacology
- Nursing

### 4. Feedback and Learning

Improve accuracy through usage feedback:

```python
# Provide feedback on resolution
resolver.provide_feedback(
    acronym="CABG",
    suggested_expansion="Coronary Artery Bypass Grafting",
    context="cardiology",
    was_correct=True
)

# System learns and updates confidence scores
# Future resolutions benefit from feedback
```

## Common Patterns

### Pattern 1: Clinical Note Interpretation

**Scenario**: Reading discharge summary with many abbreviations.

```python
note = """
Pt admitted with SOB, found to have PE. Started on AC. 
Echo showed RV strain. CT chest confirmed bilateral PE.
"""

# Expand all acronyms
expanded = resolver.expand_document(
    text=note,
    specialty="cardiology",
    output_format="annotated"
)

# Output:
# "Patient admitted with Shortness of Breath, found to have Pulmonary Embolism. 
#  Started on Anticoagulation. Echocardiogram showed Right Ventricular strain. 
#  CT chest confirmed bilateral Pulmonary Embolism."
```

**Key Acronyms Resolved:**
- SOB → Shortness of Breath (not Son of a B...)
- PE → Pulmonary Embolism (not Physical Education)
- AC → Anticoagulation (not Air Conditioning)
- RV → Right Ventricular (not Recreational Vehicle)
- CT → Computed Tomography (not Connecticut)

### Pattern 2: Manuscript Preparation

**Scenario**: Writing review article for general medical audience.

```bash
# Check if acronyms need expansion
python scripts/main.py \
  --input manuscript.txt \
  --check-clarity \
  --target-audience general \
  --output flagged_acronyms.txt

# Expand for clarity
python scripts/main.py \
  --input manuscript.txt \
  --expand-all \
  --first-use-only \
  --output expanded_manuscript.txt
```

**Guidelines:**
- Expand on first use in document
- Keep acronym in parentheses for reference
- Don't expand universally known terms (DNA, RNA, CT, MRI)
- Consider journal-specific guidelines

### Pattern 3: Patient Education Materials

**Scenario**: Converting clinical summary for patient understanding.

```python
# Patient-friendly expansion
patient_version = resolver.expand_for_patient(
    text=clinical_summary,
    reading_level=8,  # 8th grade reading level
    explain_jargon=True  # Add brief explanations
)

# Output includes:
# "You have COPD (Chronic Obstructive Pulmonary Disease - a lung condition 
#  that makes breathing difficult)..."
```

**Patient-Friendly Features:**
- Always expand medical acronyms
- Add brief parenthetical explanations
- Avoid nested acronyms (expanded form shouldn't contain more acronyms)
- Use lay terms when available ("heart attack" vs. "myocardial infarction")

### Pattern 4: Multi-Disciplinary Team Communication

**Scenario**: Tumor board with surgeons, oncologists, radiologists.

```python
# Multi-specialty context
resolver.set_context("oncology")

# Process pathology report
pathology_text = "Pt with NSCLC, s/p LLL. Path: ADC, G2, PD-L1 TPS 75%."

expanded = resolver.expand_with_specialty_notes(
    text=pathology_text,
    specialties=["pathology", "oncology", "surgery"],
    show_alternatives=True
)

# Output shows specialty-specific notes:
# "Patient with Non-Small Cell Lung Cancer [oncology: most common lung cancer type], 
#  status post Left Lower Lobectomy [surgery: removal of lower left lung lobe]. 
#  Pathology: Adenocarcinoma [pathology: cancer subtype], Grade 2 [pathology: 
#  moderately differentiated], PD-L1 Tumor Proportion Score 75% [oncology: 
#  high expression, candidate for immunotherapy]."
```

## Complete Workflow Example

**From cryptic note to readable text:**

```bash
# Step 1: Process clinical note
python scripts/main.py \
  --input note.txt \
  --specialty emergency \
  --strategy contextual \
  --output expanded_note.txt

# Step 2: Review ambiguous cases
python scripts/main.py \
  --input note.txt \
  --list-ambiguous \
  --confidence-threshold 0.70 \
  --output ambiguous_acronyms.txt

# Step 3: Generate glossary
python scripts/main.py \
  --input note.txt \
  --generate-glossary \
  --output glossary.md
```

**Python API:**

```python
from scripts.acronym_resolver import AcronymResolver
from scripts.document_processor import DocumentProcessor

# Initialize
resolver = AcronymResolver()
processor = DocumentProcessor()

# Load document
with open("clinical_note.txt", "r") as f:
    text = f.read()

# Process with high confidence threshold
results = resolver.resolve_document(
    text=text,
    specialty="cardiology",
    min_confidence=0.80
)

# Generate annotated version
annotated = processor.annotate_acronyms(
    text=text,
    resolutions=results,
    format="markdown"
)

# Save
with open("annotated_note.md", "w") as f:
    f.write(annotated)

# Review low-confidence cases
low_conf = [r for r in results if r.confidence < 0.80]
if low_conf:
    print(f"⚠ {len(low_conf)} acronyms need manual review")
    for case in low_conf:
        print(f"  - {case.acronym}: {case.top_candidates}")
```

## Quality Checklist

**Pre-Processing:**
- [ ] Document specialty identified
- [ ] Target audience clear (physician vs. patient vs. researcher)
- [ ] Ambiguous acronyms flagged for review

**During Resolution:**
- [ ] Context clues properly weighted
- [ ] Specialty-appropriate meanings prioritized
- [ ] Confidence scores calculated
- [ ] Multiple candidates ranked when ambiguous

**Post-Processing:**
- [ ] All expansions medically accurate
- [ ] Context appropriate for specialty
- [ ] Patient materials use lay terms
- [ ] No nested acronyms in expansions
- [ ] Low-confidence cases flagged for review

**Before Final Use:**
- [ ] **CRITICAL**: Review all high-stakes expansions (medication names, diagnoses)
- [ ] Verify specialty-specific meanings
- [ ] Check for false friends (same acronym, different meanings)
- [ ] Ensure patient materials are understandable

## Common Pitfalls

**Context Issues:**
- ❌ **Ignoring specialty context** → "CA" = Cancer (oncology) vs. Calcium (lab)
  - ✅ Always specify clinical specialty
  
- ❌ **Missing document-level context** → First definition establishes meaning
  - ✅ Use "first definition" strategy when explicit definitions present

**Accuracy Issues:**
- ❌ **Accepting low-confidence expansions** → Wrong diagnosis or medication
  - ✅ Set minimum confidence threshold (recommend 0.70+)
  - ✅ Flag low-confidence cases for manual review

- ❌ **Over-expanding common terms** → Expanding "DNA", "RNA", "CT" unnecessarily
  - ✅ Maintain exception list for universally understood acronyms

**Communication Issues:**
- ❌ **Nested acronyms** → Expansion contains more acronyms
  - ✅ "CABG" → "Coronary Artery Bypass Grafting" (not "CABG surgery")
  
- ❌ **Over-expanding in technical documents** → Annoying expert readers
  - ✅ Consider audience; expand less for specialist journals

## References

Available in `references/` directory:

- `medical_abbreviations_db.md` - Comprehensive abbreviation database
- `specialty_specific_terms.md` - Field-specific acronym mappings
- `context_clues_guide.md` - Semantic analysis for disambiguation
- `patient_friendly_terms.md` - Lay language equivalents
- `false_friends_list.md` - Acronyms with multiple common meanings
- `journal_specific_guidelines.md` - Expansion rules by publisher

## Scripts

Located in `scripts/` directory:

- `main.py` - CLI interface for acronym resolution
- `acronym_resolver.py` - Core disambiguation engine
- `document_processor.py` - Batch document processing
- `context_analyzer.py` - Semantic context extraction
- `specialty_loader.py` - Domain-specific database management
- `feedback_system.py` - User feedback and learning
- `glossary_generator.py` - Create abbreviation glossaries

## Limitations

- **Language**: Optimized for English medical terminology
- **Novel Acronyms**: May miss very new or institution-specific abbreviations
- **Context Dependency**: Accuracy depends on sufficient contextual text
- **Specialty Boundaries**: May struggle with multi-specialty documents
- **Patient Materials**: Requires review for appropriateness
- **Legal/Regulatory**: Not suitable for legal document interpretation

---

**🔤 Remember: Context is key. The same acronym can mean very different things in different specialties. Always verify high-stakes expansions (diagnoses, medications) in clinical contexts.**
