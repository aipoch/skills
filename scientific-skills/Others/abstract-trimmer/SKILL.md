---
name: abstract-trimmer
description: Compress academic abstracts to meet strict word limits while preserving 
  key information, scientific accuracy, and readability. Supports multiple compression 
  strategies for journal submissions, conference applications, and grant proposals.
allowed-tools: [Read, Write, Bash, Edit]
license: MIT
metadata:
    skill-author: AIPOCH
---

# Abstract Trimmer

## Overview

Precision editing tool that reduces abstract word count through intelligent compression techniques, maintaining scientific rigor while meeting strict journal and conference requirements.

**Key Capabilities:**
- **Smart Compression**: Multiple strategies (aggressive, conservative, balanced)
- **Key Information Preservation**: Retains critical findings and statistics
- **Structural Integrity**: Maintains Background-Methods-Results-Conclusion flow
- **Quantitative Safety**: Protects numbers, P-values, and confidence intervals
- **Batch Processing**: Trim multiple abstracts efficiently
- **Quality Validation**: Post-trim readability and accuracy checks

## When to Use

**✅ Use this skill when:**
- Abstract exceeds journal word limit (e.g., 300 words need to be 250)
- Conference submission has strict character/line limits
- Grant proposal abstract needs to fit page constraints
- Creating short-form summaries from detailed abstracts
- Adapting abstract for different venues (full paper → conference)
- Emergency trimming when submission deadline approaches

**❌ Do NOT use when:**
- Writing abstract from scratch → Use `abstract-summarizer`
- Need to expand short abstract → Use `abstract-summarizer`
- Substantive content changes needed → This only trims, doesn't rewrite
- Translation required → Use `medical-translation`
- Abstract is already within limit → No trimming needed

**Integration:**
- **Upstream**: `abstract-summarizer` (create initial abstract)
- **Downstream**: `journal-matchmaker` (submit to venue), `conference-abstract-adaptor` (format)

## Core Capabilities

### 1. Intelligent Word Reduction

Multiple compression strategies based on context:

```python
from scripts.trimmer import AbstractTrimmer

trimmer = AbstractTrimmer()

# Trim to specific word count
result = trimmer.trim(
    abstract=original_text,
    target_words=250,
    strategy="balanced",  # aggressive, conservative, balanced
    protect=["statistics", "sample_size", "p_values"]  # Never remove
)

print(f"Reduced from {result.original_words} to {result.final_words} words")
print(f"Reduction: {result.reduction_percent}%")
```

**Strategies:**
| Strategy | Approach | Best For |
|----------|----------|----------|
| **Conservative** | Remove filler words, simplify sentences | Minor trims (10-20 words) |
| **Balanced** | Condense phrases, merge sentences | Moderate trims (20-50 words) |
| **Aggressive** | Remove secondary details, abbreviate | Major trims (50+ words) |

### 2. Critical Information Protection

Ensure key scientific content is never lost:

```python
# Configure protected elements
protected = trimmer.configure_protection(
    quantitative=True,      # Numbers, stats, P-values
    primary_outcome=True,   # Main study endpoint
    sample_size=True,       # n= values
    effect_sizes=True,      # Cohen's d, OR, HR
    confidence_intervals=True,  # 95% CI
    safety_data=True        # Adverse events
)

# Trim with protection
trimmed = trimmer.trim_with_protection(
    abstract=abstract,
    target=250,
    protected_elements=protected
)
```

**Never Removes:**
- Sample sizes (n=128)
- P-values (p < 0.001)
- Effect sizes (OR=2.4)
- Confidence intervals (95% CI: [1.2, 3.6])
- Primary endpoint results
- Safety information (adverse events)

### 3. Section-Aware Trimming

Preserve structural integrity of academic abstracts:

```python
# Section-aware trimming
result = trimmer.trim_sections(
    abstract=abstract,
    target_words=250,
    section_limits={
        "background": 40,      # Max words per section
        "objective": 25,
        "methods": 60,
        "results": 100,
        "conclusion": 25
    }
)
```

**Section Priorities:**
1. **Results**: Highest priority (keep all primary findings)
2. **Conclusion**: High priority (main takeaway)
3. **Objective**: Medium priority (study goal)
4. **Methods**: Lower priority (can abbreviate)
5. **Background**: Lowest priority (can condense most)

### 4. Quality Validation

Verify trimmed abstract maintains scientific integrity:

```python
# Validate trimmed output
validation = trimmer.validate(
    original=original_abstract,
    trimmed=trimmed_abstract,
    checks=[
        "quantitative_accuracy",    # Numbers unchanged
        "logical_flow",             # Sentences make sense
        "scientific_rigor",         # No diluted claims
        "readability",              # Still readable
        "word_count"                # Target met
    ]
)

if validation.passed:
    save(trimmed_abstract)
else:
    print(f"Issues found: {validation.issues}")
```

## Common Patterns

### Pattern 1: Journal Submission Over-Limit

**Scenario**: Abstract is 320 words, journal limit is 250.

```bash
# Aggressive trimming for tight constraints
python scripts/main.py \
  --input abstract.txt \
  --target 250 \
  --strategy aggressive \
  --protect statistics \
  --output trimmed.txt
```

**Typical Cuts:**
- Remove: "It is widely known that..." → Background fluff
- Condense: "In this study, we investigated..." → "We investigated..."
- Merge: "We used X. We also used Y." → "We used X and Y."
- Abbreviate: "randomized controlled trial" → "RCT" (after first use)

**Result**: 320 → 248 words (22.5% reduction)

### Pattern 2: Conference Character Limit

**Scenario**: Conference has 1,500 character limit including spaces.

```python
trimmer = AbstractTrimmer()

# Character-based trimming
result = trimmer.trim_by_chars(
    abstract=abstract,
    target_chars=1500,
    strategy="balanced"
)

# Check exact character count
chars = len(result.text)
print(f"Final: {chars}/1500 characters")
```

**Approach:**
- Remove articles ("the", "a") where possible
- Use abbreviations consistently
- Convert "do not" → "don't" (saves 1 char)
- Remove redundant adjectives

### Pattern 3: Adapting Between Venues

**Scenario**: Adapt 350-word full paper abstract for 200-word conference.

```bash
# Create multiple versions
python scripts/main.py --input full_abstract.txt --target 300 --output journal_version.txt
python scripts/main.py --input full_abstract.txt --target 200 --output conference_version.txt
python scripts/main.py --input full_abstract.txt --target 150 --output short_version.txt
```

**Adaptation Strategy:**
- **Journal (300w)**: Full detail, secondary outcomes included
- **Conference (200w)**: Primary outcomes only, minimal methods
- **Short (150w)**: Key result + main conclusion only

### Pattern 4: Emergency Deadline Trimming

**Scenario**: Submission in 30 minutes, abstract is 50 words over.

```python
# Quick conservative trim
result = trimmer.quick_trim(
    abstract=abstract,
    excess_words=50,
    priority="speed"  # Fast processing
)

# One-click validation
if result.safe_to_submit:
    result.copy_to_clipboard()
    print("Ready to paste!")
```

**Emergency Tactics:**
1. Remove all adverbs ("significantly" → implied by P-value)
2. Delete "In conclusion" type phrases
3. Convert passive to active voice (saves words)
4. Remove citation numbers from abstract

## Complete Workflow Example

**From over-limit to submission-ready:**

```bash
# Step 1: Check current word count
python scripts/main.py --input my_abstract.txt --check-only
# Output: Current: 287 words | Target: 250 | Excess: 37

# Step 2: Trim with validation
python scripts/main.py \
  --input my_abstract.txt \
  --target 250 \
  --strategy balanced \
  --protect statistics \
  --validate \
  --output trimmed_abstract.txt

# Step 3: Compare versions
python scripts/compare.py \
  --original my_abstract.txt \
  --trimmed trimmed_abstract.txt \
  --format side-by-side

# Step 4: Final validation
python scripts/validate.py \
  --input trimmed_abstract.txt \
  --checks quantitative,readability \
  --output validation_report.txt
```

**Python API:**

```python
from scripts.trimmer import AbstractTrimmer
from scripts.validator import TrimValidator

# Initialize
trimmer = AbstractTrimmer()
validator = TrimValidator()

# Load abstract
with open("my_abstract.txt", "r") as f:
    original = f.read()

# Trim with protection
result = trimmer.trim(
    abstract=original,
    target_words=250,
    strategy="balanced",
    protect_quantitative=True
)

# Validate
validation = validator.compare(
    original=original,
    trimmed=result.text
)

if validation.quantitative_consistent and validation.readable:
    with open("final.txt", "w") as f:
        f.write(result.text)
    print(f"✓ Trimmed: {result.original_words} → {result.final_words} words")
else:
    print(f"⚠ Issues: {validation.issues}")
```

## Quality Checklist

**Pre-Trimming:**
- [ ] Current word count confirmed
- [ ] Target word count clear
- [ ] Protected elements identified (stats, safety data)
- [ ] Original abstract backed up

**During Trimming:**
- [ ] Protected elements preserved
- [ ] All numbers remain accurate
- [ ] P-values and CIs unchanged
- [ ] Logical flow maintained
- [ ] Grammar remains correct

**Post-Trimming:**
- [ ] Target word count achieved
- [ ] **CRITICAL**: All statistics match original
- [ ] **CRITICAL**: Primary findings preserved
- [ ] **CRITICAL**: Safety information retained (if applicable)
- [ ] Abstract still makes sense standalone
- [ ] Readability score acceptable

**Before Submission:**
- [ ] Proofread trimmed version
- [ ] Verify journal/venue requirements met
- [ ] Check for unintended meaning changes
- [ ] Confirm all authors approve changes

## Common Pitfalls

**Accuracy Issues:**
- ❌ **Rounding numbers during trim** → "67.3%" becomes "67%" → Loss of precision
  - ✅ Never modify numbers; remove surrounding words instead

- ❌ **Removing statistical significance** → "significantly (p<0.001)" → "significantly"
  - ✅ Always preserve P-values and confidence intervals

- ❌ **Losing sample size** → "n=128 patients" → "patients" → Can't assess power
  - ✅ Protect all n= values

**Quality Issues:**
- ❌ **Creating awkward sentences** → "We studied investigated..." → Grammar errors
  - ✅ Always validate readability after trimming

- ❌ **Removing key context** → "improved outcomes" (what outcomes?) → Vague claims
  - ✅ Ensure standalone comprehensibility

- ❌ **Over-trimming methods** → "We used appropriate methods" → Not informative
  - ✅ Keep key methodology descriptors

**Strategy Issues:**
- ❌ **Aggressive trim for minor excess** → 10 words over, aggressive cut
  - ✅ Match strategy to excess amount (conservative for <20 words)

- ❌ **Not validating** → Submit without checking accuracy
  - ✅ Always validate before submission

## References

Available in `references/` directory:

- `compression_strategies.md` - Detailed trimming tactics by section
- `protected_elements.md` - Critical content that must be preserved
- `journal_limits.md` - Word limits by major publishers
- `before_after_examples.md` - Successful trimming examples
- `readability_metrics.md` - Scoring compressed text quality

## Scripts

Located in `scripts/` directory:

- `main.py` - CLI for trimming operations
- `trimmer.py` - Core compression engine
- `validator.py` - Quality checking post-trim
- `comparator.py` - Side-by-side original vs. trimmed
- `batch_processor.py` - Multiple abstracts at once
- `char_counter.py` - Character-based (not word) limits

## Limitations

- **Language**: Optimized for English academic abstracts
- **Content Type**: Designed for structured abstracts (BMRC format); free-form text may need adaptation
- **Quantitative Focus**: Prioritizes preserving numbers; qualitative findings may be over-trimmed
- **No Rewriting**: Only removes/compresses; doesn't rephrase for clarity
- **Domain Specific**: Best for STEM and medical abstracts; humanities may need different strategies
- **Final Review**: Automated trimming requires human validation before submission

---

**✂️ Remember: This tool helps meet word limits, but never sacrifice scientific accuracy. Always validate that trimmed abstracts maintain the integrity of your findings.**
