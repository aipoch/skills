---
name: chemical-storage-sorter
description: Sort chemicals by compatibility for safe laboratory storage. Prevents dangerous reactions by segregating incompatible chemicals (acids, bases, oxidizers, flammables) and provides storage recommendations compliant with safety regulations.
allowed-tools: [Read, Write, Bash, Edit]
license: MIT
metadata:
  skill-author: AIPOCH
---

# Chemical Storage Sorter

Organize laboratory chemicals into safe storage groups based on chemical compatibility and hazard classification. Prevents dangerous reactions by identifying incompatible pairs and providing segregation guidelines compliant with OSHA, NFPA, and institutional safety standards.

**Key Capabilities:**
- **Automatic Chemical Classification**: Categorize chemicals into hazard groups (acids, bases, oxidizers, flammables, toxics)
- **Compatibility Checking**: Identify incompatible chemical pairs that could react dangerously if stored together
- **Storage Grouping**: Automatically sort chemical inventories into safe storage arrangements
- **Safety Warnings**: Generate warnings for incompatible storage combinations and hazardous interactions
- **Regulatory Compliance**: Follow standard chemical segregation rules per OSHA and NFPA guidelines

---

## When to Use

**✅ Use this skill when:**
- Setting up **new laboratory storage** systems and need to organize chemical inventory
- Preparing for **EHS (Environmental Health & Safety) inspections** or compliance audits
- **Relocating or reorganizing** an existing chemical storage area
- **Inventorying** chemicals and checking current storage arrangements for safety issues
- **Onboarding new lab members** and training them on chemical storage safety
- **Investigating chemical incidents** involving improper storage or reactions
- Creating **standard operating procedures** (SOPs) for chemical handling and storage

**❌ Do NOT use when:**
- Dealing with **unknown chemical compositions** or unlabeled containers → Contact EHS for proper identification first
- Needing **specific temperature requirements** for storage → Use specialized temperature monitoring tools
- Handling **radioactive materials** or biohazards → Follow specialized protocols for these materials
- Seeking **disposal instructions** for chemicals → Use `waste-disposal-guide` for disposal procedures
- Requiring **SDS (Safety Data Sheet) lookup** → Use `safety-data-sheet-reader` for detailed chemical information
- Planning **chemical inventory tracking** → Use `lab-inventory-tracker` for quantity and location tracking

**Related Skills:**
- **上游 (Upstream)**: `safety-data-sheet-reader`, `chemical-structure-converter`
- **下游 (Downstream)**: `lab-inventory-tracker`, `waste-disposal-guide`

---

## Integration with Other Skills

**Upstream Skills:**
- `safety-data-sheet-reader`: Retrieve chemical properties and hazard classifications from SDS
- `chemical-structure-converter`: Identify chemical class from structure or name for accurate categorization

**Downstream Skills:**
- `lab-inventory-tracker`: Record storage locations after chemicals are sorted and assigned
- `waste-disposal-guide`: Identify disposal requirements for incompatible chemicals that need to be removed
- `equipment-maintenance-log`: Track safety cabinet inspections and maintenance

**Complete Workflow:**
```
Chemical Inventory → safety-data-sheet-reader → chemical-storage-sorter → lab-inventory-tracker → Safe Storage
```

---

## Core Capabilities

### 1. Chemical Classification by Hazard Group

Automatically classify chemicals into standard hazard categories based on chemical name, formula, or keywords.

```python
from scripts.main import ChemicalStorageSorter

sorter = ChemicalStorageSorter()

# Classify individual chemicals
chemicals = [
    "Hydrochloric acid",
    "Sodium hydroxide",
    "Hydrogen peroxide 30%",
    "Ethanol",
    "Sodium chloride"
]

for chem in chemicals:
    group = sorter.classify_chemical(chem)
    print(f"{chem}: {group}")

# Output:
# Hydrochloric acid: acids
# Sodium hydroxide: bases
# Hydrogen peroxide 30%: oxidizers
# Ethanol: flammables
# Sodium chloride: general
```

**Hazard Groups:**

| Group | Examples | Key Hazards | Storage Requirements |
|-------|----------|-------------|---------------------|
| **Acids** | HCl, H₂SO₄, HNO₃, acetic acid | Corrosive, reactive | Acid cabinet, secondary containment |
| **Bases** | NaOH, KOH, ammonia, amines | Corrosive, caustic | Base cabinet, separate from acids |
| **Oxidizers** | H₂O₂, KMnO₄, nitrates, hypochlorites | Fire/explosion risk | Cool, dry, away from organics |
| **Flammables** | Ethanol, methanol, acetone, hexane | Fire hazard | Flammable storage cabinet |
| **Toxics** | Cyanides, mercury, arsenic compounds | Poison, bioaccumulation | Locked cabinet, limited access |
| **General** | NaCl, PBS, sucrose, glycerol | Low hazard | General storage |

**Classification Keywords:**

| Group | Keywords Triggers |
|-------|------------------|
| **Acids** | acid, hcl, sulfuric, nitric, acetic, citric, formic |
| **Bases** | hydroxide, naoh, koh, ammonia, amine, carbonate |
| **Flammables** | ethanol, methanol, acetone, ether, hexane, toluene, benzene |
| **Oxidizers** | peroxide, permanganate, hypochlorite, nitrate, chlorate, perchlorate |
| **Toxics** | cyanide, mercury, arsenic, lead, cadmium, thallium |

**Best Practices:**
- ✅ **Use full chemical names** for most accurate classification
- ✅ **Include concentration** when relevant (e.g., "hydrogen peroxide 30%" vs "3%")
- ✅ **Check ambiguous chemicals** manually if classification seems incorrect
- ✅ **Update keyword lists** for lab-specific chemicals not in default database

**Common Issues and Solutions:**

**Issue: Chemical not recognized**
- Symptom: Classified as "general" despite being hazardous
- Solution: Use more specific chemical name; add custom keywords for lab-specific compounds

**Issue: Misclassification of similar names**
- Symptom: "Sodium acetate" classified as acid due to "acet" keyword
- Solution: Check classification results; manually override if needed

### 2. Compatibility Checking Between Chemicals

Determine if two chemicals can be safely stored together without risk of dangerous reactions.

```python
from scripts.main import ChemicalStorageSorter

sorter = ChemicalStorageSorter()

# Check specific chemical pairs
pairs_to_check = [
    ("Hydrochloric acid", "Sodium hydroxide"),
    ("Ethanol", "Hydrogen peroxide"),
    ("Sodium chloride", "Potassium chloride"),
    ("Nitric acid", "Acetone")
]

for chem1, chem2 in pairs_to_check:
    compatible, message = sorter.check_compatibility(chem1, chem2)
    status = "✅ Compatible" if compatible else "❌ INCOMPATIBLE"
    print(f"{chem1} + {chem2}: {status}")
    if not compatible:
        print(f"   Warning: {message}")

# Output:
# Hydrochloric acid + Sodium hydroxide: ❌ INCOMPATIBLE
#    Warning: INCOMPATIBLE: acids cannot be stored with bases
# Ethanol + Hydrogen peroxide: ❌ INCOMPATIBLE
#    Warning: INCOMPATIBLE: flammables cannot be stored with oxidizers
# Sodium chloride + Potassium chloride: ✅ Compatible
# Nitric acid + Acetone: ❌ INCOMPATIBLE
```

**Incompatibility Matrix:**

| Chemical Group | Incompatible With | Reaction Risk |
|----------------|------------------|---------------|
| **Acids** | Bases, oxidizers, cyanides, sulfides | Violent neutralization, toxic gas generation |
| **Bases** | Acids, oxidizers, halogenated compounds | Heat generation, decomposition |
| **Oxidizers** | Flammables, acids, bases, reducing agents | Fire, explosion, violent reactions |
| **Flammables** | Oxidizers, acids | Fire, combustion enhancement |
| **Toxics** | Acids, oxidizers | Toxic gas release, increased hazard |

**Best Practices:**
- ✅ **Check all new chemicals** against existing storage before placement
- ✅ **Use minimum 3-foot separation** for incompatible groups
- ✅ **Consider secondary containment** for highly reactive pairs
- ✅ **Document exceptions** with engineering controls in place

**Common Issues and Solutions:**

**Issue: False positive compatibility**
- Symptom: Tool says compatible but chemicals actually react
- Causes: Missing specific incompatibility not in general rules
- Solution: Always consult SDS for specific incompatibilities; use this as first check only

**Issue: Ambiguous compatibility**
- Symptom: "Compatible with precautions" message for borderline cases
- Solution: Err on side of caution; store separately or consult EHS

### 3. Automated Storage Grouping

Sort an entire chemical inventory into safe storage groups based on hazard classifications.

```python
from scripts.main import ChemicalStorageSorter

sorter = ChemicalStorageSorter()

# Example lab inventory
inventory = [
    "Hydrochloric acid (conc.)",
    "Sodium hydroxide pellets",
    "Ethanol 95%",
    "Acetone",
    "Hydrogen peroxide 30%",
    "Potassium permanganate",
    "Sodium chloride",
    "PBS buffer",
    "Glycerol",
    "Sulfuric acid",
    "Ammonium hydroxide",
    "Methanol",
    "Hexane",
    "Mercury(II) chloride"
]

# Sort into storage groups
groups = sorter.sort_chemicals(inventory)

# Display results
for group, chemicals in groups.items():
    if chemicals:
        print(f"\n{group.upper()} STORAGE:")
        for chem in chemicals:
            print(f"  • {chem}")
```

**Storage Group Output:**
```
ACIDS STORAGE:
  • Hydrochloric acid (conc.)
  • Sulfuric acid

BASES STORAGE:
  • Sodium hydroxide pellets
  • Ammonium hydroxide

OXIDIZERS STORAGE:
  • Hydrogen peroxide 30%
  • Potassium permanganate

FLAMMABLES STORAGE:
  • Ethanol 95%
  • Acetone
  • Methanol
  • Hexane

TOXICS STORAGE:
  • Mercury(II) chloride

GENERAL STORAGE:
  • Sodium chloride
  • PBS buffer
  • Glycerol
```

**Best Practices:**
- ✅ **Sort alphabetically within groups** for easier location
- ✅ **Include concentration** in labels for diluted vs concentrated chemicals
- ✅ **Group by frequency of use** - most used chemicals most accessible
- ✅ **Reserve general storage** for the bulk of inventory (typically 60-70%)

**Common Issues and Solutions:**

**Issue: Chemical fits multiple categories**
- Symptom: Chemical has multiple hazards (e.g., concentrated HNO₃ is both acid and oxidizer)
- Solution: Store in most restrictive group (oxidizer cabinet for this example); check all incompatibilities

**Issue: Large inventory processing**
- Symptom: Hundreds of chemicals to sort
- Solution: Process in batches by lab area; export to spreadsheet for manual review

### 4. Storage Plan Generation with Safety Warnings

Generate a complete storage plan with specific warnings and segregation requirements.

```python
from scripts.main import ChemicalStorageSorter

sorter = ChemicalStorageSorter()

# Generate full storage plan
demo_inventory = [
    "HCl (concentrated)",
    "NaOH pellets",
    "Ethanol",
    "Hydrogen peroxide",
    "Sodium cyanide",
    "PBS",
    "Acetone"
]

groups = sorter.sort_chemicals(demo_inventory)
sorter.print_storage_plan(groups)
```

**Sample Output:**
```
============================================================
CHEMICAL STORAGE PLAN
============================================================

ACIDS STORAGE:
----------------------------------------
  • HCl (concentrated)
  ⚠️  Keep away from: bases, oxidizers, cyanides, sulfides

BASES STORAGE:
----------------------------------------
  • NaOH pellets
  ⚠️  Keep away from: acids, oxidizers, halogenated

OXIDIZERS STORAGE:
----------------------------------------
  • Hydrogen peroxide
  ⚠️  Keep away from: flammables, acids, bases, reducing

FLAMMABLES STORAGE:
----------------------------------------
  • Ethanol
  • Acetone
  ⚠️  Keep away from: oxidizers, acids

TOXICS STORAGE:
----------------------------------------
  • Sodium cyanide
  ⚠️  Keep away from: acids, oxidizers

GENERAL STORAGE:
----------------------------------------
  • PBS

============================================================
```

**Storage Requirements by Group:**

| Group | Cabinet Type | Ventilation | Special Requirements |
|-------|-------------|-------------|---------------------|
| **Acids** | Acid cabinet | Fume hood access | Secondary containment, corrosion-resistant |
| **Bases** | Base cabinet | Standard | Keep separate from acids (minimum 3 feet) |
| **Oxidizers** | Standard/oxidizer | Cool, dry location | Away from ignition sources |
| **Flammables** | Flammable cabinet | Explosion-proof | Bonding/grounding for dispensing |
| **Toxics** | Locked cabinet | Standard | Access log, limited quantities |
| **General** | Standard shelving | Standard | Standard lab storage |

**Best Practices:**
- ✅ **Post plan visibly** near storage areas
- ✅ **Update when chemicals are added/removed**
- ✅ **Include emergency contact info** on storage plan
- ✅ **Review quarterly** for accuracy

**Common Issues and Solutions:**

**Issue: Insufficient storage space**
- Symptom: Multiple groups needing same cabinet type
- Solution: Prioritize by hazard level; obtain additional cabinets if needed

**Issue: Chemicals with multiple incompatibilities**
- Symptom: One chemical incompatible with many others
- Solution: Isolate in separate location; consider reducing inventory

### 5. Batch Inventory Processing

Process large chemical inventories from files for comprehensive storage organization.

```python
from scripts.main import ChemicalStorageSorter
import json

def process_inventory_file(file_path: str) -> dict:
    """
    Process chemical inventory from text file.
    
    Expected format: One chemical per line
    """
    sorter = ChemicalStorageSorter()
    
    # Read inventory
    with open(file_path, 'r') as f:
        chemicals = [line.strip() for line in f if line.strip()]
    
    # Sort into groups
    groups = sorter.sort_chemicals(chemicals)
    
    # Calculate statistics
    stats = {
        'total_chemicals': len(chemicals),
        'groups': {group: len(items) for group, items in groups.items() if items},
        'hazardous_chemicals': sum(len(items) for group, items in groups.items() 
                                   if group != 'general' and items)
    }
    
    # Check for incompatibilities within current storage
    incompatibilities = []
    all_groups = list(groups.keys())
    
    for i, group1 in enumerate(all_groups):
        for group2 in all_groups[i+1:]:
            if group2 in sorter.COMPATIBILITY_GROUPS[group1]['incompatible']:
                if groups[group1] and groups[group2]:
                    incompatibilities.append({
                        'group1': group1,
                        'chemicals1': groups[group1],
                        'group2': group2,
                        'chemicals2': groups[group2]
                    })
    
    return {
        'groups': groups,
        'statistics': stats,
        'incompatibilities': incompatibilities
    }

# Example usage
# results = process_inventory_file('lab_inventory.txt')
# print(json.dumps(results, indent=2))
```

**Input File Format:**
```
# lab_inventory.txt
Hydrochloric acid (37%)
Sodium hydroxide
Ethanol (95%)
Acetone
Hydrogen peroxide (30%)
Potassium permanganate
Sodium chloride
Phosphate buffered saline
Glycerol
Sulfuric acid (conc.)
```

**Best Practices:**
- ✅ **Use standardized naming** in inventory files
- ✅ **Include concentrations** for diluted chemicals
- ✅ **Date the inventory** for tracking changes
- ✅ **Archive old versions** for historical reference

**Common Issues and Solutions:**

**Issue: Typos and inconsistent naming**
- Symptom: Same chemical listed multiple ways
- Solution: Standardize naming convention; use CAS numbers for ambiguous cases

**Issue: Concentration variations**
- Symptom: Multiple entries for "ethanol" at different concentrations
- Solution: Include concentration in name; store according to most hazardous concentration

### 6. Custom Classification Rules

Extend the classification system with lab-specific chemicals and custom rules.

```python
from scripts.main import ChemicalStorageSorter

class CustomChemicalSorter(ChemicalStorageSorter):
    """Extended sorter with lab-specific chemicals."""
    
    def __init__(self):
        super().__init__()
        # Add custom chemicals to groups
        self.COMPATIBILITY_GROUPS['acids']['examples'].extend([
            'trifluoroacetic acid',
            'trichloroacetic acid'
        ])
        
        self.COMPATIBILITY_GROUPS['flammables']['examples'].extend([
            'isopropanol',
            'isopropyl alcohol',
            '2-propanol'
        ])
        
        # Add custom keyword mappings
        self.custom_keywords = {
            'acids': ['tfa', 'tca'],
            'flammables': ['ipa', 'propanol']
        }
    
    def classify_chemical(self, name):
        """Override with custom keyword checking."""
        name_lower = name.lower()
        
        # Check custom keywords first
        for group, keywords in self.custom_keywords.items():
            if any(kw in name_lower for kw in keywords):
                return group
        
        # Fall back to parent classification
        return super().classify_chemical(name)

# Use custom sorter
custom_sorter = CustomChemicalSorter()
print(custom_sorter.classify_chemical("TFA"))  # Will classify as acid
print(custom_sorter.classify_chemical("IPA"))  # Will classify as flammable
```

**Best Practices:**
- ✅ **Document custom rules** in lab SOP
- ✅ **Share with all lab members** for consistency
- ✅ **Review periodically** for completeness
- ✅ **Update when new chemicals** are introduced

**Common Issues and Solutions:**

**Issue: Custom rules conflict with defaults**
- Symptom: Chemical classified differently than expected
- Solution: Check rule priority; custom rules should typically override defaults

**Issue: Too many custom chemicals**
- Symptom: Most chemicals need custom classification
- Solution: Update default database instead; contribute improvements upstream

---

## Complete Workflow Example

**From chemical inventory to organized storage:**

```bash
# Step 1: List current chemicals
python scripts/main.py --chemicals "HCl,NaOH,ethanol,acetone,H2O2,PBS"

# Step 2: Check compatibility of specific pair
python scripts/main.py --chemicals "HCl" --check "NaOH"

# Step 3: View storage groups
python scripts/main.py --list-groups

# Step 4: Process full inventory file
python scripts/main.py --chemicals "$(cat inventory.txt | tr '\n' ',')"
```

**Python API Usage:**

```python
from scripts.main import ChemicalStorageSorter

def organize_lab_storage(chemical_inventory: list) -> dict:
    """
    Complete workflow for organizing laboratory chemical storage.
    
    Returns:
        Dictionary with storage groups, warnings, and recommendations
    """
    sorter = ChemicalStorageSorter()
    
    # Sort chemicals into groups
    groups = sorter.sort_chemicals(chemical_inventory)
    
    # Generate storage plan
    print("\n" + "="*60)
    print("LABORATORY CHEMICAL STORAGE ORGANIZATION")
    print("="*60)
    
    sorter.print_storage_plan(groups)
    
    # Identify potential issues
    issues = []
    
    # Check for high-hazard concentrations
    hazardous_chemicals = []
    for group in ['acids', 'bases', 'oxidizers']:
        for chem in groups[group]:
            if 'conc' in chem.lower() or 'concentrated' in chem.lower():
                hazardous_chemicals.append((chem, group))
    
    if hazardous_chemicals:
        issues.append({
            'type': 'concentrated_hazard',
            'chemicals': hazardous_chemicals,
            'recommendation': 'Ensure secondary containment and fume hood access'
        })
    
    # Check storage space distribution
    total_chemicals = len(chemical_inventory)
    general_percentage = len(groups['general']) / total_chemicals * 100
    
    if general_percentage < 50:
        issues.append({
            'type': 'high_hazard_ratio',
            'message': f'Only {general_percentage:.1f}% chemicals are general storage',
            'recommendation': 'Review if all hazardous classifications are necessary'
        })
    
    # Compile results
    results = {
        'storage_groups': groups,
        'statistics': {
            'total_chemicals': total_chemicals,
            'hazardous_chemicals': total_chemicals - len(groups['general']),
            'general_percentage': general_percentage
        },
        'issues': issues,
        'recommendations': [
            'Label all storage cabinets with group names',
            'Post incompatibility matrix near storage area',
            'Schedule quarterly storage inspections',
            'Train all lab members on chemical segregation'
        ]
    }
    
    return results

# Execute workflow
inventory = [
    "Hydrochloric acid (conc.)",
    "Sulfuric acid",
    "Sodium hydroxide",
    "Potassium hydroxide",
    "Ethanol 95%",
    "Methanol",
    "Acetone",
    "Hydrogen peroxide 30%",
    "Nitric acid",
    "Sodium chloride",
    "PBS",
    "Tris buffer",
    "EDTA",
    "Glycerol"
]

results = organize_lab_storage(inventory)

print("\n" + "="*60)
print("SUMMARY")
print("="*60)
print(f"Total chemicals: {results['statistics']['total_chemicals']}")
print(f"Hazardous: {results['statistics']['hazardous_chemicals']}")
print(f"General storage: {results['statistics']['general_percentage']:.1f}%")

if results['issues']:
    print("\n⚠️  Issues identified:")
    for issue in results['issues']:
        print(f"  - {issue['type']}: {issue.get('recommendation', '')}")

print("\n📋 Recommendations:")
for rec in results['recommendations']:
    print(f"  • {rec}")
```

**Expected Output Files:**

```
storage_organization/
├── storage_plan.txt         # Human-readable storage layout
├── chemical_groups.json     # Machine-readable group assignments
├── incompatibilities.csv    # List of incompatible pairs
└── recommendations.md       # Safety recommendations
```

---

## Common Patterns

### Pattern 1: New Lab Setup

**Scenario**: Setting up chemical storage for a new laboratory from scratch.

```json
{
  "setup_type": "new_lab",
  "space": "2 fume hoods, 3 acid cabinets, 2 flammable cabinets",
  "inventory_size": "~200 chemicals expected",
  "special_requirements": [
    "Cell culture focus - many biological buffers",
    "Molecular biology - EtBr, acrylamide",
    "Some organic synthesis - various solvents"
  ],
  "compliance": "OSHA, university EHS"
}
```

**Workflow:**
1. Inventory all chemicals before they arrive
2. Classify each chemical using this tool
3. Assign storage locations based on groups
4. Purchase appropriate cabinets (acid, flammable, etc.)
5. Label all storage areas clearly
6. Train all lab members on the system
7. Post emergency procedures and contact numbers

**Output Example:**
```
New Lab Storage Plan:

CABINET ASSIGNMENTS:
  Acid Cabinet #1: 12 acids
  Acid Cabinet #2: 8 oxidizers (also acids)
  Base Cabinet: 6 bases
  Flammable Cabinet #1: 15 solvents (ethanol, methanol, etc.)
  Flammable Cabinet #2: 8 other flammables
  Toxic Cabinet: 3 chemicals (EtBr, acrylamide, mercury salts)
  General Storage: 148 buffers, salts, reagents

SPACE UTILIZATION:
  Acid cabinets: 20/30 capacity (67%)
  Flammable: 23/40 capacity (58%)
  General: 148/200 capacity (74%)

RECOMMENDATION: Current space adequate for planned inventory
```

### Pattern 2: Safety Inspection Preparation

**Scenario**: Preparing for annual EHS safety inspection.

```json
{
  "inspection_type": "annual_ehs",
  "focus_areas": [
    "Chemical segregation compliance",
    "Incompatible storage checks",
    "Labeling and signage",
    "Secondary containment"
  ],
  "documentation_required": [
    "Chemical inventory",
    "Storage plan",
    "Incompatibility records"
  ]
}
```

**Workflow:**
1. Run full inventory through storage sorter
2. Check current storage against recommendations
3. Identify any incompatibilities in current arrangement
4. Move misplaced chemicals to proper storage
5. Update storage plan documentation
6. Print and post current storage map
7. Verify all cabinets properly labeled
8. Check secondary containment systems

**Output Example:**
```
Pre-Inspection Report:

✅ COMPLIANT STORAGE: 187/195 chemicals (95.9%)

⚠️  ISSUES IDENTIFIED:
  1. Acetic acid (glacial) stored with general chemicals
     → Move to acid cabinet
  2. Hydrogen peroxide near ethanol shelf
     → Move to oxidizer section
  3. Missing secondary containment for HCl
     → Add acid tray

📋 DOCUMENTATION READY:
  ✓ Chemical inventory (195 items)
  ✓ Storage plan (updated 2026-02-09)
  ✓ Incompatibility matrix (posted)
  ✓ Emergency contacts (current)

INSPECTION READINESS: 95% (2 chemicals need moving)
```

### Pattern 3: Chemical Relocation

**Scenario**: Moving chemicals to a new location or different lab.

```json
{
  "relocation_type": "lab_move",
  "from": "Building A, Room 301",
  "to": "Building B, Room 205",
  "chemicals_to_move": 150,
  "special_considerations": [
    "Some chemicals expire soon",
    "Unknown origin of 5 chemicals",
    "Need to dispose of 20 chemicals"
  ]
}
```

**Workflow:**
1. Inventory all chemicals at current location
2. Classify and sort all chemicals
3. Identify chemicals for disposal (expired, unknown, unneeded)
4. Plan packing by storage groups (pack together)
5. Ensure proper segregation during transport
6. Design storage layout at new location
7. Unpack directly into appropriate storage
8. Update inventory with new locations

**Output Example:**
```
Relocation Plan:

CHEMICALS TO MOVE: 130 items
  - Acids: 8 (pack together, upright)
  - Bases: 5 (pack together, separate from acids)
  - Flammables: 22 (DOT-approved containers)
  - Oxidizers: 6 (separate transport)
  - Toxics: 2 (locked container, manifest required)
  - General: 87 (standard boxes)

CHEMICALS TO DISPOSE: 20 items
  - Expired: 12
  - Unknown: 5
  - Unneeded: 3
  → Schedule waste pickup before move

PACKING SEQUENCE:
  Day 1: Dispose of waste chemicals
  Day 2: Pack general chemicals
  Day 3: Pack flammables and toxics
  Day 4: Transport and unpack
  Day 5: Final inventory at new location
```

### Pattern 4: Training New Lab Members

**Scenario**: Training new graduate students or technicians on chemical safety.

```json
{
  "training_type": "new_member_safety",
  "trainees": 3,
  "duration": "2 hours",
  "topics": [
    "Chemical hazard recognition",
    "Storage segregation rules",
    "Emergency procedures",
    "Finding chemicals in lab"
  ]
}
```

**Workflow:**
1. Introduce chemical hazard groups using this tool
2. Show real examples from lab inventory
3. Demonstrate compatibility checking
4. Practice classifying unknown chemicals
5. Tour actual storage areas
6. Quiz on incompatible pairs
7. Provide storage plan as reference
8. Document training completion

**Output Example:**
```
Training Session: Chemical Storage Safety

DEMONSTRATION EXAMPLES:
  1. Show classification: "ethanol" → flammable
  2. Show incompatibility: HCl + NaOH → violent reaction
  3. Show safe storage: PBS + NaCl → general storage together

INTERACTIVE QUIZ:
  Q: Can you store acetone near hydrogen peroxide?
  A: No - flammable + oxidizer = fire risk ✅
  
  Q: Where should concentrated HCl go?
  A: Acid cabinet with secondary containment ✅

HANDOUTS PROVIDED:
  ✓ Storage plan (current)
  ✓ Incompatibility matrix
  ✓ Emergency contact card
  ✓ SDS access instructions

TRAINING COMPLETE: 3/3 trainees passed quiz (100%)
```

---

## Quality Checklist

**Pre-Organization:**
- [ ] **CRITICAL**: Ensure all chemical containers are properly labeled
- [ ] Obtain complete chemical inventory (including concentrations)
- [ ] Review SDS for chemicals with unclear classifications
- [ ] Measure available storage space (cabinets, shelves)
- [ ] Identify existing storage equipment (acid cabinets, flammable cabinets)
- [ ] Check for expired or unneeded chemicals to dispose
- [ ] Verify emergency equipment availability (eyewash, shower, spill kit)
- [ ] Review institutional EHS requirements and restrictions

**During Classification:**
- [ ] Classify each chemical using full chemical name
- [ ] Note concentrations for diluted vs concentrated forms
- [ ] **CRITICAL**: Verify classification of borderline chemicals manually
- [ ] Check all chemicals with multiple hazards (e.g., oxidizing acid)
- [ ] Document any custom classifications or exceptions
- [ ] Flag chemicals requiring special storage (temperature, light-sensitive)
- [ ] Identify chemicals needing secondary containment
- [ ] Note any chemicals with expiration dates

**Storage Assignment:**
- [ ] **CRITICAL**: Ensure incompatible groups are physically separated (minimum 3 feet)
- [ ] Verify adequate space in each storage category
- [ ] Place most hazardous chemicals in most secure locations
- [ ] Ensure frequently used chemicals are easily accessible
- [ ] Check that cabinet ventilation is appropriate for contents
- [ ] Verify flammable cabinet is properly grounded
- [ ] Ensure acid cabinet has corrosion-resistant construction
- [ ] Confirm toxic chemicals are in locked storage

**Post-Organization Verification:**
- [ ] **CRITICAL**: Check no incompatible chemicals stored together
- [ ] Verify all containers properly labeled with chemical name and hazards
- [ ] Confirm storage plan is posted near chemical area
- [ ] Check emergency procedures are posted and visible
- [ ] Verify spill kits are appropriate for stored chemicals
- [ ] Ensure SDS binder is accessible and current
- [ ] Test that all lab members can locate chemicals easily
- [ ] Schedule first quarterly inspection

**Documentation:**
- [ ] **CRITICAL**: Update chemical inventory with new storage locations
- [ ] Document any exceptions to standard storage rules
- [ ] Record training completion for all lab members
- [ ] File storage plan with lab notebook or ELN
- [ ] Share storage map with EHS coordinator
- [ ] Set calendar reminder for next inspection
- [ ] Archive old storage plans for reference
- [ ] Update lab SOP with storage procedures

---

## Common Pitfalls

**Classification Errors:**
- ❌ **Assuming dilute = safe** → Even dilute acids/bases need proper storage
  - ✅ Classify by chemical identity, not just concentration
  
- ❌ **Ignoring chemical name keywords** → Missing hazards in complex names
  - ✅ Check for multiple hazard indicators in chemical names
  
- ❌ **Not considering mixtures** → Commercial reagents may have multiple components
  - ✅ Check SDS for mixture compositions and store by most hazardous component
  
- ❌ **Classifying by use rather than hazard** → Storing buffer salts with acids
  - ✅ Always use hazard-based classification for storage

**Storage Arrangement Errors:**
- ❌ **Inadequate separation** → 6-inch gap instead of 3-foot minimum
  - ✅ Use physical barriers (cabinets) for incompatible groups
  
- ❌ **Storing by alphabetical order** → Acetic acid next to acetone
  - ✅ Always prioritize chemical compatibility over alphabetical order
  
- ❌ **Ignoring spill containment** → No secondary containment for liquids
  - ✅ Use trays or bunds for liquid chemicals, especially corrosives
  
- ❌ **Overcrowding cabinets** → Blocking access to emergency equipment
  - ✅ Maintain clear access to all chemicals and safety equipment

**Documentation Errors:**
- ❌ **Outdated storage plans** → Chemicals moved but map not updated
  - ✅ Update storage documentation whenever chemicals are relocated
  
- ❌ **Missing hazard warnings** → No incompatibility matrix posted
  - ✅ Post storage plan with clear hazard warnings
  
- ❌ **No training records** → Cannot prove safety training occurred
  - ✅ Document all safety training with signatures
  
- ❌ **Incomplete inventories** → Missing chemicals from tracking system
  - ✅ Maintain complete, up-to-date chemical inventory

**Operational Errors:**
- ❌ **Using food containers** → Chemicals stored in drink bottles
  - ✅ Use only appropriate chemical storage containers
  
- ❌ **No expiration monitoring** → Old peroxides or other degradables
  - ✅ Track expiration dates; dispose of expired chemicals promptly
  
- ❌ **Improper labeling** → Abbreviations or formulas only
  - ✅ Use full chemical names plus hazard symbols
  
- ❌ **Blocking access** → Storage in front of eyewash or shower
  - ✅ Maintain 3-foot clearance around all safety equipment

---

## Troubleshooting

**Problem: Chemical cannot be classified**
- Symptoms: Tool returns "general" for obviously hazardous chemical
- Causes:
  - Chemical name not in keyword database
  - Unusual or proprietary chemical name
  - Mixture with complex name
- Solutions:
  - Check SDS for proper chemical name and hazards
  - Use CAS number to look up chemical class
  - Consult with EHS for unusual chemicals
  - Add custom classification rule for lab-specific chemicals

**Problem: Too many "incompatible" pairs identified**
- Symptoms: Hundreds of incompatibilities in a small lab
- Causes:
  - Overly broad incompatibility rules
  - Chemicals already properly separated but flagged
  - Concentration not considered (dilute vs concentrated)
- Solutions:
  - Focus on actual storage arrangements, not theoretical incompatibilities
  - Check if chemicals are already properly segregated
  - Consider concentration exemptions (very dilute solutions)
  - Prioritize by hazard severity

**Problem: Storage space insufficient**
- Symptoms: More chemicals than available cabinet space
- Causes:
  - Inventory has grown over time
  - Improper disposal of old chemicals
  - Over-purchasing of chemicals
- Solutions:
  - Dispose of expired or unneeded chemicals
  - Share chemicals between labs when possible
  - Request additional storage equipment from EHS
  - Implement "just-in-time" purchasing for expensive chemicals
  - Consider chemical inventory reduction program

**Problem: Lab members resist new storage system**
- Symptoms: Chemicals found in wrong locations after reorganization
- Causes:
  - Inadequate training
  - System too complex
  - Old habits hard to break
- Solutions:
  - Provide clear, hands-on training
  - Make storage locations intuitive and convenient
  - Post visible storage maps at point of use
  - Gentle reminders and positive reinforcement
  - Regular audits with feedback

**Problem: Chemical reactions in storage**
- Symptoms: Evidence of reaction (discoloration, gas, heat, fumes)
- Causes:
  - Incompatible chemicals stored together
  - Degradation of unstable chemicals
  - Contamination during storage
- Solutions:
  - **Immediately evacuate area if fumes or heat**
  - Contact EHS for safe cleanup
  - Review storage arrangements to prevent recurrence
  - Check for other potentially affected chemicals
  - Document incident and lessons learned

**Problem: Cannot find chemical when needed**
- Symptoms: Chemical in inventory but not in expected location
- Causes:
  - Chemical moved but inventory not updated
  - Mislabeling or unclear labels
  - Inconsistent naming (acetic acid vs ethanoic acid)
- Solutions:
  - Update inventory immediately when chemicals are moved
  - Use standardized, full chemical names
  - Implement barcode or QR code tracking
  - Keep storage plan current and accessible
  - Regular inventory reconciliation

---

## References

Available in `references/` directory:

- (No reference files currently available for this skill)

**External Resources:**
- OSHA Chemical Storage Guidelines: https://www.osha.gov/chemical-storage
- NFPA 45: Fire Protection for Laboratories Using Chemicals
- Prudent Practices in the Laboratory (National Research Council)
- SDS Search (MSDSOnline): https://www.msdsonline.com

---

## Scripts

Located in `scripts/` directory:

- `main.py` - Chemical classification and storage sorting engine

---

## Chemical Storage Quick Reference

**General Rules:**
1. **Separate incompatible chemicals** by at least 3 feet or physical barrier
2. **Store acids and bases** in separate cabinets
3. **Keep oxidizers away** from flammables and organics
4. **Lock toxic chemicals** and limit access
5. **Use secondary containment** for liquid corrosives
6. **Label all containers** with chemical name and hazards
7. **Never store chemicals** in food containers or near food areas
8. **Maintain access** to safety equipment (eyewash, shower, exits)

**Emergency Contacts:**
- Fire: 911
- Poison Control: 1-800-222-1222
- Campus EHS: [Insert local number]
- Chemical Spill Hotline: [Insert local number]

## Parameters

| Parameter | Type | Default | Required | Description |
|-----------|------|---------|----------|-------------|
| `--chemicals`, `-c` | string | - | No | Comma-separated chemical list |
| `--check` | string | - | No | Check compatibility with another chemical |
| `--list-groups`, `-l` | flag | - | No | List storage groups |

## Usage

### Basic Usage

```bash
# Sort list of chemicals
python scripts/main.py --chemicals "HCl,NaOH,ethanol,H2O2"

# Check compatibility between two chemicals
python scripts/main.py --chemicals "HCl" --check "NaOH"

# List all storage groups
python scripts/main.py --list-groups
```

## Risk Assessment

| Risk Indicator | Assessment | Level |
|----------------|------------|-------|
| Code Execution | Python script executed locally | Low |
| Network Access | No external API calls | Low |
| File System Access | No file access | Low |
| Data Exposure | No sensitive data | Low |
| Safety Risk | Provides chemical safety guidance | Medium |

## Security Checklist

- [x] No hardcoded credentials or API keys
- [x] No file system access
- [x] Input validation for chemical names
- [x] Output does not expose sensitive information
- [x] Error messages sanitized
- [x] Script execution in sandboxed environment

## Prerequisites

```bash
# Python 3.7+
# No additional packages required (uses standard library)
```

## Evaluation Criteria

### Success Metrics
- [x] Successfully classifies chemicals into storage groups
- [x] Identifies incompatible chemical pairs
- [x] Provides storage recommendations
- [x] Lists all available storage groups

### Test Cases
1. **Chemical List**: Input list → Sorted by compatibility groups
2. **Compatibility Check**: Two chemicals → Compatible/Incompatible result
3. **Unknown Chemical**: Unrecognized name → General group assignment

## Lifecycle Status

- **Current Stage**: Active
- **Next Review Date**: 2026-03-09
- **Known Issues**: None
- **Planned Improvements**:
  - Expand chemical database
  - Add SDS integration
  - Support for custom storage rules

---

**Last Updated**: 2026-02-09  
**Skill ID**: 184  
**Version**: 2.0 (K-Dense Standard)
