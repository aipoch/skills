---
name: 3d-molecule-ray-tracer
description: Generate photorealistic rendering scripts for PyMOL and UCSF ChimeraX 
  to create publication-quality molecular visualizations. Supports ray-tracing, 
  depth of field, ambient occlusion, and cinematic lighting for journal covers 
  and high-impact figures. Use for creating visually compelling protein structure 
  images; not for routine structural analysis or measurements.
allowed-tools: [Read, Write, Bash, Edit]
license: MIT
metadata:
    skill-author: AIPOCH
---

# 3D Molecule Ray Tracer

## Overview

Advanced molecular visualization tool that generates professional-grade rendering scripts with cinematic effects for creating publication-quality and cover-worthy molecular images.

**Key Capabilities:**
- **Multi-Software Support**: Generate scripts for PyMOL and UCSF ChimeraX
- **Photorealistic Rendering**: Ray-tracing, depth of field, ambient occlusion
- **Cinematic Lighting**: Studio, outdoor, and dramatic lighting presets
- **Publication Presets**: Pre-configured settings for journals, covers, and presentations
- **Customizable Scenes**: Fine control over camera, materials, and atmosphere

## When to Use

**✅ Use this skill when:**
- Creating cover images for high-impact journals (Nature, Science, Cell)
- Generating publication-quality figures for manuscripts
- Preparing visually compelling slides for conference presentations
- Producing marketing materials or grant proposal visuals
- Creating molecular animations with professional lighting
- Teaching materials requiring engaging 3D visualizations

**❌ Do NOT use when:**
- Routine structural analysis or measurements needed → Use standard PyMOL/ChimeraX GUI
- Quick preview of protein structure → Use web-based viewers (Mol*, NGLView)
- Comparing multiple structures quantitatively → Use `protein-struct-viz` for structural alignment
- Generating molecular dynamics trajectories → Use `molecular-dynamics-viz`
- Creating simple 2D diagrams → Use `chemical-structure-converter`
- Need real-time interactive exploration → Use standalone PyMOL/ChimeraX

**Related Skills:**
- **上游**: `protein-struct-viz` (structure preparation), `protein-docking-configurator` (pose optimization)
- **下游**: `multi-panel-figure-assembler` (combining images), `journal-cover-prompter` (AI cover design)

## Integration with Other Skills

**Upstream Skills:**
- `protein-struct-viz`: Prepare and optimize protein structures before rendering
- `protein-docking-configurator`: Generate ligand-protein complexes for visualization
- `chemical-structure-converter`: Convert SMILES to 3D structures for small molecules
- `alphafold-database`: Download predicted structures for visualization

**Downstream Skills:**
- `multi-panel-figure-assembler`: Combine multiple renderings into composite figures
- `journal-cover-prompter`: Design AI-generated covers incorporating molecular images
- `dpi-upscaler-checker`: Enhance resolution of rendered images if needed
- `graphical-abstract-wizard`: Create comprehensive graphical abstracts

**Complete Workflow:**
```
Protein Structure Viz (prepare structure) → 
  3D Molecule Ray Tracer (this skill, generate render scripts) → 
    Multi-panel Figure Assembler (combine images) → 
      Journal Cover Prompter (AI-enhanced cover design)
```

## Core Capabilities

### 1. Multi-Software Script Generation

Generate rendering scripts for both major molecular visualization platforms:

```python
from scripts.renderer import MoleculeRenderer

renderer = MoleculeRenderer()

# Generate PyMOL script
pymol_script = renderer.generate_pymol_script(
    pdb_file="protein.pdb",
    preset="cover",
    output="render.pml"
)

# Generate ChimeraX script
chimerax_script = renderer.generate_chimerax_script(
    pdb_file="protein.pdb",
    preset="cinematic",
    output="render.cxc"
)
```

**Supported Software:**
| Software | Best For | Features |
|----------|----------|----------|
| **PyMOL** | Traditional rendering, ease of use | Ray tracing, shadows, AO |
| **ChimeraX** | Modern effects, large structures | PBR lighting, ambient occlusion, VR |

**Parameters:**
| Parameter | Type | Required | Description | Options |
|-----------|------|----------|-------------|---------|
| `pdb_file` | str | Yes | PDB file path or 4-letter ID | "1mbn", "path/to/file.pdb" |
| `software` | str | Yes | Target software | "pymol", "chimerax" |
| `preset` | str | No | Rendering preset | "standard", "cover", "publication", "cinematic" |
| `output` | str | No | Output script filename | "render.pml" |

**Best Practices:**
- ✅ Choose PyMOL for simpler scenes and faster setup
- ✅ Choose ChimeraX for advanced lighting and large complexes
- ✅ Always test render at low resolution before final high-res output
- ✅ Save camera position and lighting as separate commands for reproducibility

**Common Issues and Solutions:**

**Issue: Software not found**
- Symptom: "PyMOL command not found" when running script
- Solution: Ensure PyMOL/ChimeraX installed and in PATH; use full path if needed

**Issue: Structure loading errors**
- Symptom: "Error loading PDB file"
- Solution: Verify PDB ID exists; check file format (should be standard PDB)

### 2. Advanced Rendering Effects

Configure photorealistic effects for stunning visuals:

```python
# Configure rendering effects
effects = renderer.configure_effects(
    ray_tracing=True,
    depth_of_field=True,
    dof_focus="residue 50",  # Focus on specific residue
    dof_aperture=2.0,        # Higher = more blur
    ambient_occlusion=True,
    shadows=True,
    fog=0.3,                 # Atmospheric depth
    reflections=True
)

# Apply to scene
script = renderer.apply_effects(script, effects)
```

**Rendering Effects:**
| Effect | Description | Best Use Case |
|--------|-------------|---------------|
| **Ray Tracing** | Physically accurate light simulation | All publication images |
| **Depth of Field** | Focus on subject, blur background | Cover images, highlighting active sites |
| **Ambient Occlusion** | Soft shadows in crevices | Enhancing surface detail |
| **Shadows** | Cast shadows for depth | Dramatic lighting scenes |
| **Fog/Volumetrics** | Atmospheric depth | Large complexes, membranes |
| **Reflections** | Glossy surface highlights | Membrane proteins, lipid bilayers |

**Best Practices:**
- ✅ Enable ray tracing for all publication-quality images
- ✅ Use depth of field to draw attention to key features (active sites, binding pockets)
- ✅ Apply ambient occlusion for surface-heavy representations (surface, mesh)
- ✅ Use fog sparingly; can obscure detail if too dense
- ✅ Match effect intensity to journal style (Science/Nature = dramatic; Cell = clean)

**Common Issues and Solutions:**

**Issue: Render time too long**
- Symptom: "Ray tracing taking 30+ minutes"
- Solution: Reduce resolution for preview; disable expensive effects (fog, reflections); use faster AO settings

**Issue: Depth of field looks unnatural**
- Symptom: "Background blur is distracting or too subtle"
- Solution: Adjust aperture (1.0 = subtle, 3.0 = very blurry); ensure focus point is clearly defined

### 3. Lighting and Atmosphere

Professional lighting setups for different moods:

```python
# Apply lighting preset
lighting = renderer.configure_lighting(
    preset="cinematic",  # Options: default, studio, cinematic, outdoor
    key_light_angle=(45, 45),    # Directional light
    fill_light_intensity=0.3,    # Reduce harsh shadows
    rim_light=True,              # Edge highlighting
    ambient_intensity=0.2        # Base illumination
)

# Customize background
background = renderer.set_background(
    color="white",           # white, black, gradient
    gradient_colors=["#1a365d", "#2b6cb0"],  # For gradient
    transparent=False        # For PNG with alpha
)
```

**Lighting Presets:**
| Preset | Characteristics | Best For |
|--------|----------------|----------|
| **Default** | Even, neutral lighting | Standard publication figures |
| **Studio** | Soft shadows, clean | Methodology figures |
| **Cinematic** | Dramatic, high contrast | Journal covers, presentations |
| **Outdoor** | Natural sunlight simulation | Educational materials |
| **Three-Point** | Professional photo studio setup | Product-style renders |

**Best Practices:**
- ✅ Use studio lighting for clear, unambiguous figures
- ✅ Use cinematic lighting for cover submissions (creates visual interest)
- ✅ Match background to journal preference (Nature = white; Science = neutral)
- ✅ Add rim light for complex structures to separate from background
- ✅ Consider color blindness when choosing molecular color schemes

**Common Issues and Solutions:**

**Issue: Washed out or too dark renders**
- Symptom: "Structure barely visible or overexposed"
- Solution: Adjust ambient intensity and key/fill light balance; use auto-exposure if available

**Issue: Harsh shadows obscure features**
- Symptom: "Active site in deep shadow"
- Solution: Add fill light from opposite direction; reduce key light intensity; use AO instead of hard shadows

### 4. Preset Configurations

Pre-optimized settings for common use cases:

```python
# Use presets for quick setup
script = renderer.apply_preset(
    pdb_file="protein.pdb",
    preset="cover",  # Pre-configured for journal covers
    software="pymol"
)

# Available presets:
# - "standard": Basic ray-tracing, 2400px
# - "cover": High-res (3000px), DOF, AO, enhanced lighting
# - "publication": Journal-standard (2400px), clean, no DOF
# - "cinematic": 4K, dramatic lighting, fog, reflections
```

**Preset Specifications:**
| Preset | Resolution | Ray Trace | DOF | AO | Shadows | Use Case |
|--------|------------|-----------|-----|-----|---------|----------|
| **Standard** | 2400px | ✓ | ✗ | ✗ | ✗ | Quick high-quality |
| **Cover** | 3000px | ✓ | ✓ | ✓ | ✓ | Journal covers |
| **Publication** | 2400px | ✓ | ✗ | ✓ | ✗ | Manuscript figures |
| **Cinematic** | 3840px | ✓ | ✓ | ✓ | ✓ | Presentations |
| **Draft** | 1200px | ✗ | ✗ | ✗ | ✗ | Preview/testing |

**Best Practices:**
- ✅ Start with preset closest to your goal, then customize
- ✅ Use "draft" preset for rapid iteration; switch to "cover" for final
- ✅ Check journal guidelines for resolution requirements (usually 300 DPI at print size)
- ✅ Consider file size: 3000px+ images can be 10MB+ per figure

**Common Issues and Solutions:**

**Issue: Preset settings don't match journal requirements**
- Symptom: "Journal requires 300 DPI at 8 inches = 2400px, but preset is 3000px"
- Solution: Calculate required resolution: (print width in inches) × 300 = pixels needed; adjust preset or use "publication" preset

### 5. Camera Positioning and Composition

Control viewpoint for optimal composition:

```python
# Set camera position
camera = renderer.configure_camera(
    position=(0, 0, -50),      # X, Y, Z coordinates
    look_at="center",          # Or specific residue: "residue 100"
    zoom=1.5,                  # Zoom level
    rotation=(30, 45, 0)       # X, Y, Z rotation in degrees
)

# Apply rule of thirds for composition
renderer.apply_composition_guide(
    style="rule_of_thirds",    # rule_of_thirds, golden_ratio, center
    focal_point="binding_site" # What should be at intersection
)
```

**Composition Techniques:**
| Technique | Description | Best For |
|-----------|-------------|----------|
| **Rule of Thirds** | Place subject at 1/3 intersection | Dynamic, modern look |
| **Center** | Symmetrical, subject in middle | Traditional, formal |
| **Diagonal** | Tilted view along protein axis | Emphasizing length/structure |
| **Close-up** | Fill frame with subject | Active sites, interfaces |
| **Context** | Show protein in membrane/assembly | Large complexes |

**Best Practices:**
- ✅ Orient protein with N-to-C terminus or functional axis as diagonal
- ✅ Position active site or key feature at focal point
- ✅ Leave breathing room (don't crop too tightly)
- ✅ Match orientation to figure layout (portrait vs landscape)
- ✅ Save camera positions as named views for consistency across figures

**Common Issues and Solutions:**

**Issue: Important features cut off or obscured**
- Symptom: "Active site at edge of frame" or "ligand behind protein"
- Solution: Adjust camera distance and rotation; use clipping planes to see interior; split into multiple views

**Issue: Flat, uninteresting composition**
- Symptom: "Looks like default PyMOL view"
- Solution: Add slight rotation (not perfectly front-on); use depth of field; consider unusual angles (from below, extreme close-up)

### 6. Batch Rendering and Automation

Generate multiple views or structures efficiently:

```python
# Batch process multiple structures
structures = ["1mbn", "1hho", "2c0k", "3fgu"]

for pdb in structures:
    script = renderer.generate_pymol_script(
        pdb_file=pdb,
        preset="publication",
        output=f"renders/{pdb}_figure.pml"
    )

# Generate multiple views of same structure
views = ["front", "back", "side", "top", "close_up"]
for view in views:
    script = renderer.generate_view(
        pdb_file="protein.pdb",
        view_name=view,
        preset="publication",
        output=f"renders/{view}.pml"
    )
```

**Best Practices:**
- ✅ Use batch processing for figure supplements (many similar views)
- ✅ Maintain consistent lighting and style across batch
- ✅ Name files descriptively: "protein_ligand_bound_front_3000px.png"
- ✅ Create contact sheet (thumbnail overview) of all renders
- ✅ Render at final resolution only after approving previews

**Common Issues and Solutions:**

**Issue: Inconsistent style across batch**
- Symptom: "Figure 1A and 1B have different lighting"
- Solution: Save complete scene setup (lighting, colors, camera) as template; load template for each new render

## Complete Workflow Example

**From PDB to journal cover:**

```bash
# Step 1: Generate cover-quality render script
python scripts/main.py \
  --software pymol \
  --pdb 1mbn \
  --preset cover \
  --dof-focus "residue 64" \
  --dof-aperture 2.0 \
  --lighting cinematic \
  --style surface \
  --resolution 3000 \
  --output cover_render.pml

# Step 2: Run script in PyMOL
pymol cover_render.pml
# This will: load structure, apply settings, ray trace, save image

# Step 3: Check output
cover_render.png  # 3000px, ray-traced, with DOF and AO
```

**Python API Usage:**

```python
from scripts.renderer import MoleculeRenderer
from scripts.composition import CompositionGuide

# Initialize
renderer = MoleculeRenderer()
composer = CompositionGuide()

# Generate scene
script = renderer.generate_pymol_script(
    pdb_file="complex.pdb",
    preset="cover"
)

# Customize effects
script = renderer.add_depth_of_field(
    script, 
    focus="ligand",
    aperture=2.5
)

script = renderer.configure_lighting(
    script,
    preset="cinematic",
    rim_light=True
)

# Optimize composition
script = composer.apply_rule_of_thirds(
    script,
    focal_point="binding_site"
)

# Save and render
renderer.save_script(script, "final_cover.pml")
renderer.render("final_cover.pml", output="cover_image.png")
```

**Expected Output Files:**
```
output/
├── cover_render.pml          # PyMOL script
├── cover_render.png          # Final rendered image (3000px)
├── cover_render_settings.txt # Documentation of settings used
└── preview_1200px.png        # Low-res preview (if generated)
```

## Quality Checklist

**Pre-Rendering Checks:**
- [ ] Structure is complete (no missing residues in key regions)
- [ ] Ligands and cofactors properly oriented
- [ ] Hydrogens added if showing surface representation
- [ ] Chain colors distinguishable for multimeric proteins
- [ ] pH-appropriate protonation states

**During Scene Setup:**
- [ ] Camera angle shows key features clearly
- [ ] Lighting highlights important regions (active sites, interfaces)
- [ ] Color scheme is colorblind-friendly (avoid red-green combinations)
- [ ] Background contrasts well with protein colors
- [ ] Scale bar or dimension reference included if needed

**Rendering Settings:**
- [ ] Resolution meets journal requirements (usually 300 DPI × print width)
- [ ] Ray tracing enabled for publication quality
- [ ] **CRITICAL**: Depth of field focused on intended subject
- [ ] Antialiasing enabled for smooth edges
- [ ] Oversampling set appropriately (2× for final renders)

**Post-Rendering Review:**
- [ ] Image is sharp (no blur from camera shake or low resolution)
- [ ] No artifacts (banding, moiré patterns, jagged edges)
- [ ] Colors accurate (check on multiple monitors)
- [ ] Important features not obscured by labels or effects
- [ ] File format appropriate (PNG for web/PowerPoint; TIFF for print)

**Before Submission:**
- [ ] **CRITICAL**: Image reviewed by co-authors for accuracy
- [ ] Figure legend describes what is shown
- [ ] Source data (PDB ID, structure preparation) documented
- [ ] Rendering parameters documented for reproducibility
- [ ] Image compressed appropriately for journal submission guidelines

## Common Pitfalls

**Technical Issues:**
- ❌ **Insufficient resolution** (1200px for journal cover) → Pixelated when printed
  - ✅ Calculate required pixels: (print width in inches) × 300 DPI; minimum 2400px for single-column figures

- ❌ **Wrong file format** (JPEG with compression artifacts) → Loss of detail
  - ✅ Use PNG for digital; TIFF for print; avoid JPEG for line art

- ❌ **Missing ray tracing** → Looks like screenshot, not professional render
  - ✅ Always enable ray tracing for publication; accept longer render times

**Composition Issues:**
- ❌ **Dead-center composition** → Boring, static image
  - ✅ Use rule of thirds; place key features at intersections

- ❌ **Cluttered background** → Distracts from main subject
  - ✅ Use depth of field or plain background; remove unnecessary molecules

- ❌ **Poor color choices** (red protein on red background) → Invisible subject
  - ✅ Ensure contrast between subject and background; check colorblind simulations

**Scientific Accuracy Issues:**
- ❌ **Misleading representations** (surface hides important details) → Readers miss key features
  - ✅ Use insets or multiple views to show detail; label important residues

- ❌ **Incorrect coloring** (rainbow by chain when should be by B-factor) → Wrong information conveyed
  - ✅ Choose coloring scheme that highlights relevant feature (conservation, B-factor, electrostatics)

**Lighting and Effects Issues:**
- ❌ **Overuse of effects** (max fog + max DOF + reflections) → Tacky, obscures structure
  - ✅ Use effects purposefully; subtlety often more effective

- ❌ **Inconsistent lighting** across figure panels → Looks unprofessional
  - ✅ Use same lighting setup for all related figures; document and reuse

**Workflow Issues:**
- ❌ **No preview renders** → Final render has errors, wasted time
  - ✅ Always render at low resolution first; check composition and effects

- ❌ **Not saving scene setup** → Cannot reproduce figure for revision
  - ✅ Save complete script with all settings; version control scene files

## Troubleshooting

**Problem: Rendered image is blurry or pixelated**
- Symptoms: "Edges look jagged" or "details not crisp"
- Causes: Resolution too low; antialiasing disabled; zoomed in too far
- Solutions:
  - Increase resolution to 3000px+ for covers
  - Enable antialiasing (2× or 4× oversampling)
  - Don't zoom beyond 100% of original structure size

**Problem: Depth of field not working**
- Symptoms: "Everything in focus despite DOF settings"
- Causes: Software version doesn't support DOF; focal point not set correctly
- Solutions:
  - Update to PyMOL 2.5+ or ChimeraX 1.5+
  - Verify focal point syntax (residue number, selection)
  - Check if ray tracing enabled (DOF requires ray tracing)

**Problem: Colors look different in final render**
- Symptoms: "Bright colors in preview, muted in final image"
- Causes: Color space conversion; gamma correction
- Solutions:
  - Use consistent color space (sRGB for web, Adobe RGB for print)
  - Check monitor calibration
  - Export with embedded color profile

**Problem: Structure doesn't load**
- Symptoms: "Error: Cannot load PDB file"
- Causes: File path incorrect; PDB format invalid; network issues (if fetching)
- Solutions:
  - Use absolute paths or ensure relative paths are correct
  - Validate PDB file format (should start with HEADER or ATOM)
  - For PDB IDs, check internet connection; download manually if needed

**Problem: Render crashes or hangs**
- Symptoms: "PyMOL freezes" or "out of memory error"
- Causes: Resolution too high; insufficient RAM; large structure
- Solutions:
  - Reduce resolution for preview; render final at high res overnight
  - Close other applications to free RAM
  - Split large structures into domains; render separately
  - Use CPU instead of GPU rendering if GPU memory limited

**Problem: Lighting looks flat or artificial**
- Symptoms: "No sense of depth; looks like clip art"
- Causes: Single light source; no shadows; uniform ambient light
- Solutions:
  - Enable three-point lighting setup (key, fill, rim)
  - Add ambient occlusion for contact shadows
  - Use gradient background instead of solid color
  - Slight rotation (not perfectly head-on) creates depth

**Problem: Can't match journal's style**
- Symptoms: "Editor says figure doesn't match journal aesthetic"
- Causes: Lighting too dramatic or too plain; wrong color scheme
- Solutions:
  - Study recent covers/figures from target journal
  - Use journal's official style guide if available
  - Contact editor for specific requirements
  - Submit multiple options (conservative and dramatic)

## References

Available in `references/` directory:

- `pymol_rendering_guide.md` - PyMOL-specific rendering techniques
- `chimerax_lighting.md` - Advanced lighting in ChimeraX
- `color_schemes.md` - Colorblind-friendly palettes and journal preferences
- `composition_principles.md` - Photography and design for molecular graphics
- `journal_requirements.md` - Submission guidelines for major journals
- `rendering_math.md` - Technical details of ray tracing algorithms

## Scripts

Located in `scripts/` directory:

- `main.py` - CLI interface for script generation
- `renderer.py` - Core rendering script generator
- `effects.py` - Advanced effects configuration (DOF, AO, etc.)
- `lighting.py` - Lighting preset and customization
- `composition.py` - Camera positioning and composition guides
- `batch_processor.py` - Batch rendering multiple structures
- `utils.py` - Helper functions for PDB handling and validation

## Performance and Resources

**Rendering Time Estimates:**
| Resolution | Effects | PyMOL | ChimeraX |
|------------|---------|-------|----------|
| 1200px | Basic ray trace | 30 sec | 20 sec |
| 2400px | Ray trace + AO | 3-5 min | 2-4 min |
| 3000px | Full effects | 10-15 min | 8-12 min |
| 3840px (4K) | Full effects | 20-30 min | 15-25 min |

**Hardware Requirements:**
- **RAM**: 8 GB minimum; 16 GB+ for large structures (>10,000 atoms)
- **CPU**: Multi-core speeds up rendering; ray tracing is CPU-bound
- **GPU**: Optional for preview; most ray tracers use CPU
- **Storage**: 500 MB+ for output images (high-res PNGs are large)

**Optimization Tips:**
- Render final images overnight or during breaks
- Use draft settings (1200px, no ray trace) for layout testing
- Batch render similar scenes to amortize setup time
- Consider render farm or cloud for large batches

## Limitations

- **Static Images Only**: Generates scripts for still images, not animations (use `molecular-dynamics-viz` for animations)
- **Software Dependency**: Requires separately installed PyMOL or ChimeraX
- **Rendering Time**: High-quality renders can take 10-30 minutes per image
- **Learning Curve**: Advanced effects require understanding of photography concepts
- **File Sizes**: High-res images can be 10-50 MB each
- **Limited Editing**: Generated scripts create new scenes; complex edits require manual PyMOL/ChimeraX knowledge
- **No Automatic Layout**: Creates single images; figure assembly requires separate tools (use `multi-panel-figure-assembler`)

## Version History

- **v1.0.0** (Current): Initial release with PyMOL and ChimeraX support, 4 presets, comprehensive effects
- Planned: Blender integration for cinematic animations, AI-assisted composition suggestions, real-time preview mode

---

**💡 Tip: For creating multiple related figures, save your complete scene setup (lighting, camera, colors) as a PyMOL session file (.pse) or ChimeraX session (.cxs), then modify only the specific elements needed for each figure. This ensures consistency across figure panels.**
