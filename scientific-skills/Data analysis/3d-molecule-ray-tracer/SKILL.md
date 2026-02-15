---
name: 3d-molecule-ray-tracer
description: "Generate photorealistic rendering scripts for PyMOL and UCSF ChimeraX 
  to create publication-quality molecular visualizations. Supports ray-tracing, 
  depth of field, ambient occlusion, and cinematic lighting for journal covers 
  and high-impact figures. Use for creating visually compelling protein structure 
  images; not for routine structural analysis or measurements."
version: 1.0.0
category: Visual
tags: ["pymol", "chimerax", "rendering", "visualization", "ray-tracing"]
author: AIPOCH
license: MIT
status: Draft
risk_level: High
skill_type: Hybrid (Tool/Script + Network/API)
owner: AIPOCH
reviewer: ''
last_updated: '2026-02-15'
---

# 3D Molecule Ray Tracer

Advanced molecular visualization tool that generates professional-grade rendering scripts with cinematic effects for creating publication-quality and cover-worthy molecular images.

## Features

- **Multi-Software Support**: Generate scripts for PyMOL and UCSF ChimeraX
- **Photorealistic Rendering**: Ray-tracing, depth of field, ambient occlusion
- **Cinematic Lighting**: Studio, outdoor, and dramatic lighting presets
- **Publication Presets**: Pre-configured settings for journals, covers, and presentations
- **Customizable Scenes**: Fine control over camera, materials, and atmosphere

## Usage

### Basic Usage

```bash
# Generate PyMOL script with default settings
python scripts/main.py --pdb 1mbn

# Generate cover-quality render script
python scripts/main.py --pdb 1mbn --preset cover

# Generate ChimeraX script
python scripts/main.py --software chimerax --pdb 1abc --preset publication
```

### Parameters

| Parameter | Type | Default | Required | Description |
|-----------|------|---------|----------|-------------|
| `--software` | str | pymol | No | Target rendering software (pymol/chimerax) |
| `--pdb` | str | None | Yes | PDB file path or 4-letter PDB ID |
| `--preset` | str | standard | No | Rendering preset (standard/cover/publication/cinematic) |
| `--style` | str | cartoon | No | Molecular representation style |
| `--resolution` | int | from preset | No | Output resolution in pixels |
| `--bg-color` | str | white | No | Background color |
| `--ao-on` | flag | False | No | Enable ambient occlusion |
| `--shadows` | flag | False | No | Enable shadow casting |
| `--fog` | float | from preset | No | Fog density (0-1) |
| `--dof-on` | flag | False | No | Enable depth of field |
| `--dof-focus` | str | center | No | DOF focus point |
| `--dof-aperture` | float | from preset | No | Aperture size (higher = more blur) |
| `--lighting` | str | from preset | No | Lighting preset |
| `--output` | str | auto | No | Output script filename |

### Advanced Usage

```bash
# Cover-quality render with depth of field
python scripts/main.py \
  --software pymol \
  --pdb 1mbn \
  --preset cover \
  --dof-on \
  --dof-focus "A:64" \
  --dof-aperture 2.0 \
  --style surface \
  --output cover_render.pml

# Cinematic 4K render
python scripts/main.py \
  --software pymol \
  --pdb complex.pdb \
  --preset cinematic \
  --resolution 3840 \
  --ao-on \
  --shadows \
  --lighting cinematic
```

## Rendering Presets

| Preset | Resolution | Ray Trace | DOF | AO | Shadows | Use Case |
|--------|------------|-----------|-----|-----|---------|----------|
| **Standard** | 2400px | ✓ | ✗ | ✗ | ✗ | Quick high-quality |
| **Cover** | 3000px | ✓ | ✓ | ✓ | ✓ | Journal covers |
| **Publication** | 2400px | ✓ | ✗ | ✓ | ✗ | Manuscript figures |
| **Cinematic** | 3840px | ✓ | ✓ | ✓ | ✓ | Presentations |

## Supported Software

| Software | Best For | Features |
|----------|----------|----------|
| **PyMOL** | Traditional rendering, ease of use | Ray tracing, shadows, AO |
| **ChimeraX** | Modern effects, large structures | PBR lighting, ambient occlusion, VR |

## Technical Difficulty: **MEDIUM**

⚠️ **AI自主验收状态**: 需人工检查

This skill requires:
- Python 3.8+ environment
- PyMOL 2.5+ or ChimeraX 1.5+ installed separately
- Understanding of molecular visualization concepts

## Dependencies

### Required Python Packages

```bash
pip install -r requirements.txt
```

### External Software

- **PyMOL**: https://pymol.org/
- **UCSF ChimeraX**: https://www.cgl.ucsf.edu/chimerax/

## Risk Assessment

| Risk Indicator | Assessment | Level |
|----------------|------------|-------|
| Code Execution | Python scripts executed locally | Medium |
| Network Access | Fetches PDB structures from RCSB (optional) | Low |
| File System Access | Writes rendering scripts | Low |
| Instruction Tampering | Standard prompt guidelines | Low |
| Data Exposure | No sensitive data exposure | Low |

## Security Checklist

- [x] No hardcoded credentials or API keys
- [x] No unauthorized file system access (../)
- [x] Output does not expose sensitive information
- [x] Prompt injection protections in place
- [x] Input file paths validated
- [x] Output directory restricted to workspace
- [x] Script execution in sandboxed environment
- [x] Error messages sanitized
- [x] Dependencies audited

## Prerequisites

```bash
# Python dependencies
pip install -r requirements.txt

# Install PyMOL or ChimeraX separately
```

## Output Example

```
✓ Rendering script generated: /path/to/cover_render.pml

Configuration:
  Software: pymol
  Preset: cover
  Style: cartoon
  Resolution: 3000px
  Depth of Field: ON
  Ambient Occlusion: ON
  Shadows: ON
  Lighting: cinematic

To render:
  pymol cover_render.pml
  # Or within PyMOL:
  @ cover_render.pml
```

## Evaluation Criteria

### Success Metrics
- [ ] Successfully generates valid PyMOL/ChimeraX scripts
- [ ] Scripts execute without errors in target software
- [ ] Output images meet quality standards
- [ ] Handles edge cases gracefully

### Test Cases
1. **Basic Functionality**: Generate script for PDB ID → Valid script created
2. **File Input**: Generate script from PDB file → Valid script created
3. **Preset Override**: Custom parameters override preset → Correct settings applied
4. **Both Software**: Generate for PyMOL and ChimeraX → Both scripts valid

## Lifecycle Status

- **Current Stage**: Draft
- **Next Review Date**: 2026-03-15
- **Known Issues**: None
- **Planned Improvements**: 
  - Blender integration
  - AI-assisted composition suggestions
  - Real-time preview mode

## References

See `references/` for:
- PyMOL-specific rendering techniques
- ChimeraX lighting documentation
- Colorblind-friendly palettes
- Journal submission guidelines

## Limitations

- **Static Images Only**: Generates scripts for still images, not animations
- **Software Dependency**: Requires separately installed PyMOL or ChimeraX
- **Rendering Time**: High-quality renders can take 10-30 minutes per image
- **Learning Curve**: Advanced effects require understanding of photography concepts
- **File Sizes**: High-res images can be 10-50 MB each
- **No Automatic Layout**: Creates single images; figure assembly requires separate tools

---

**💡 Tip: For creating multiple related figures, save your complete scene setup (lighting, camera, colors) as a PyMOL session file (.pse) or ChimeraX session (.cxs), then modify only the specific elements needed for each figure. This ensures consistency across figure panels.**
