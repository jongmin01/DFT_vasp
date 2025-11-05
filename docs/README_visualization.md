# POSCAR Visualization Tool v2 - User Guide

## 📖 Overview

`visualize_poscar.py` is a Python script that reads VASP POSCAR files and automatically generates high-quality VESTA-style visualizations.

## 🆕 v2 Major Changes

- **3-Axis Side Views**: Creates views along a, b, and c axes
- **3 Separate Files**: Side views, Structure schematic, and Information saved as individual PNG files
- **File Naming Convention**: 
  - `{POSCAR_name}_sideview.png`
  - `{POSCAR_name}_structure_schematic.png`
  - `{POSCAR_name}_information.png`

---

## 🚀 Basic Usage

### 1. Simplest Use Case
```bash
python visualize_poscar.py POSCAR
```
- Creates 3 PNG files:
  - `POSCAR_sideview.png` (a, b, c axis views)
  - `POSCAR_structure_schematic.png`
  - `POSCAR_information.png`

### 2. Custom Title
```bash
python visualize_poscar.py POSCAR --title "MoSSe/WSSe Heterostructure"
```

### 3. Specify Output Filename Prefix
```bash
python visualize_poscar.py POSCAR --output-prefix my_structure
```
- Generated files:
  - `my_structure_sideview.png`
  - `my_structure_structure_schematic.png`
  - `my_structure_information.png`

### 4. High Resolution Output
```bash
python visualize_poscar.py POSCAR --dpi 600 --output-prefix high_res
```

### 5. Complete Example with All Options
```bash
python visualize_poscar.py MoSSe_WSSe/POSCAR \
    --title "Janus TMDC Heterostructure" \
    --output-prefix MoSSe_WSSe \
    --dpi 300
```

---

## 📊 Visualization Content

### File 1: `{prefix}_sideview.png` (3 Panels)

**a-axis View (Left)**
- VESTA-style side view
- a-axis pointing upward
- b-c plane projection

**b-axis View (Center)**
- VESTA-style side view
- b-axis pointing upward
- a-c plane projection

**c-axis View (Right)**
- VESTA-style side view
- c-axis pointing upward (default view)
- a-b plane projection

All views include:
- Atoms displayed as spheres
- Color-coded boxes for each layer
- vdW gap indication (for heterostructures)

### File 2: `{prefix}_structure_schematic.png`

- Layers **aligned to actual z-coordinates**
- Each element shown as a box
- Accurately reflects interlayer distances
- Mo and W vertically aligned
- Explicit vdW gap labels

### File 3: `{prefix}_information.png`

- File information (total atoms, cell size, etc.)
- All a, b, c axis lengths displayed
- Element order and counts
- **Auto-generated POTCAR order**
  - Mo, W, Ti, Zr, Hf → `_pv` suffix added
  - S, Se, etc. → as is

---

## 🎨 Element Colors and Sizes

### Supported Elements:
| Element | Color | Use Case |
|---------|-------|----------|
| Mo | Red (#FF6B6B) | TMDC |
| W | Purple (#9D7CD8) | TMDC |
| S | Yellow (#FFD93D) | Chalcogen |
| Se | Green (#6BCB77) | Chalcogen |
| Ti | Gray (#BFC2C7) | TMDC |
| Zr | Cyan (#73E6E6) | TMDC |
| Hf | Blue (#4DC2FF) | TMDC |

### Adding New Elements:
Add to `ELEMENT_COLORS` and `ELEMENT_SIZES` dictionaries at the top of the code:

```python
ELEMENT_COLORS = {
    # ... existing elements ...
    'V': '#A6CEE3',  # Example: Adding Vanadium
}

ELEMENT_SIZES = {
    # ... existing elements ...
    'V': 190,
}
```

---

## 💡 Advanced Usage

### Direct Use in Python:
```python
from visualize_poscar import visualize_poscar

structure, layers, layer_groups = visualize_poscar(
    'POSCAR',
    title='My Custom Title',
    output_prefix='my_output',
    dpi=600
)

# structure: structure information dictionary
# layers: c-axis layer analysis results
# layer_groups: heterostructure layer groups
```

### Returned Data Structures:

**structure:**
```python
{
    'comment': 'first line of file',
    'lattice': np.array([[...], [...], [...]]),  # 3x3 lattice vectors
    'elements': ['Mo', 'W', 'S', 'Se'],
    'counts': [25, 25, 50, 50],
    'positions': np.array(...),  # fractional coordinates
    'positions_cart': np.array(...),  # Cartesian coordinates
    'symbols': ['Mo', 'Mo', ..., 'Se'],  # element symbol for each atom
}
```

**layers:**
```python
[
    {
        'coord': 0.0,  # Cartesian coordinate along axis (Å)
        'elements': {'Se': 25},  # elements at this level
        'total': 25
    },
    {
        'coord': 1.56,
        'elements': {'Mo': 25},
        'total': 25
    },
    # ...
]
```

---

## 🔧 Troubleshooting

### 1. "ModuleNotFoundError: No module named 'matplotlib'"
```bash
pip install matplotlib numpy
```

### 2. Certain elements appear in gray (#CCCCCC)
→ Need to add the element to `ELEMENT_COLORS`

### 3. Axis view looks strange
- Layer analysis is performed independently for each axis
- Adjust `gap_threshold` in `find_layer_groups` function if needed (default: 3.0 Å)

### 4. Only 1 file created instead of 3
- You might be using the old version
- Check v2: Code header should mention "a, b, c axis views"

### 5. Visualization is too small or large
- Adjust `--dpi` value (default: 300)
- Change `figsize` in code:
  - Side views: (20, 8)
  - Schematic/Info: (10, 10)

---

## 📝 Usage Examples

### Example 1: Single TMDC Layer
```bash
python visualize_poscar.py MoS2.vasp \
    --title "MoS2 Monolayer" \
    --output-prefix MoS2
```
Generated files:
- `MoS2_sideview.png`
- `MoS2_structure_schematic.png`
- `MoS2_information.png`

### Example 2: Heterostructure
```bash
python visualize_poscar.py MoSSe_WSSe_hetero.vasp \
    --title "MoSSe/WSSe Janus Heterostructure" \
    --output-prefix MoSSe_WSSe \
    --dpi 600
```

### Example 3: Batch Processing Multiple Files
```bash
# Bash script
for poscar in */POSCAR; do
    dir=$(dirname "$poscar")
    python visualize_poscar.py "$poscar" \
        --title "$dir" \
        --output-prefix "$dir"
done
```

Each directory will have 3 PNG files:
- `dirname_sideview.png`
- `dirname_structure_schematic.png`
- `dirname_information.png`

---

## 🎯 Key Features

### ✅ Automatic Recognition:
- Automatic Direct/Cartesian coordinate conversion
- Independent layer structure analysis for each axis
- Automatic heterostructure vdW gap detection
- Automatic POTCAR order generation

### ✅ VESTA Style:
- Side view with each axis pointing upward
- Atoms displayed as spheres
- Reflects actual crystal structure

### ✅ 3-Axis Views:
- **a-axis view**: Structure seen from b-c plane
- **b-axis view**: Structure seen from a-c plane
- **c-axis view**: Structure seen from a-b plane (default)

### ✅ Separate Files:
- Side views in one file
- Structure schematic in separate file
- Information in separate file
- Each can be used independently

---

## 🆚 Differences from v1

| Feature | v1 | v2 |
|---------|----|----|
| Side views | c-axis only | a, b, c axes |
| Output files | 1 (combined) | 3 (separate) |
| Filename option | `--output` | `--output-prefix` |
| Cell info | c-axis only | a, b, c all |

---

## 📚 Related Files

- `visualize_poscar.py` - Main script (v2)
- `README_visualization_v2.md` - This document

---

## 🆘 Help

```bash
python visualize_poscar.py --help
```

Output:
```
usage: visualize_poscar.py [-h] [--title TITLE] 
                           [--output-prefix OUTPUT_PREFIX] 
                           [--dpi DPI]
                           poscar

Visualize POSCAR file in VESTA style (a, b, c axis views + 3 separate files)

positional arguments:
  poscar                Path to POSCAR file

optional arguments:
  -h, --help            show this help message and exit
  --title TITLE         Custom title (default: POSCAR comment line)
  --output-prefix OUTPUT_PREFIX
                        Output filename prefix (default: POSCAR filename)
  --dpi DPI             Resolution (default: 300)

Examples:
  python visualize_poscar.py POSCAR
  python visualize_poscar.py POSCAR --title "My Structure"
  python visualize_poscar.py POSCAR --output-prefix my_structure --dpi 600

Output files:
  {prefix}_sideview.png            - a, b, c axis side views
  {prefix}_structure_schematic.png - layer structure schematic
  {prefix}_information.png         - POSCAR/POTCAR information
```

---

## 📮 Questions and Bug Reports

If you have issues or need new features:
1. Modify the code (comments are detailed)
2. Customize `ELEMENT_COLORS`/`ELEMENT_SIZES` dictionaries
3. Directly modify `plot_*` functions
4. Use `axis` parameter for analyzing different directions

---

**Made with ❤️ for VASP users**

**Version 2.0** - Multiple axis views & Separate file outputs
