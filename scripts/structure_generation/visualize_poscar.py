#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
visualize_poscar.py
===================
Universal code to read POSCAR files and visualize them in VESTA style

Usage:
    python visualize_poscar.py POSCAR [options]
    
Options:
    --title "Custom Title"     : Custom title (default: POSCAR comment line)
    --prefix output            : Output file prefix (default: POSCAR name)
    --dpi 300                  : Resolution (default: 300)
    --show-bonds               : Show bonds (default: False)

Output files:
    {prefix}_sideview.png            : c-axis side view (a-b plane)
    {prefix}_structure_schematic.png : Structure schematic
    {prefix}_poscar_info.png         : POSCAR information
"""

import sys
import argparse
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Circle, FancyBboxPatch

# ============================================================
# Configuration
# ============================================================

# Element colors (VESTA style)
ELEMENT_COLORS = {
    'H': '#FFFFFF',
    'C': '#909090',
    'N': '#3050F8',
    'O': '#FF0D0D',
    'S': '#FFD93D',
    'Se': '#6BCB77',
    'Mo': '#FF6B6B',
    'W': '#9D7CD8',
    'Ti': '#BFC2C7',
    'Zr': '#73E6E6',
    'Hf': '#4DC2FF',
    # More can be added
}

# Element sizes (relative)
ELEMENT_SIZES = {
    'H': 80,
    'C': 100,
    'N': 100,
    'O': 100,
    'S': 150,
    'Se': 180,
    'Mo': 200,
    'W': 200,
    'Ti': 180,
    'Zr': 200,
    'Hf': 200,
}

# ============================================================
# POSCAR Reading
# ============================================================

def read_poscar(filename):
    """Read POSCAR file"""
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    # Comment line
    comment = lines[0].strip()
    
    # Scaling factor
    scale = float(lines[1].strip())
    
    # Lattice vectors
    lattice = []
    for i in range(2, 5):
        vec = [float(x) for x in lines[i].split()]
        lattice.append(vec)
    lattice = np.array(lattice) * scale
    
    # Element names
    elements = lines[5].split()
    
    # Element counts
    counts = [int(x) for x in lines[6].split()]
    
    # Coordinate type
    coord_type = lines[7].strip()[0].upper()
    
    # Positions
    start_line = 8
    total_atoms = sum(counts)
    positions = []
    
    for i in range(start_line, start_line + total_atoms):
        pos = [float(x) for x in lines[i].split()[:3]]
        positions.append(pos)
    
    positions = np.array(positions)
    
    # Convert to Cartesian if needed
    if coord_type == 'D':  # Direct
        positions_cart = positions @ lattice
    else:  # Cartesian
        positions_cart = positions
    
    # Create atom list with element names
    atom_symbols = []
    for elem, count in zip(elements, counts):
        atom_symbols.extend([elem] * count)
    
    return {
        'comment': comment,
        'lattice': lattice,
        'elements': elements,
        'counts': counts,
        'positions': positions,
        'positions_cart': positions_cart,
        'symbols': atom_symbols,
        'coord_type': coord_type
    }


def analyze_layers(structure, axis=2):
    """
    Analyze layer structure along specified axis
    
    Parameters:
    -----------
    structure : dict
        Structure information
    axis : int
        Axis along which to analyze layers (0=x, 1=y, 2=z)
    """
    coords = structure['positions_cart'][:, axis]
    symbols = structure['symbols']
    
    # Find unique levels (tolerance: 0.1 Å)
    unique_levels = []
    tolerance = 0.1
    
    for coord in sorted(set(coords)):
        if not any(abs(coord - ul) < tolerance for ul in unique_levels):
            unique_levels.append(coord)
    
    # Element information for each level
    layers = []
    for level in unique_levels:
        mask = np.abs(coords - level) < tolerance
        layer_symbols = [symbols[i] for i, m in enumerate(mask) if m]
        
        from collections import Counter
        counts = Counter(layer_symbols)
        layers.append({
            'coord': level,
            'elements': dict(counts),
            'total': len(layer_symbols)
        })
    
    return layers


def find_layer_groups(layers, gap_threshold=3.0):
    """Group layers together (heterostructure detection)"""
    if len(layers) < 2:
        return [list(range(len(layers)))]
    
    groups = []
    current_group = [0]
    
    for i in range(1, len(layers)):
        gap = layers[i]['coord'] - layers[i-1]['coord']
        
        if gap > gap_threshold:
            # Start new group
            groups.append(current_group)
            current_group = [i]
        else:
            current_group.append(i)
    
    groups.append(current_group)
    return groups


# ============================================================
# Visualization
# ============================================================

def get_element_color(element):
    """Get element color"""
    return ELEMENT_COLORS.get(element, '#CCCCCC')


def get_element_size(element):
    """Get element size"""
    return ELEMENT_SIZES.get(element, 120)


def plot_side_view_axis(ax, structure, layers, layer_groups, axis=2, 
                        axis_name='c', show_bonds=False):
    """
    Side view along specified axis (VESTA style)
    
    Parameters:
    -----------
    axis : int
        Axis to view along (0=a/x, 1=b/y, 2=c/z)
    axis_name : str
        Name of the axis for labeling
    """
    ax.set_title(f'Side View along {axis_name}-axis', fontsize=14, fontweight='bold', pad=15)
    
    # Set labels based on axis
    if axis == 0:  # a-axis view
        ax.set_xlabel('b-c plane', fontsize=11)
        ax.set_ylabel('a-axis (Angstrom)', fontsize=11)
        height = structure['lattice'][0, 0]
        plane_name = 'b-c'
    elif axis == 1:  # b-axis view
        ax.set_xlabel('a-c plane', fontsize=11)
        ax.set_ylabel('b-axis (Angstrom)', fontsize=11)
        height = structure['lattice'][1, 1]
        plane_name = 'a-c'
    else:  # c-axis view (default)
        ax.set_xlabel('a-b plane', fontsize=11)
        ax.set_ylabel('c-axis (Angstrom)', fontsize=11)
        height = structure['lattice'][2, 2]
        plane_name = 'a-b'
    
    ax.set_xlim(-2, 8)
    ax.set_ylim(-0.5, height + 0.5)
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.2, axis='y')
    
    # Axis arrow (moved outside plot area)
    arrow_x = -1.5
    arrow_y_start = max(-0.3, 0)
    arrow_length = min(2.5, height * 0.15)
    ax.arrow(arrow_x, arrow_y_start, 0, arrow_length, 
            head_width=0.15, head_length=0.2,
            fc='blue', ec='blue', linewidth=2, clip_on=False)
    ax.text(arrow_x - 0.15, arrow_y_start + arrow_length + 0.3, axis_name, 
           fontsize=14, fontweight='bold', color='blue', clip_on=False)
    
    # Draw each layer
    positions_cart = structure['positions_cart']
    symbols = structure['symbols']
    
    # Determine x positions (arrange atoms in a row)
    x_positions = np.linspace(1, 6, 5)
    
    # Track elements at each y-level to avoid overlap
    from collections import defaultdict
    elements_at_level = defaultdict(list)
    
    for layer_idx, layer in enumerate(layers):
        coord = layer['coord']
        
        # Skip if outside visible range
        if coord < -0.5 or coord > height + 0.5:
            continue
        
        # Atoms at this level
        mask = np.abs(positions_cart[:, axis] - coord) < 0.1
        layer_atoms = [(symbols[i], positions_cart[i]) for i, m in enumerate(mask) if m]
        
        # Collect elements at this level
        for elem in layer['elements']:
            elements_at_level[coord].append(elem)
        
        # Draw by element
        for elem in layer['elements']:
            elem_count = layer['elements'][elem]
            n_draw = min(elem_count, len(x_positions))
            
            for i in range(n_draw):
                x = x_positions[i]
                circle = Circle((x, coord), 0.35, 
                              facecolor=get_element_color(elem),
                              edgecolor='black', linewidth=2, zorder=5)
                ax.add_patch(circle)
    
    # Add element labels - fixed position by element type to avoid overlap
    sorted_coords = sorted(elements_at_level.keys())
    
    for coord in sorted_coords:
        elems = elements_at_level[coord]
        
        for elem in sorted(set(elems)):
            # Only show if within visible range
            if not (-0.5 <= coord <= height + 0.5):
                continue
            
            # Position based on element type
            if elem in ['Mo', 'W', 'Ti', 'Zr', 'Hf']:  # Transition metals - left side
                label_x = 0.15
                ha = 'right'
            else:  # S, Se, etc. (Chalcogens) - right side
                label_x = 0.55
                ha = 'left'
            
            ax.text(label_x, coord, elem, fontsize=10, fontweight='bold',
                   ha=ha, va='center', clip_on=True, zorder=10,
                   bbox=dict(boxstyle='round,pad=0.15', facecolor='white', 
                            alpha=0.9, edgecolor='none'))
    
    # Draw layer group boxes
    colors_box = ['red', 'blue', 'green', 'orange']
    labels_box = ['Layer 1', 'Layer 2', 'Layer 3', 'Layer 4']
    
    for group_idx, group in enumerate(layer_groups):
        if len(group) == 0:
            continue
        
        coord_min = layers[group[0]]['coord']
        coord_max = layers[group[-1]]['coord']
        
        # Only draw if within visible range
        if coord_max < -0.5 or coord_min > height + 0.5:
            continue
        
        # Clip to visible range
        coord_min = max(coord_min, -0.5)
        coord_max = min(coord_max, height + 0.5)
        box_height = coord_max - coord_min + 0.5
        
        color = colors_box[group_idx % len(colors_box)]
        label = labels_box[group_idx % len(labels_box)]
        
        ax.add_patch(FancyBboxPatch((0.8, coord_min - 0.2), 5.7, box_height,
                                    boxstyle="round,pad=0.1",
                                    edgecolor=color, facecolor='none',
                                    linewidth=2.5, linestyle='-', alpha=0.6))
        
        ax.text(7.2, (coord_min + coord_max) / 2, label, fontsize=11, fontweight='bold',
               bbox=dict(boxstyle='round', facecolor='white', 
                        edgecolor=color, linewidth=2), clip_on=False)
    
    # Show vdW gap
    if len(layer_groups) > 1:
        for i in range(len(layer_groups) - 1):
            coord1 = layers[layer_groups[i][-1]]['coord']
            coord2 = layers[layer_groups[i+1][0]]['coord']
            
            # Only show if within visible range
            if coord1 < -0.5 or coord2 > height + 0.5:
                continue
            
            gap = coord2 - coord1
            
            ax.annotate('', xy=(-0.3, coord1), xytext=(-0.3, coord2),
                       arrowprops=dict(arrowstyle='<->', color='green', lw=3),
                       clip_on=False)
            ax.text(0.0, (coord1 + coord2) / 2, f'{gap:.2f} Å\nvdW gap',
                   fontsize=9, color='green', fontweight='bold', va='center',
                   clip_on=False)


def plot_schematic(ax, structure, layers, layer_groups):
    """Schematic view - layers aligned by z-coordinate"""
    ax.set_title('Structure Schematic (to scale)', fontsize=14, fontweight='bold')
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 16)
    ax.axis('off')
    
    # Calculate total z range
    z_min = layers[0]['coord']
    z_max = layers[-1]['coord']
    z_range = z_max - z_min
    
    # Scale y coordinates (fit between 2~14)
    y_scale = 10 / z_range if z_range > 0 else 1
    y_offset = 2
    
    # Draw each layer group at different x positions
    x_positions = [2.5, 7.0]
    colors_group = ['red', 'blue', 'green', 'orange']
    
    for group_idx, group in enumerate(layer_groups[:2]):
        x_base = x_positions[group_idx]
        color_group = colors_group[group_idx % len(colors_group)]
        
        # Layers in this group
        for layer_idx in group:
            layer = layers[layer_idx]
            z_real = layer['coord']
            
            # Convert real z coordinate to y coordinate
            y_pos = y_offset + (z_real - z_min) * y_scale
            
            # Display each element as a layer
            for elem in sorted(layer['elements'].keys()):
                # Determine height (fixed for visibility)
                if elem in ['Mo', 'W', 'Ti', 'Zr', 'Hf']:
                    height = 0.3
                    label_color = 'white'
                else:
                    height = 0.25
                    label_color = 'black'
                
                rect = mpatches.Rectangle((x_base - 1.5, y_pos), 3, height,
                                         facecolor=get_element_color(elem),
                                         edgecolor='black', linewidth=2)
                ax.add_patch(rect)
                
                ax.text(x_base, y_pos + height/2, elem,
                       ha='center', va='center', fontsize=10,
                       fontweight='bold', color=label_color)
                
                # Show z coordinate
                ax.text(x_base + 1.8, y_pos + height/2, f'z={z_real:.2f}',
                       ha='left', va='center', fontsize=7, style='italic')
        
        # Group label
        group_y_min = y_offset + (layers[group[0]]['coord'] - z_min) * y_scale
        group_y_max = y_offset + (layers[group[-1]]['coord'] - z_min) * y_scale
        
        ax.text(x_base - 2.2, (group_y_min + group_y_max) / 2,
               f'Layer {group_idx + 1}',
               rotation=90, ha='right', va='center', fontsize=10,
               fontweight='bold',
               bbox=dict(boxstyle='round', facecolor='white',
                        edgecolor=color_group, linewidth=2))
    
    # Show vdW gap (between layer groups)
    if len(layer_groups) > 1:
        for i in range(len(layer_groups) - 1):
            z1_real = layers[layer_groups[i][-1]]['coord']
            z2_real = layers[layer_groups[i+1][0]]['coord']
            gap_real = z2_real - z1_real
            
            y1 = y_offset + (z1_real - z_min) * y_scale
            y2 = y_offset + (z2_real - z_min) * y_scale
            
            ax.annotate('', xy=(5, y1), xytext=(5, y2),
                       arrowprops=dict(arrowstyle='<->', color='green', lw=2))
            ax.text(5.3, (y1 + y2) / 2, f'vdW gap\n{gap_real:.2f} Å',
                   ha='left', va='center', fontsize=9, color='green',
                   fontweight='bold')


def plot_info(ax, structure):
    """POSCAR/POTCAR information"""
    ax.set_title('POSCAR Information', fontsize=14, fontweight='bold')
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 12)
    ax.axis('off')
    
    # POSCAR information
    y = 10
    ax.text(5, y, 'POSCAR Format', ha='center', fontsize=13, fontweight='bold',
           bbox=dict(boxstyle='round', facecolor='lightblue', linewidth=2))
    
    y -= 1.2
    total_atoms = len(structure['symbols'])
    lattice = structure['lattice']
    a_axis = np.linalg.norm(lattice[0])
    b_axis = np.linalg.norm(lattice[1])
    c_axis = np.linalg.norm(lattice[2])
    
    info_text = f"""File: {structure['comment'][:30]}...
Total atoms: {total_atoms}
Cell a-axis: {a_axis:.4f} Å
Cell b-axis: {b_axis:.4f} Å
Cell c-axis: {c_axis:.4f} Å"""
    
    ax.text(5, y, info_text, ha='center', va='top', fontsize=9, family='monospace',
           bbox=dict(boxstyle='round', facecolor='white', edgecolor='gray', linewidth=1.5))
    
    # Element order
    y -= 4
    ax.text(5, y, 'Element Order', ha='center', fontsize=12, fontweight='bold',
           bbox=dict(boxstyle='round', facecolor='#FFE6CC', linewidth=2))
    
    y -= 0.8
    for i, (elem, count) in enumerate(zip(structure['elements'], structure['counts']), 1):
        rect = mpatches.Rectangle((1.5, y - 0.35), 0.4, 0.3,
                                  facecolor=get_element_color(elem),
                                  edgecolor='black', linewidth=2)
        ax.add_patch(rect)
        ax.text(2.2, y - 0.2, f'{i}. {elem:3s} ({count:2d} atoms)',
               va='center', fontsize=11, fontweight='bold', family='monospace')
        y -= 0.6
    
    # POTCAR information
    y -= 0.5
    ax.text(5, y, 'POTCAR Order:', ha='center', fontsize=11, fontweight='bold')
    y -= 0.5
    
    potcar_lines = []
    for elem in structure['elements']:
        if elem in ['Mo', 'W', 'Ti', 'Zr', 'Hf']:
            potcar_lines.append(f"{elem}_pv")
        else:
            potcar_lines.append(f"{elem}")
    
    potcar_text = '\n'.join(potcar_lines)
    ax.text(5, y, potcar_text, ha='center', va='top', fontsize=10, family='monospace',
           bbox=dict(boxstyle='round', facecolor='#E6FFE6', edgecolor='green', linewidth=2))


def visualize_poscar(poscar_file, title=None, output_prefix=None,
                     dpi=300, show_bonds=False):
    """
    Visualize POSCAR file with multiple views and save separate files
    
    Parameters:
    -----------
    poscar_file : str
        Path to POSCAR file
    title : str, optional
        Custom title (uses POSCAR comment line if None)
    output_prefix : str, optional
        Output file prefix (uses POSCAR filename if None)
    dpi : int
        Resolution
    show_bonds : bool
        Whether to show bonds
    
    Returns:
    --------
    structure : dict
        Structure information
    """
    # Read POSCAR
    structure = read_poscar(poscar_file)
    
    # Determine output prefix
    if output_prefix is None:
        base_name = os.path.splitext(os.path.basename(poscar_file))[0]
        output_prefix = base_name
    
    # Set title
    if title is None:
        title = structure['comment']
    
    # ============================================================
    # 1. Side View (c-axis only)
    # ============================================================
    fig_side, ax = plt.subplots(1, 1, figsize=(10, 8))
    fig_side.suptitle(f'{title} - Side View', fontsize=16, fontweight='bold', y=0.98)
    
    # Analyze layers for c-axis (z-direction)
    layers_c = analyze_layers(structure, axis=2)
    layer_groups_c = find_layer_groups(layers_c)
    plot_side_view_axis(ax, structure, layers_c, layer_groups_c, 
                       axis=2, axis_name='c', show_bonds=show_bonds)
    
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    sideview_file = f'{output_prefix}_sideview.png'
    plt.savefig(sideview_file, dpi=dpi, bbox_inches='tight')
    plt.close()
    print(f"✓ Side view saved: {sideview_file}")
    
    # ============================================================
    # 2. Structure Schematic
    # ============================================================
    fig_schem, ax = plt.subplots(1, 1, figsize=(10, 12))
    fig_schem.suptitle(f'{title} - Structure Schematic', 
                       fontsize=16, fontweight='bold', y=0.98)
    
    # Use c-axis for schematic (z-direction)
    layers_z = analyze_layers(structure, axis=2)
    layer_groups_z = find_layer_groups(layers_z)
    plot_schematic(ax, structure, layers_z, layer_groups_z)
    
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    schematic_file = f'{output_prefix}_structure_schematic.png'
    plt.savefig(schematic_file, dpi=dpi, bbox_inches='tight')
    plt.close()
    print(f"✓ Structure schematic saved: {schematic_file}")
    
    # ============================================================
    # 3. POSCAR Information
    # ============================================================
    fig_info, ax = plt.subplots(1, 1, figsize=(10, 12))
    fig_info.suptitle(f'{title} - POSCAR Information', 
                      fontsize=16, fontweight='bold', y=0.98)
    
    plot_info(ax, structure)
    
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    info_file = f'{output_prefix}_poscar_info.png'
    plt.savefig(info_file, dpi=dpi, bbox_inches='tight')
    plt.close()
    print(f"✓ POSCAR info saved: {info_file}")
    
    print(f"\n✓ All visualizations complete!")
    print(f"  Generated files:")
    print(f"    - {sideview_file}")
    print(f"    - {schematic_file}")
    print(f"    - {info_file}")
    
    return structure


# ============================================================
# Main
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description='Visualize POSCAR file in VESTA style',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python visualize_poscar.py POSCAR
  python visualize_poscar.py POSCAR --title "My Structure"
  python visualize_poscar.py POSCAR --prefix my_structure --dpi 600
  
Output files:
  {prefix}_sideview.png            : c-axis side view (a-b plane)
  {prefix}_structure_schematic.png : Structure schematic
  {prefix}_poscar_info.png         : POSCAR information
        """
    )
    
    parser.add_argument('poscar', help='Path to POSCAR file')
    parser.add_argument('--title', type=str, default=None,
                       help='Custom title (default: POSCAR comment line)')
    parser.add_argument('--prefix', type=str, default=None,
                       help='Output file prefix (default: POSCAR filename)')
    parser.add_argument('--dpi', type=int, default=300,
                       help='Resolution (default: 300)')
    parser.add_argument('--show-bonds', action='store_true',
                       help='Show bonds')
    
    args = parser.parse_args()
    
    # Execute visualization
    try:
        visualize_poscar(
            args.poscar,
            title=args.title,
            output_prefix=args.prefix,
            dpi=args.dpi,
            show_bonds=args.show_bonds
        )
    except Exception as e:
        print(f"✗ Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
