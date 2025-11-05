#!/bin/bash
# clean_vasp_outputs.sh
# Script to clean VASP output files safely

echo "========================================"
echo "VASP Output Files Cleanup Script"
echo "========================================"
echo ""

# Target directories
DIRS="MoSSe_Se_up__WSSe_S_up MoSSe_S_up__WSSe_Se_up MoSSe_Se_up__WSSe_Se_up MoSSe_S_up__WSSe_S_up"

# Check if any directory exists
DIR_EXISTS=0
for dir in $DIRS; do
    if [ -d "$dir" ]; then
        DIR_EXISTS=1
        break
    fi
done

if [ $DIR_EXISTS -eq 0 ]; then
    echo "ERROR: No target directories found!"
    echo "Expected directories:"
    for dir in $DIRS; do
        echo "  - $dir"
    done
    exit 1
fi

# Backup option
echo "This script will delete VASP output files from:"
for dir in $DIRS; do
    if [ -d "$dir" ]; then
        echo "  - $dir"
    fi
done
echo ""
echo "Files to be deleted:"
echo "  - OUTCAR, OSZICAR, CONTCAR, XDATCAR"
echo "  - WAVECAR, CHGCAR, CHG"
echo "  - vasprun.xml, EIGENVAL, DOSCAR, PROCAR"
echo "  - *.out, *.log, IBZKPT, PCDAT"
echo ""
echo "Input files will be PRESERVED:"
echo "  - POSCAR, INCAR, KPOINTS, POTCAR"
echo ""

read -p "Create backup before deletion? (y/n): " backup_choice
echo ""

if [ "$backup_choice" = "y" ] || [ "$backup_choice" = "Y" ]; then
    BACKUP_DIR="old_calculations_backup_$(date +%Y%m%d_%H%M%S)"
    mkdir -p "$BACKUP_DIR"
    echo "Backup directory created: $BACKUP_DIR"
    echo ""
fi

# Process each directory
for dir in $DIRS; do
    if [ -d "$dir" ]; then
        echo "========================================"
        echo "Processing: $dir"
        echo "========================================"
        cd "$dir"
        
        # Check if output files exist
        OUTPUT_EXISTS=0
        for file in OUTCAR OSZICAR CONTCAR WAVECAR; do
            if [ -f "$file" ]; then
                OUTPUT_EXISTS=1
                break
            fi
        done
        
        if [ $OUTPUT_EXISTS -eq 0 ]; then
            echo "  No output files found - skipping"
            cd ..
            continue
        fi
        
        # Create backup if requested
        if [ "$backup_choice" = "y" ] || [ "$backup_choice" = "Y" ]; then
            echo "  Creating backup..."
            tar -czf "../${BACKUP_DIR}/${dir}.tar.gz" \
                OUTCAR OSZICAR CONTCAR XDATCAR \
                WAVECAR CHGCAR CHG \
                vasprun.xml EIGENVAL DOSCAR PROCAR \
                *.out *.log IBZKPT PCDAT 2>/dev/null
            
            if [ $? -eq 0 ]; then
                echo "  ✓ Backup created: ${BACKUP_DIR}/${dir}.tar.gz"
            else
                echo "  ✗ Backup failed (some files may not exist)"
            fi
        fi
        
        # Delete output files
        echo "  Deleting output files..."
        rm -f OUTCAR OSZICAR CONTCAR XDATCAR
        rm -f WAVECAR CHGCAR CHG CHGCAR.*
        rm -f vasprun.xml EIGENVAL DOSCAR PROCAR
        rm -f *.out *.log PCDAT IBZKPT
        echo "  ✓ Output files deleted"
        
        # Verify input files still exist
        echo "  Checking input files..."
        ALL_INPUTS_OK=1
        for file in POSCAR INCAR KPOINTS POTCAR; do
            if [ -f "$file" ]; then
                echo "    ✓ $file"
            else
                echo "    ✗ $file MISSING!"
                ALL_INPUTS_OK=0
            fi
        done
        
        if [ $ALL_INPUTS_OK -eq 1 ]; then
            echo "  ✓ All input files intact"
        else
            echo "  ✗ WARNING: Some input files are missing!"
        fi
        
        # Show disk usage
        echo "  Current directory size: $(du -sh . | cut -f1)"
        echo ""
        
        cd ..
    else
        echo "Directory not found: $dir - skipping"
    fi
done

echo "========================================"
echo "Cleanup Complete!"
echo "========================================"
echo ""

if [ "$backup_choice" = "y" ] || [ "$backup_choice" = "Y" ]; then
    echo "Backups saved in: $BACKUP_DIR"
    echo "Backup total size: $(du -sh $BACKUP_DIR | cut -f1)"
    echo ""
fi

echo "Summary of cleaned directories:"
for dir in $DIRS; do
    if [ -d "$dir" ]; then
        echo "  $dir: $(du -sh $dir | cut -f1)"
    fi
done
echo ""
echo "You can now copy new PBS and INCAR files to each directory."
