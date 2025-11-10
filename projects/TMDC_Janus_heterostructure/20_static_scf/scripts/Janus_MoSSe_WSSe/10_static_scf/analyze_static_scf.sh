#!/bin/bash
# analyze_static_scf.sh
# Quick analysis of static SCF results

echo "========================================"
echo "Static SCF Results Analysis"
echo "Time: $(date)"
echo "========================================"
echo ""

DIRS="MoSSe_Se_up__WSSe_Se_up MoSSe_Se_up__WSSe_S_up MoSSe_S_up__WSSe_Se_up MoSSe_S_up__WSSe_S_up"

for dir in $DIRS; do
    SCF_DIR="$dir/static_scf"
    
    if [ ! -d "$SCF_DIR" ]; then
        continue
    fi
    
    echo "========================================="
    echo "$dir"
    echo "========================================="
    
    cd "$SCF_DIR"
    
    # Check if calculation completed
    if [ ! -f "OUTCAR" ]; then
        echo "Status: Not started or running"
        cd ../..
        echo ""
        continue
    fi
    
    # Check convergence
    if grep -q "reached required accuracy" OUTCAR; then
        echo "Status: ✓ Converged"
    else
        echo "Status: ✗ Not converged"
    fi
    
    # Final energy
    if [ -f "OUTCAR" ]; then
        energy=$(grep "free  energy   TOTEN" OUTCAR | tail -1 | awk '{print $5}')
        echo "Final energy: $energy eV"
    fi
    
    # Fermi level
    if grep -q "E-fermi" OUTCAR; then
        efermi=$(grep "E-fermi" OUTCAR | tail -1 | awk '{print $3}')
        echo "Fermi level: $efermi eV"
    fi
    
    # Magnetization
    if grep -q "number of electron" OUTCAR; then
        mag=$(grep "number of electron" OUTCAR | tail -1 | awk '{print $NF}')
        echo "Magnetization: $mag μB"
    fi
    
    # Electronic iterations
    if [ -f "OSZICAR" ]; then
        iterations=$(grep "DAV:" OSZICAR | wc -l)
        echo "Electronic iterations: $iterations"
    fi
    
    # Output files
    echo ""
    echo "Output files:"
    for file in OUTCAR DOSCAR EIGENVAL CHGCAR WAVECAR; do
        if [ -f "$file" ]; then
            size=$(du -h "$file" | cut -f1)
            echo "  ✓ $file ($size)"
        else
            echo "  ✗ $file (missing)"
        fi
    done
    
    # Band gap estimate (rough)
    if [ -f "DOSCAR" ]; then
        echo ""
        echo "DOS available - use visualization tools for band gap"
        # Simple check around Fermi level
        # (This is very rough - proper analysis needs python/plotting)
    fi
    
    cd ../..
    echo ""
done

echo "========================================="
echo "Summary"
echo "========================================="
echo ""

# Count completed
completed=0
total=0
for dir in $DIRS; do
    if [ -f "$dir/static_scf/OUTCAR" ]; then
        total=$((total + 1))
        if grep -q "reached required accuracy" "$dir/static_scf/OUTCAR"; then
            completed=$((completed + 1))
        fi
    fi
done

echo "Completed: $completed / 4"
echo ""

if [ $completed -eq 4 ]; then
    echo "✓ All calculations complete!"
    echo ""
    echo "Next steps:"
    echo "  1. Extract DOS: vaspkit (option 11)"
    echo "  2. Band structure: Need band calculation"
    echo "  3. Charge analysis: Bader, etc."
    echo ""
else
    echo "Some calculations still running or not started"
    echo ""
    echo "Check status:"
    echo "  qstat -u \$USER"
    echo "  tail -f <directory>/static_scf/*.o*"
fi

echo "========================================"
