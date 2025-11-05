#!/bin/bash
# diagnose_vasp_runs.sh
# Script to diagnose VASP calculation problems

echo "=========================================="
echo "VASP Calculation Diagnostics"
echo "=========================================="
echo ""

# Check problematic calculations
problem_dirs=(
    "MoSSe_Se_up__WSSe_S_up"
    "MoSSe_S_up__WSSe_Se_up"
)

for dir in "${problem_dirs[@]}"; do
    if [ ! -d "$dir" ]; then
        echo "Directory $dir not found, skipping..."
        continue
    fi
    
    echo "====== $dir ======"
    cd "$dir" || continue
    
    # 1. Check if files exist
    echo "1. File check:"
    for file in POSCAR INCAR KPOINTS POTCAR OUTCAR OSZICAR; do
        if [ -f "$file" ]; then
            size=$(stat -f%z "$file" 2>/dev/null || stat -c%s "$file" 2>/dev/null)
            echo "  $file: ${size} bytes"
        else
            echo "  $file: MISSING"
        fi
    done
    echo ""
    
    # 2. OSZICAR last lines
    echo "2. OSZICAR (last 5 lines):"
    if [ -f "OSZICAR" ]; then
        tail -5 OSZICAR | awk '{printf "  %s\n", $0}'
    else
        echo "  No OSZICAR found"
    fi
    echo ""
    
    # 3. Energy from OSZICAR
    echo "3. Energy evolution:"
    if [ -f "OSZICAR" ]; then
        grep "F=" OSZICAR | awk '{printf "  Step %s: E= %s eV\n", $1, $3}'
    fi
    echo ""
    
    # 4. Errors and warnings
    echo "4. Errors/Warnings in OUTCAR:"
    if [ -f "OUTCAR" ]; then
        grep -i "error\|warning\|fatal" OUTCAR | head -10 | awk '{printf "  %s\n", $0}'
        if [ $(grep -c -i "error\|warning\|fatal" OUTCAR) -eq 0 ]; then
            echo "  No errors/warnings found"
        fi
    else
        echo "  No OUTCAR found"
    fi
    echo ""
    
    # 5. ZBRENT error check
    echo "5. ZBRENT error check:"
    if [ -f "OUTCAR" ]; then
        if grep -q "ZBRENT" OUTCAR; then
            echo "  ZBRENT ERROR FOUND!"
            grep -A 2 "ZBRENT" OUTCAR | awk '{printf "  %s\n", $0}'
        else
            echo "  No ZBRENT error"
        fi
    fi
    echo ""
    
    # 6. Atomic distances from first ionic step
    echo "6. Initial atomic positions and forces:"
    if [ -f "OUTCAR" ]; then
        grep -A 10 "POSITION.*TOTAL-FORCE" OUTCAR | head -14 | awk '{printf "  %s\n", $0}'
    fi
    echo ""
    
    # 7. SCF convergence
    echo "7. SCF convergence (first ionic step):"
    if [ -f "OUTCAR" ]; then
        awk '/ITERATION/{p=1} p&&/DAV/{print "  "$0} /F=/{if(p)exit}' OUTCAR | head -15
    fi
    echo ""
    
    # 8. Check actual runtime
    echo "8. Calculation timing:"
    if [ -f "OUTCAR" ]; then
        echo -n "  Started: "
        head -20 OUTCAR | grep "executed on" | awk '{print $6, $7, $8}'
        echo -n "  Last update: "
        stat -f "%Sm" OUTCAR 2>/dev/null || stat -c "%y" OUTCAR 2>/dev/null
    fi
    echo ""
    
    cd ..
    echo "=========================================="
    echo ""
done

# Summary
echo ""
echo "====== Quick Summary ======"
for dir in "${problem_dirs[@]}"; do
    if [ -d "$dir" ]; then
        echo -n "$dir: "
        if [ -f "$dir/OSZICAR" ]; then
            last_line=$(tail -1 "$dir/OSZICAR")
            echo "$last_line"
        else
            echo "No OSZICAR"
        fi
    fi
done

echo ""
echo "====== Recommendations ======"
echo "If you see:"
echo "  - 'ZBRENT error' → Try ALGO=Normal, smaller AMIX"
echo "  - 'High energy (+400 eV)' → Check atom overlap in POSCAR"
echo "  - 'Only 1 step' → Structure problem or SCF not converging"
echo "  - 'WARNING: distance' → Atoms too close, increase delta_z"
echo ""
echo "Next steps:"
echo "  1. Check POSCAR z-coordinates manually"
echo "  2. Try with modified build_relax0_inputs_v2_fixed.py"
echo "  3. Consider 2-stage relaxation (coarse → fine)"
