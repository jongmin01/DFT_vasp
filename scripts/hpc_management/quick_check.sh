#!/bin/bash
# quick_check.sh
# Quick validation before submission (minimal version)

echo "========================================"
echo "QUICK PRE-SUBMISSION CHECK"
echo "========================================"
echo ""

# 1. Check z-coordinates (MOST CRITICAL)
echo "1. Checking z-coordinates (critical bug check)..."
for d in MoSSe_*__WSSe_*/; do
    nz=$(awk 'NR>8 && NF==3 {print $3}' "$d/POSCAR" | sort -u | wc -l)
    if [ "$nz" -lt 3 ]; then
        echo "  ✗ $d: ONLY $nz z-values - BUG DETECTED!"
    else
        echo "  ✓ $d: $nz z-values (OK)"
    fi
done
echo ""

# 2. Check species order
echo "2. Checking POSCAR species order..."
for d in MoSSe_*__WSSe_*/; do
    species=$(sed -n '6p' "$d/POSCAR")
    if echo "$species" | grep -q "Mo.*W.*S.*Se"; then
        echo "  ✓ $d: $species"
    else
        echo "  ✗ $d: Wrong order - $species"
    fi
done
echo ""

# 3. Check atom counts
echo "3. Checking atom counts..."
for d in MoSSe_*__WSSe_*/; do
    count=$(sed -n '7p' "$d/POSCAR" | awk '{sum=0; for(i=1;i<=NF;i++) sum+=$i; print sum}')
    if [ "$count" -eq 150 ]; then
        echo "  ✓ $d: $count atoms"
    else
        echo "  ! $d: $count atoms (expected 150)"
    fi
done
echo ""

# 4. Check all files exist
echo "4. Checking required files..."
all_ok=true
for d in MoSSe_*__WSSe_*/; do
    missing=""
    for f in POSCAR INCAR KPOINTS POTCAR run.pbs; do
        if [ ! -f "$d/$f" ]; then
            missing="$missing $f"
        fi
    done
    if [ -z "$missing" ]; then
        echo "  ✓ $d: All files present"
    else
        echo "  ✗ $d: Missing:$missing"
        all_ok=false
    fi
done
echo ""

# 5. Check POTCAR
echo "5. Checking POTCAR..."
if [ -f "POTCAR" ]; then
    ntitel=$(grep -c TITEL POTCAR)
    if [ "$ntitel" -eq 4 ]; then
        echo "  ✓ POTCAR has 4 elements"
        grep TITEL POTCAR | awk '{print "    ", $4}'
    else
        echo "  ✗ POTCAR has $ntitel elements (expected 4)"
    fi
else
    echo "  ✗ POTCAR not found in current directory"
fi
echo ""

# Summary
echo "========================================"
if [ "$all_ok" = true ]; then
    echo "READY TO SUBMIT!"
    echo ""
    echo "Next steps:"
    echo "  1. Test one job:"
    echo "     cd MoSSe_Se_up__WSSe_S_up"
    echo "     qsub run.pbs"
    echo ""
    echo "  2. After 5-10 min, check:"
    echo "     tail -20 MoSSe_Se_up__WSSe_S_up/OSZICAR"
    echo ""
    echo "  3. If OK, submit all:"
    echo "     ./submit_all_jobs.sh"
else
    echo "ISSUES FOUND - REVIEW ABOVE"
    echo ""
    echo "Run detailed check:"
    echo "  ./final_check_before_submit.sh"
fi
echo "========================================"
