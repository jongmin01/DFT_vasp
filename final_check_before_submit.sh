#!/bin/bash
# final_check_before_submit.sh
# Comprehensive pre-submission validation script

echo "============================================================"
echo "FINAL PRE-SUBMISSION CHECK"
echo "$(date)"
echo "============================================================"
echo ""

# Color codes for output (if terminal supports it)
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Counters
total_checks=0
passed_checks=0
failed_checks=0
warning_checks=0

# Function to print status
print_status() {
    local status=$1
    local message=$2
    
    total_checks=$((total_checks + 1))
    
    if [ "$status" = "PASS" ]; then
        echo -e "${GREEN}[✓]${NC} $message"
        passed_checks=$((passed_checks + 1))
    elif [ "$status" = "FAIL" ]; then
        echo -e "${RED}[✗]${NC} $message"
        failed_checks=$((failed_checks + 1))
    elif [ "$status" = "WARN" ]; then
        echo -e "${YELLOW}[!]${NC} $message"
        warning_checks=$((warning_checks + 1))
    fi
}

# ============================================================
# CHECK 1: DIRECTORY STRUCTURE
# ============================================================
echo "1. Checking directory structure..."
echo ""

calc_dirs=(MoSSe_*__WSSe_*/)

if [ ${#calc_dirs[@]} -eq 0 ]; then
    print_status "FAIL" "No calculation directories found (MoSSe_*__WSSe_*)"
    echo ""
    echo "CRITICAL: Cannot proceed without calculation directories"
    exit 1
else
    print_status "PASS" "Found ${#calc_dirs[@]} calculation directories"
fi

expected_dirs=("MoSSe_S_up__WSSe_S_up" "MoSSe_S_up__WSSe_Se_up" "MoSSe_Se_up__WSSe_S_up" "MoSSe_Se_up__WSSe_Se_up")
for dir in "${expected_dirs[@]}"; do
    if [ -d "$dir" ]; then
        print_status "PASS" "Directory exists: $dir"
    else
        print_status "WARN" "Expected directory not found: $dir"
    fi
done

echo ""

# ============================================================
# CHECK 2: REQUIRED INPUT FILES
# ============================================================
echo "2. Checking required input files in each directory..."
echo ""

required_files=("POSCAR" "INCAR" "KPOINTS" "POTCAR" "run.pbs")

for d in "${calc_dirs[@]}"; do
    dirname="${d%/}"
    echo "  Checking: $dirname"
    
    all_files_present=true
    for file in "${required_files[@]}"; do
        if [ -f "$d/$file" ]; then
            size=$(stat -f%z "$d/$file" 2>/dev/null || stat -c%s "$d/$file" 2>/dev/null)
            if [ "$size" -gt 100 ]; then
                echo "    ✓ $file (${size} bytes)"
            else
                print_status "FAIL" "$dirname: $file exists but is too small (${size} bytes)"
                all_files_present=false
            fi
        else
            print_status "FAIL" "$dirname: Missing $file"
            all_files_present=false
        fi
    done
    
    if [ "$all_files_present" = true ]; then
        print_status "PASS" "$dirname: All required files present"
    fi
    echo ""
done

# ============================================================
# CHECK 3: POSCAR VALIDATION (CRITICAL - Previous Bug)
# ============================================================
echo "3. Validating POSCAR files (checking for z-coordinate bug)..."
echo ""

for d in "${calc_dirs[@]}"; do
    dirname="${d%/}"
    
    if [ ! -f "$d/POSCAR" ]; then
        print_status "FAIL" "$dirname: POSCAR not found"
        continue
    fi
    
    # Check atom counts
    species_line=$(sed -n '6p' "$d/POSCAR")
    counts_line=$(sed -n '7p' "$d/POSCAR")
    
    echo "  $dirname:"
    echo "    Species: $species_line"
    echo "    Counts: $counts_line"
    
    # Check if species order is correct
    if echo "$species_line" | grep -q "Mo.*W.*S.*Se"; then
        print_status "PASS" "$dirname: Species order correct (Mo W S Se)"
    else
        print_status "FAIL" "$dirname: Species order incorrect: $species_line"
    fi
    
    # CRITICAL: Check for z-coordinate bug (all atoms same z)
    # Extract z-coordinates (3rd column, skip first 8 lines)
    unique_z=$(awk 'NR>8 && NF==3 {print $3}' "$d/POSCAR" | sort -u | wc -l)
    
    if [ "$unique_z" -lt 3 ]; then
        print_status "FAIL" "$dirname: CRITICAL BUG - Only $unique_z unique z-coordinates (atoms overlapping!)"
        echo "    This is the bug we fixed! All atoms have same z-coordinate."
        echo "    Sample z-values:"
        awk 'NR>8 && NR<13 && NF==3 {print "      ", $3}' "$d/POSCAR"
    elif [ "$unique_z" -lt 4 ]; then
        print_status "WARN" "$dirname: Only $unique_z unique z-coordinates (expected 4-6)"
    else
        print_status "PASS" "$dirname: $unique_z unique z-coordinates (structure looks good)"
        echo "    Sample z-values (first 5):"
        awk 'NR>8 && NF==3 {print $3}' "$d/POSCAR" | sort -u | head -5 | sed 's/^/      /'
    fi
    
    # Check atom count (should be 150 for 5x5 supercell)
    total_atoms=$(awk 'NR==7 {sum=0; for(i=1;i<=NF;i++) sum+=$i; print sum}' "$d/POSCAR")
    if [ "$total_atoms" -eq 150 ]; then
        print_status "PASS" "$dirname: Correct atom count (150)"
    else
        print_status "WARN" "$dirname: Unexpected atom count ($total_atoms, expected 150)"
    fi
    
    echo ""
done

# ============================================================
# CHECK 4: POTCAR VALIDATION
# ============================================================
echo "4. Validating POTCAR files..."
echo ""

for d in "${calc_dirs[@]}"; do
    dirname="${d%/}"
    
    if [ ! -f "$d/POTCAR" ]; then
        print_status "FAIL" "$dirname: POTCAR not found"
        continue
    fi
    
    # Check POTCAR species order
    potcar_species=$(grep TITEL "$d/POTCAR" | awk '{print $4}' | tr '\n' ' ')
    echo "  $dirname POTCAR species: $potcar_species"
    
    # Count TITEL lines (should be 4: Mo, W, S, Se)
    titel_count=$(grep -c TITEL "$d/POTCAR")
    if [ "$titel_count" -eq 4 ]; then
        print_status "PASS" "$dirname: POTCAR has 4 elements"
    else
        print_status "FAIL" "$dirname: POTCAR has $titel_count elements (expected 4)"
    fi
    
    # Check if order matches POSCAR
    poscar_species=$(sed -n '6p' "$d/POSCAR" | tr -d ' ')
    potcar_order=$(grep TITEL "$d/POTCAR" | awk '{print $4}' | sed 's/_pv//g' | sed 's/_sv//g' | tr '\n' ' ' | tr -d ' ')
    
    # Simplified check (just ensure Mo, W, S, Se are present)
    if echo "$potcar_species" | grep -q "Mo" && \
       echo "$potcar_species" | grep -q "W" && \
       echo "$potcar_species" | grep -q "S" && \
       echo "$potcar_species" | grep -q "Se"; then
        print_status "PASS" "$dirname: POTCAR contains all required elements"
    else
        print_status "FAIL" "$dirname: POTCAR missing required elements"
    fi
    
    echo ""
done

# ============================================================
# CHECK 5: INCAR SETTINGS
# ============================================================
echo "5. Checking INCAR settings..."
echo ""

critical_incar_tags=("ENCUT" "IBRION" "NSW" "EDIFFG" "PREC" "ALGO")

for d in "${calc_dirs[@]}"; do
    dirname="${d%/}"
    
    if [ ! -f "$d/INCAR" ]; then
        print_status "FAIL" "$dirname: INCAR not found"
        continue
    fi
    
    echo "  $dirname INCAR:"
    
    all_tags_present=true
    for tag in "${critical_incar_tags[@]}"; do
        if grep -q "^[[:space:]]*$tag" "$d/INCAR"; then
            value=$(grep "^[[:space:]]*$tag" "$d/INCAR" | head -1 | awk -F= '{print $2}' | tr -d ' ')
            echo "    $tag = $value"
        else
            print_status "WARN" "$dirname: $tag not found in INCAR"
            all_tags_present=false
        fi
    done
    
    if [ "$all_tags_present" = true ]; then
        print_status "PASS" "$dirname: All critical INCAR tags present"
    fi
    
    # Check for problematic settings
    if grep -q "LREAL.*Auto" "$d/INCAR"; then
        print_status "WARN" "$dirname: LREAL=Auto (consider .FALSE. for accuracy)"
    fi
    
    echo ""
done

# ============================================================
# CHECK 6: KPOINTS
# ============================================================
echo "6. Checking KPOINTS..."
echo ""

for d in "${calc_dirs[@]}"; do
    dirname="${d%/}"
    
    if [ ! -f "$d/KPOINTS" ]; then
        print_status "FAIL" "$dirname: KPOINTS not found"
        continue
    fi
    
    # Check k-mesh
    kmesh=$(sed -n '4p' "$d/KPOINTS")
    echo "  $dirname k-mesh: $kmesh"
    
    # For 5x5 supercell, 2x2x1 or 3x3x1 is reasonable
    if echo "$kmesh" | grep -qE "[2-3][[:space:]]+[2-3][[:space:]]+1"; then
        print_status "PASS" "$dirname: K-mesh looks reasonable for supercell"
    else
        print_status "WARN" "$dirname: Unusual k-mesh for 5x5 supercell: $kmesh"
    fi
    
    echo ""
done

# ============================================================
# CHECK 7: PBS SCRIPT VALIDATION
# ============================================================
echo "7. Checking PBS scripts (run.pbs)..."
echo ""

for d in "${calc_dirs[@]}"; do
    dirname="${d%/}"
    
    if [ ! -f "$d/run.pbs" ]; then
        print_status "FAIL" "$dirname: run.pbs not found"
        continue
    fi
    
    # Check if executable
    if [ -x "$d/run.pbs" ]; then
        print_status "PASS" "$dirname: run.pbs is executable"
    else
        print_status "WARN" "$dirname: run.pbs not executable (will still work with qsub)"
    fi
    
    # Check for critical PBS directives
    pbs_checks=("#PBS -N" "#PBS -A" "#PBS -l nodes" "#PBS -l walltime" "mpirun")
    
    all_pbs_ok=true
    for check in "${pbs_checks[@]}"; do
        if ! grep -q "$check" "$d/run.pbs"; then
            print_status "WARN" "$dirname: Missing '$check' in run.pbs"
            all_pbs_ok=false
        fi
    done
    
    if [ "$all_pbs_ok" = true ]; then
        print_status "PASS" "$dirname: run.pbs has all critical PBS directives"
    fi
    
    # Extract and show PBS settings
    echo "    PBS settings:"
    grep "#PBS -N" "$d/run.pbs" | head -1 | sed 's/^/      /'
    grep "#PBS -l nodes" "$d/run.pbs" | head -1 | sed 's/^/      /'
    grep "#PBS -l walltime" "$d/run.pbs" | head -1 | sed 's/^/      /'
    
    echo ""
done

# ============================================================
# CHECK 8: PREVIOUS CALCULATION CLEANUP
# ============================================================
echo "8. Checking for previous calculation outputs..."
echo ""

for d in "${calc_dirs[@]}"; do
    dirname="${d%/}"
    
    # Check for output files that would interfere
    conflict_files=("WAVECAR" "CHGCAR" "OUTCAR" "OSZICAR" "vasprun.xml")
    found_conflicts=false
    
    for file in "${conflict_files[@]}"; do
        if [ -f "$d/$file" ]; then
            if [ ! "$found_conflicts" = true ]; then
                echo "  $dirname has previous output files:"
                found_conflicts=true
            fi
            size=$(stat -f%z "$d/$file" 2>/dev/null || stat -c%s "$d/$file" 2>/dev/null)
            echo "    - $file (${size} bytes)"
        fi
    done
    
    if [ "$found_conflicts" = true ]; then
        print_status "WARN" "$dirname: Previous outputs exist (may want to backup/clean)"
    else
        print_status "PASS" "$dirname: No previous outputs (clean start)"
    fi
    
    echo ""
done

# ============================================================
# CHECK 9: DISK SPACE
# ============================================================
echo "9. Checking available disk space..."
echo ""

# Check current directory disk space
df_output=$(df -h . | tail -1)
available=$(echo "$df_output" | awk '{print $4}')
percent_used=$(echo "$df_output" | awk '{print $5}' | tr -d '%')

echo "  Current directory: $df_output"
echo "  Available: $available"

if [ "$percent_used" -gt 90 ]; then
    print_status "FAIL" "Disk usage over 90% - may run out of space during calculation"
elif [ "$percent_used" -gt 80 ]; then
    print_status "WARN" "Disk usage over 80% - monitor space during calculation"
else
    print_status "PASS" "Sufficient disk space available"
fi

echo ""

# ============================================================
# CHECK 10: MODULE AVAILABILITY
# ============================================================
echo "10. Checking VASP module availability..."
echo ""

# Check if module command exists
if command -v module &> /dev/null; then
    print_status "PASS" "Module system available"
    
    # Check for VASP modules
    echo "  Available VASP modules:"
    module avail vasp 2>&1 | grep -i vasp | sed 's/^/    /' || echo "    (none found - check module name)"
    echo ""
else
    print_status "WARN" "Module system not available (may not be on compute node)"
fi

# ============================================================
# SUMMARY
# ============================================================
echo ""
echo "============================================================"
echo "VALIDATION SUMMARY"
echo "============================================================"
echo "Total checks: $total_checks"
echo "  Passed:  $passed_checks"
echo "  Failed:  $failed_checks"
echo "  Warnings: $warning_checks"
echo ""

if [ "$failed_checks" -eq 0 ]; then
    echo -e "${GREEN}✓ ALL CRITICAL CHECKS PASSED${NC}"
    echo ""
    echo "You are ready to submit!"
    echo ""
    echo "Next steps:"
    echo "  1. Review this output carefully"
    echo "  2. Submit test job: cd MoSSe_Se_up__WSSe_S_up && qsub run.pbs"
    echo "  3. Monitor: qstat -u \$USER"
    echo "  4. If test successful, submit all: ./submit_all_jobs.sh"
    echo ""
    exit_code=0
else
    echo -e "${RED}✗ $failed_checks CRITICAL ISSUE(S) FOUND${NC}"
    echo ""
    echo "Please fix the failed checks before submitting!"
    echo ""
    exit_code=1
fi

if [ "$warning_checks" -gt 0 ]; then
    echo -e "${YELLOW}Note: $warning_checks warning(s) found${NC}"
    echo "Warnings may be acceptable but should be reviewed."
    echo ""
fi

echo "============================================================"
exit $exit_code
