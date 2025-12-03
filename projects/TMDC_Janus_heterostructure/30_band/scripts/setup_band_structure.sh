#!/bin/bash
# setup_band_structure.sh
# Setup band structure calculations from static SCF results

echo "========================================"
echo "Band Structure Setup (SOC OFF)"
echo "========================================"
echo ""

# Target structure (most stable)
STRUCTURE="MoSSe_Se_up__WSSe_S_up"

# Check if required files exist
echo "Checking setup files..."
REQUIRED_FILES="INCAR_band_soc_off KPOINTS_band run_band_soc_off.pbs"
ALL_FILES_OK=1

for file in $REQUIRED_FILES; do
    if [ -f "$file" ]; then
        echo "  ✓ $file"
    else
        echo "  ✗ $file NOT FOUND!"
        ALL_FILES_OK=0
    fi
done
echo ""

if [ $ALL_FILES_OK -eq 0 ]; then
    echo "ERROR: Required files missing!"
    exit 1
fi

# Check if structure directory exists
if [ ! -d "$STRUCTURE" ]; then
    echo "✗ Directory not found: $STRUCTURE"
    exit 1
fi

echo "========================================="
echo "Setting up: $STRUCTURE"
echo "========================================="

# Check if static_scf exists
if [ ! -d "$STRUCTURE/static_scf" ]; then
    echo "  ✗ ERROR: static_scf directory not found!"
    echo "    Band structure needs static SCF results"
    exit 1
fi

# Check for required files from static SCF
echo "  Checking static SCF results..."
REQUIRED_SCF="POSCAR POTCAR CHGCAR"
SCF_OK=1

for file in $REQUIRED_SCF; do
    if [ -f "$STRUCTURE/static_scf/$file" ]; then
        echo "    ✓ $file"
    else
        echo "    ✗ $file missing!"
        SCF_OK=0
    fi
done

if [ $SCF_OK -eq 0 ]; then
    echo "  ✗ ERROR: Required static SCF files missing!"
    exit 1
fi

# Create band_soc_off subdirectory
mkdir -p "$STRUCTURE/band_soc_off"
cd "$STRUCTURE/band_soc_off"

echo ""
echo "  Creating band_soc_off directory..."

# Copy files from static SCF
echo "  Copying files from static_scf..."
cp ../static_scf/POSCAR .
cp ../static_scf/POTCAR .
cp ../static_scf/CHGCAR .

# Copy new input files
echo "  Setting up band structure inputs..."
cp ../../INCAR_band_soc_off INCAR
cp ../../KPOINTS_band KPOINTS
cp ../../run_band_soc_off.pbs .

# Verify all files
echo ""
echo "  Verifying files..."
ALL_OK=1
for file in POSCAR INCAR KPOINTS POTCAR CHGCAR run_band_soc_off.pbs; do
    if [ -f "$file" ]; then
        size=$(stat -c%s "$file" 2>/dev/null || stat -f%z "$file" 2>/dev/null)
        echo "    ✓ $file ($size bytes)"
    else
        echo "    ✗ $file MISSING!"
        ALL_OK=0
    fi
done

cd ../..

if [ $ALL_OK -eq 1 ]; then
    echo ""
    echo "  ✓ Setup complete for $STRUCTURE/band_soc_off"
else
    echo ""
    echo "  ✗ Setup FAILED for $STRUCTURE"
    exit 1
fi

echo ""
echo "========================================="
echo "Setup Complete!"
echo "========================================="
echo ""

# Summary
echo "Directory structure:"
echo "$STRUCTURE/"
echo "  ├── stage1_results/     (relaxation)"
echo "  ├── static_scf/         (DOS, charge density)"
echo "  └── band_soc_off/       (NEW - band structure)"
echo "      ├── POSCAR          (from static_scf)"
echo "      ├── INCAR           (ICHARG=11, k-path)"
echo "      ├── KPOINTS         (Γ-M-K-Γ path)"
echo "      ├── POTCAR          (from static_scf)"
echo "      ├── CHGCAR          (from static_scf)"
echo "      └── run_band_soc_off.pbs"
echo ""

# Ask about job submission
read -p "Submit band structure job now? (y/n): " submit_choice
echo ""

if [ "$submit_choice" = "y" ] || [ "$submit_choice" = "Y" ]; then
    echo "Submitting job..."
    cd "$STRUCTURE/band_soc_off"
    
    JOB_ID=$(qsub run_band_soc_off.pbs)
    
    if [ $? -eq 0 ]; then
        echo "  ✓ Job submitted successfully!"
        echo "  Job ID: $JOB_ID"
        echo ""
        echo "Monitor with:"
        echo "  qstat -u \$USER"
        echo "  tail -f $STRUCTURE/band_soc_off/*.o*"
    else
        echo "  ✗ Job submission FAILED"
    fi
    
    cd ../..
else
    echo "Job NOT submitted."
    echo ""
    echo "To submit manually:"
    echo "  cd $STRUCTURE/band_soc_off"
    echo "  qsub run_band_soc_off.pbs"
fi

echo ""
echo "========================================="
echo "Expected calculation time: 4-6 hours"
echo "========================================="
echo ""

# Print k-path info
echo "K-path information:"
echo "  Γ-M-K-Γ (standard for hexagonal 2D)"
echo "  40 points per segment"
echo "  Total: ~120 k-points"
echo ""

echo "After completion:"
echo "  1. Download EIGENVAL file"
echo "  2. Plot band structure"
echo "  3. Analyze band gap (direct/indirect)"
echo "  4. Compare with DOS"
echo ""

echo "Next step: SOC ON calculation"
echo "  (After SOC OFF completes successfully)"
