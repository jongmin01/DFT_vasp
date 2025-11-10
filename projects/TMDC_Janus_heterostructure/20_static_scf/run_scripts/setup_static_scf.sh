#!/bin/bash
# setup_static_scf.sh
# Setup static SCF calculations from relaxed structures

echo "========================================"
echo "Static SCF Calculation Setup"
echo "========================================"
echo ""

DIRS="MoSSe_Se_up__WSSe_Se_up MoSSe_Se_up__WSSe_S_up MoSSe_S_up__WSSe_Se_up MoSSe_S_up__WSSe_S_up"

# Check if required files exist
echo "Checking required setup files..."
REQUIRED_FILES="INCAR_static_scf KPOINTS_static_scf run_static_scf.pbs"
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

# Process each directory
for dir in $DIRS; do
    if [ ! -d "$dir" ]; then
        echo "✗ Directory not found: $dir"
        continue
    fi
    
    echo "========================================="
    echo "Setting up: $dir"
    echo "========================================="
    
    # Create static_scf subdirectory
    mkdir -p "$dir/static_scf"
    cd "$dir/static_scf"
    
    # Check if relaxed structure exists
    if [ -f "../CONTCAR_final" ]; then
        echo "  ✓ Found relaxed structure"
        cp ../CONTCAR_final POSCAR
    elif [ -f "../stage1_results/CONTCAR" ]; then
        echo "  ✓ Using Stage 1 CONTCAR"
        cp ../stage1_results/CONTCAR POSCAR
    else
        echo "  ✗ ERROR: No relaxed structure found!"
        cd ../..
        continue
    fi
    
    # Copy input files
    echo "  Copying input files..."
    cp ../../INCAR_static_scf INCAR
    cp ../../KPOINTS_static_scf KPOINTS
    cp ../../run_static_scf.pbs .
    
    # Copy POTCAR from parent directory
    if [ -f "../POTCAR" ]; then
        cp ../POTCAR .
        echo "    ✓ POTCAR"
    else
        echo "    ✗ POTCAR not found!"
        cd ../..
        continue
    fi
    
    # Verify all files
    echo "  Verifying files..."
    ALL_OK=1
    for file in POSCAR INCAR KPOINTS POTCAR run_static_scf.pbs; do
        if [ -f "$file" ]; then
            echo "    ✓ $file"
        else
            echo "    ✗ $file MISSING!"
            ALL_OK=0
        fi
    done
    
    if [ $ALL_OK -eq 1 ]; then
        echo "  ✓ Setup complete for $dir/static_scf"
    else
        echo "  ✗ Setup FAILED for $dir"
    fi
    
    cd ../..
    echo ""
done

echo "========================================="
echo "Setup Complete!"
echo "========================================="
echo ""

# Ask about job submission
read -p "Submit jobs now? (y/n): " submit_choice
echo ""

if [ "$submit_choice" = "y" ] || [ "$submit_choice" = "Y" ]; then
    echo "Submitting jobs..."
    echo ""
    
    for dir in $DIRS; do
        if [ -f "$dir/static_scf/run_static_scf.pbs" ]; then
            cd "$dir/static_scf"
            
            JOB_ID=$(qsub run_static_scf.pbs)
            
            if [ $? -eq 0 ]; then
                echo "  ✓ $dir: Job submitted (ID: $JOB_ID)"
            else
                echo "  ✗ $dir: Job submission FAILED"
            fi
            
            cd ../..
        fi
    done
    
    echo ""
    echo "========================================="
    echo "Jobs submitted!"
    echo "========================================="
    echo ""
    echo "Check job status with:"
    echo "  qstat -u \$USER"
    echo ""
    echo "Monitor progress:"
    echo "  tail -f <directory>/static_scf/*.o*"
    echo ""
else
    echo "Jobs NOT submitted."
    echo ""
    echo "To submit manually:"
    echo "  cd <directory>/static_scf"
    echo "  qsub run_static_scf.pbs"
    echo ""
fi

echo "========================================="
echo "Directory structure:"
echo "========================================="
for dir in $DIRS; do
    if [ -d "$dir/static_scf" ]; then
        echo "$dir/"
        echo "  ├── stage1_results/     (relaxation)"
        echo "  ├── CONTCAR_final       (relaxed structure)"
        echo "  └── static_scf/         (NEW - static calculation)"
        echo "      ├── POSCAR"
        echo "      ├── INCAR"
        echo "      ├── KPOINTS"
        echo "      ├── POTCAR"
        echo "      └── run_static_scf.pbs"
        echo ""
    fi
done

echo "Expected calculation time: 8-12 hours per structure"
echo "Total: ~12 hours (parallel execution)"
