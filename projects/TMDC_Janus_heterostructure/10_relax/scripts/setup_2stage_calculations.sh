#!/bin/bash
# setup_2stage_calculations.sh
# Complete setup script for 2-stage VASP calculations

echo "========================================"
echo "2-Stage VASP Calculation Setup Script"
echo "========================================"
echo ""

# Target directories
DIRS="MoSSe_Se_up__WSSe_S_up MoSSe_S_up__WSSe_Se_up MoSSe_Se_up__WSSe_Se_up MoSSe_S_up__WSSe_S_up"

# Check if required files exist in current directory
echo "Step 1: Checking required setup files..."
REQUIRED_FILES="run_2stage.pbs INCAR_step1 INCAR_step2"
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
    echo "Please ensure these files are in the current directory:"
    echo "  - run_2stage.pbs"
    echo "  - INCAR_step1"
    echo "  - INCAR_step2"
    exit 1
fi

# Check if directories exist
echo "Step 2: Checking target directories..."
DIR_EXISTS=0
for dir in $DIRS; do
    if [ -d "$dir" ]; then
        echo "  ✓ $dir"
        DIR_EXISTS=1
    else
        echo "  ✗ $dir not found"
    fi
done
echo ""

if [ $DIR_EXISTS -eq 0 ]; then
    echo "ERROR: No target directories found!"
    exit 1
fi

# Ask for confirmation
echo "This script will:"
echo "  1. Copy run_2stage.pbs, INCAR_step1, INCAR_step2 to each directory"
echo "  2. Optionally clean old output files"
echo "  3. Submit jobs to the queue"
echo ""
read -p "Continue? (y/n): " continue_choice

if [ "$continue_choice" != "y" ] && [ "$continue_choice" != "Y" ]; then
    echo "Setup cancelled."
    exit 0
fi
echo ""

# Ask about cleaning old files
read -p "Clean old VASP output files? (y/n): " clean_choice
echo ""

# Process each directory
for dir in $DIRS; do
    if [ ! -d "$dir" ]; then
        continue
    fi
    
    echo "========================================"
    echo "Setting up: $dir"
    echo "========================================"
    cd "$dir"
    
    # Clean old files if requested
    if [ "$clean_choice" = "y" ] || [ "$clean_choice" = "Y" ]; then
        echo "  Cleaning old output files..."
        
        # Create backup directory
        if [ -f OUTCAR ]; then
            BACKUP_DIR="old_backup_$(date +%Y%m%d_%H%M%S)"
            mkdir -p "$BACKUP_DIR"
            
            # Move old files to backup
            mv OUTCAR OSZICAR CONTCAR XDATCAR \
               WAVECAR CHGCAR CHG vasprun.xml \
               EIGENVAL DOSCAR PROCAR \
               *.out *.log IBZKPT PCDAT \
               "$BACKUP_DIR/" 2>/dev/null
            
            echo "    ✓ Old files moved to: $BACKUP_DIR"
        else
            echo "    No old output files found"
        fi
    fi
    
    # Copy new files
    echo "  Copying setup files..."
    cp ../run_2stage.pbs .
    cp ../INCAR_step1 .
    cp ../INCAR_step2 .
    echo "    ✓ run_2stage.pbs"
    echo "    ✓ INCAR_step1"
    echo "    ✓ INCAR_step2"
    
    # Verify input files
    echo "  Verifying input files..."
    INPUT_OK=1
    for file in POSCAR KPOINTS POTCAR; do
        if [ -f "$file" ]; then
            echo "    ✓ $file"
        else
            echo "    ✗ $file MISSING!"
            INPUT_OK=0
        fi
    done
    
    if [ $INPUT_OK -eq 0 ]; then
        echo "  ✗ ERROR: Missing input files in $dir"
        cd ..
        continue
    fi
    
    echo "  ✓ Setup complete for $dir"
    echo ""
    cd ..
done

echo "========================================"
echo "Setup Complete!"
echo "========================================"
echo ""

# Ask about job submission
read -p "Submit jobs now? (y/n): " submit_choice
echo ""

if [ "$submit_choice" = "y" ] || [ "$submit_choice" = "Y" ]; then
    echo "Submitting jobs..."
    echo ""
    
    for dir in $DIRS; do
        if [ -d "$dir" ] && [ -f "$dir/run_2stage.pbs" ]; then
            cd "$dir"
            
            # Update job name in PBS script
            sed -i "s/#PBS -N .*/#PBS -N ${dir}/" run_2stage.pbs
            
            # Submit job
            JOB_ID=$(qsub run_2stage.pbs)
            
            if [ $? -eq 0 ]; then
                echo "  ✓ $dir: Job submitted (ID: $JOB_ID)"
            else
                echo "  ✗ $dir: Job submission FAILED"
            fi
            
            cd ..
        fi
    done
    
    echo ""
    echo "========================================"
    echo "Jobs submitted!"
    echo "========================================"
    echo ""
    echo "Check job status with:"
    echo "  qstat -u \$USER"
    echo ""
    echo "Monitor progress with:"
    echo "  tail -f <directory>/<jobname>.o*"
    echo "  grep 'DAV:' <directory>/OSZICAR | tail -20"
    echo ""
else
    echo "Jobs NOT submitted."
    echo ""
    echo "To submit manually, run in each directory:"
    echo "  cd <directory>"
    echo "  qsub run_2stage.pbs"
    echo ""
fi

echo "========================================"
echo "Next steps:"
echo "========================================"
echo "1. Monitor jobs: qstat -u \$USER"
echo "2. Check Stage 1 progress: grep 'Stage 1' */run_2stage.pbs.o*"
echo "3. Check convergence: grep 'CONVERGENCE' */run_2stage.pbs.o*"
echo ""
echo "Expected timeline:"
echo "  Stage 1: 12-24 hours"
echo "  Stage 2: 24-48 hours"
echo "  Total: 36-72 hours"
echo ""
