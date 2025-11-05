#!/bin/bash
# submit_all_jobs.sh
# ============================================================
# Batch submission script for VASP relaxation jobs
# Usage: ./submit_all_jobs.sh [base_directory]
# ============================================================

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Get base directory
BASE_DIR="${1:-.}"

if [ ! -d "$BASE_DIR" ]; then
    echo -e "${RED}Error: Directory not found: $BASE_DIR${NC}"
    echo "Usage: $0 [base_directory]"
    exit 1
fi

echo "============================================================"
echo "BATCH JOB SUBMISSION FOR VASP RELAXATION"
echo "============================================================"
echo "Base directory: $BASE_DIR"
echo ""

# Find all subdirectories with submit.pbs
CALC_DIRS=()
for dir in "$BASE_DIR"/*; do
    if [ -d "$dir" ] && [ -f "$dir/submit.pbs" ]; then
        CALC_DIRS+=("$dir")
    fi
done

echo "Found ${#CALC_DIRS[@]} calculation directories"
echo ""

if [ ${#CALC_DIRS[@]} -eq 0 ]; then
    echo -e "${RED}No calculation directories found!${NC}"
    echo "Looking for directories with submit.pbs"
    exit 1
fi

# List directories
echo "Directories to submit:"
echo "------------------------------------------------------------"
for i in "${!CALC_DIRS[@]}"; do
    dir_name=$(basename "${CALC_DIRS[$i]}")
    echo "  $((i+1)). $dir_name"
done
echo "------------------------------------------------------------"
echo ""

# Ask for confirmation
read -p "Submit all jobs? (y/n): " -n 1 -r
echo ""

if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "Submission cancelled"
    exit 0
fi

# Submit jobs
echo ""
echo "Submitting jobs..."
echo "============================================================"

SUBMITTED=0
FAILED=0
JOB_IDS=()

for calc_dir in "${CALC_DIRS[@]}"; do
    dir_name=$(basename "$calc_dir")
    pbs_file="$calc_dir/submit.pbs"
    
    echo -n "Submitting $dir_name... "
    
    # Change to directory and submit
    cd "$calc_dir" || { 
        echo -e "${RED}Failed (cannot cd)${NC}"
        ((FAILED++))
        continue
    }
    
    # Submit job
    output=$(qsub submit.pbs 2>&1)
    exit_code=$?
    
    if [ $exit_code -eq 0 ]; then
        job_id=$(echo "$output" | awk '{print $1}')
        echo -e "${GREEN}OK${NC} (Job ID: $job_id)"
        JOB_IDS+=("$job_id:$dir_name")
        ((SUBMITTED++))
    else
        echo -e "${RED}Failed${NC}"
        echo "  Error: $output"
        ((FAILED++))
    fi
    
    cd - > /dev/null || exit
done

echo "============================================================"
echo "SUBMISSION SUMMARY"
echo "============================================================"
echo "Submitted: $SUBMITTED"
echo "Failed: $FAILED"
echo ""

if [ $SUBMITTED -gt 0 ]; then
    echo "Job IDs:"
    echo "------------------------------------------------------------"
    for job_info in "${JOB_IDS[@]}"; do
        job_id="${job_info%%:*}"
        job_dir="${job_info##*:}"
        echo "  $job_id  $job_dir"
    done
    echo "============================================================"
    echo ""
    
    # Create job tracking file
    tracking_file="$BASE_DIR/submitted_jobs.txt"
    echo "# Submitted jobs - $(date)" > "$tracking_file"
    echo "# Job_ID Directory" >> "$tracking_file"
    for job_info in "${JOB_IDS[@]}"; do
        job_id="${job_info%%:*}"
        job_dir="${job_info##*:}"
        echo "$job_id $job_dir" >> "$tracking_file"
    done
    echo -e "${GREEN}Job tracking file created: $tracking_file${NC}"
    echo ""
    
    # Show queue status
    echo "Current queue status:"
    echo "------------------------------------------------------------"
    qstat -u $USER
    echo "============================================================"
    echo ""
    
    # Monitoring tips
    echo -e "${BLUE}MONITORING COMMANDS:${NC}"
    echo "------------------------------------------------------------"
    echo "Check all jobs:"
    echo "  qstat -u \$USER"
    echo ""
    echo "Check specific job:"
    echo "  qstat -f <job_id>"
    echo "  qstat -n <job_id>  # show nodes"
    echo ""
    echo "Check progress:"
    echo "  tail -f [directory]/vasp.log"
    echo "  grep 'F=' [directory]/OSZICAR | tail -10"
    echo ""
    echo "Quick convergence check:"
    echo "  grep 'reached required accuracy' */OUTCAR"
    echo ""
    echo "Kill a job:"
    echo "  qdel <job_id>"
    echo "============================================================"
fi

echo ""
echo "Done!"
