#!/bin/bash
# submit_all_jobs.sh
# Submit all VASP jobs for heterostructure calculations

echo "============================================================"
echo "Batch Job Submission Script"
echo "============================================================"
echo ""

# Find all directories with run.pbs
job_dirs=(MoSSe_*__WSSe_*/)

if [ ${#job_dirs[@]} -eq 0 ]; then
    echo "ERROR: No MoSSe_*__WSSe_* directories found!"
    exit 1
fi

echo "Found ${#job_dirs[@]} directories:"
for d in "${job_dirs[@]}"; do
    echo "  $d"
done
echo ""

# Check if run.pbs exists in each directory
echo "Checking for run.pbs scripts..."
all_exist=true
for d in "${job_dirs[@]}"; do
    if [ -f "$d/run.pbs" ]; then
        echo "  OK: $d/run.pbs"
    else
        echo "  MISSING: $d/run.pbs"
        all_exist=false
    fi
done
echo ""

if [ "$all_exist" = false ]; then
    echo "ERROR: Some directories are missing run.pbs"
    echo "Run: python generate_pbs_scripts.py"
    exit 1
fi

# Ask for confirmation
echo "Ready to submit ${#job_dirs[@]} jobs"
echo ""
read -p "Submit all jobs? (y/n) " -n 1 -r
echo ""

if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "Cancelled"
    exit 0
fi

# Submit jobs
echo ""
echo "Submitting jobs..."
echo ""

job_ids=()
for d in "${job_dirs[@]}"; do
    # Remove trailing slash
    dirname="${d%/}"
    
    echo "Submitting: $dirname"
    cd "$dirname" || continue
    
    # Submit job and capture job ID
    job_output=$(qsub run.pbs 2>&1)
    
    if [ $? -eq 0 ]; then
        job_id=$(echo "$job_output" | grep -oE '[0-9]+\.')
        job_ids+=("$job_id")
        echo "  Job ID: $job_output"
    else
        echo "  ERROR: $job_output"
    fi
    
    cd ..
    echo ""
done

echo "============================================================"
echo "Submission complete"
echo "============================================================"
echo ""
echo "Submitted ${#job_ids[@]} jobs:"
for jid in "${job_ids[@]}"; do
    echo "  $jid"
done
echo ""
echo "Check status:"
echo "  qstat -u \$USER"
echo ""
echo "Monitor specific job:"
echo "  qstat -f <job_id>"
echo ""
echo "Check output:"
echo "  tail -f MoSSe_Se_up__WSSe_S_up/*.o<jobid>"
echo "============================================================"
