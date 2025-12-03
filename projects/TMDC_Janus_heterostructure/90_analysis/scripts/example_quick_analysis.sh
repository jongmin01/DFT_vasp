#!/bin/bash
#
# example_quick_analysis.sh - Quick analysis example for MoSSe/WSSe heterostructure
#
# This script demonstrates how to run PROCAR analysis after copying VASP output files
# from HPC calculation to local machine.
#
# Usage:
#   1. Copy VASP output files to a data directory
#   2. Run this script with the data directory as argument
#
# Example:
#   bash example_quick_analysis.sh /path/to/vasp/outputs

# Exit on error
set -e

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}  PROCAR Analysis Quick Start${NC}"
echo -e "${BLUE}  MoSSe/WSSe Janus Heterostructure${NC}"
echo -e "${BLUE}========================================${NC}"
echo ""

# Check if data directory is provided
if [ $# -eq 0 ]; then
    echo -e "${RED}Error: No data directory specified${NC}"
    echo ""
    echo "Usage: $0 <data-directory>"
    echo ""
    echo "Example:"
    echo "  $0 ../MoSSe_Se_up_WSSe_S_up"
    echo ""
    echo "The data directory should contain:"
    echo "  - PROCAR"
    echo "  - DOSCAR"
    echo "  - POSCAR or CONTCAR"
    exit 1
fi

DATA_DIR=$1

# Check if data directory exists
if [ ! -d "$DATA_DIR" ]; then
    echo -e "${RED}Error: Directory '$DATA_DIR' not found${NC}"
    exit 1
fi

echo -e "${YELLOW}Data directory: $DATA_DIR${NC}"
echo ""

# Check for required files
echo -e "${BLUE}Checking for required files...${NC}"

REQUIRED_FILES=("PROCAR" "DOSCAR")
STRUCTURE_FILES=("CONTCAR" "POSCAR")

MISSING=0

for file in "${REQUIRED_FILES[@]}"; do
    if [ -f "$DATA_DIR/$file" ]; then
        echo -e "  ${GREEN}✓${NC} Found $file"
    else
        echo -e "  ${RED}✗${NC} Missing $file"
        MISSING=1
    fi
done

# Check for structure file (either CONTCAR or POSCAR)
STRUCTURE_FOUND=0
for file in "${STRUCTURE_FILES[@]}"; do
    if [ -f "$DATA_DIR/$file" ]; then
        echo -e "  ${GREEN}✓${NC} Found $file"
        STRUCTURE_FOUND=1
        break
    fi
done

if [ $STRUCTURE_FOUND -eq 0 ]; then
    echo -e "  ${RED}✗${NC} Missing POSCAR or CONTCAR"
    MISSING=1
fi

if [ $MISSING -eq 1 ]; then
    echo ""
    echo -e "${RED}Error: Required files missing${NC}"
    echo ""
    echo "Please ensure the following files are in $DATA_DIR:"
    echo "  - PROCAR (orbital projections)"
    echo "  - DOSCAR (density of states)"
    echo "  - CONTCAR or POSCAR (structure)"
    echo ""
    echo "Note: PROCAR is typically excluded from git (.gitignore)"
    echo "      You need to copy it from your HPC calculation results"
    exit 1
fi

echo ""
echo -e "${GREEN}All required files found!${NC}"
echo ""

# Create output directory
OUTPUT_DIR="$DATA_DIR/procar_analysis"
echo -e "${BLUE}Creating output directory: $OUTPUT_DIR${NC}"
mkdir -p "$OUTPUT_DIR"
echo ""

# Run comprehensive analysis
echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}Running Comprehensive PROCAR Analysis${NC}"
echo -e "${BLUE}========================================${NC}"
echo ""

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

python3 "$SCRIPT_DIR/comprehensive_procar_analysis.py" \
    --data-dir "$DATA_DIR" \
    --output-dir "$OUTPUT_DIR"

# Check if analysis completed successfully
if [ $? -eq 0 ]; then
    echo ""
    echo -e "${GREEN}========================================${NC}"
    echo -e "${GREEN}  Analysis Complete!${NC}"
    echo -e "${GREEN}========================================${NC}"
    echo ""
    echo -e "${YELLOW}Results saved to: $OUTPUT_DIR${NC}"
    echo ""
    echo "Generated files:"
    echo "  - COMPREHENSIVE_ANALYSIS_SUMMARY.txt (start here!)"
    echo "  - band_character_analysis.txt"
    echo "  - valley_analysis_report.txt"
    echo "  - Various PNG plots"
    echo ""
    echo "Next steps:"
    echo "  1. Read COMPREHENSIVE_ANALYSIS_SUMMARY.txt for overview"
    echo "  2. Check band_character_analysis.txt for Type-I/Type-II alignment"
    echo "  3. Check valley_analysis_report.txt for valleytronics potential"
    echo "  4. View PNG plots for visualization"
    echo ""

    # Open summary file if possible
    if command -v cat &> /dev/null; then
        echo -e "${BLUE}Summary Report:${NC}"
        echo -e "${BLUE}----------------------------------------${NC}"
        cat "$OUTPUT_DIR/COMPREHENSIVE_ANALYSIS_SUMMARY.txt" | head -n 50
        echo ""
        echo -e "${YELLOW}(Showing first 50 lines, see full report in output directory)${NC}"
    fi
else
    echo ""
    echo -e "${RED}========================================${NC}"
    echo -e "${RED}  Analysis Failed${NC}"
    echo -e "${RED}========================================${NC}"
    echo ""
    echo "Please check the error messages above."
    echo ""
    echo "Common issues:"
    echo "  - Missing Python packages (numpy, matplotlib)"
    echo "  - Incorrect file format"
    echo "  - Corrupted VASP output files"
    exit 1
fi
