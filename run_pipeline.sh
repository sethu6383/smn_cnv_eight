#!/bin/bash

# run_pipeline.sh - Enhanced SMN CNV detection pipeline focusing on exon 8
# Usage: ./run_pipeline.sh <input_bam_dir> [OPTIONS]

set -euo pipefail

# Default paths - UPDATED for exon 8 focus
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="$(dirname "$SCRIPT_DIR")"
CONFIG_DIR="/data/SMN/cnv_pipeline/config"
RESULTS_DIR="$PIPELINE_DIR/data/SMN/results_CNV"
BIN_DIR="/data/SMN/cnv_pipeline/bin"
LOG_DIR="$RESULTS_DIR/logs"

# Configuration files - UPDATED for exon 8 only
BED_FILE="/data/SMN/cnv_pipeline/config/smn_exons.bed"
SNP_FILE="/data/SMN/cnv_pipeline/config/discriminating_snps.txt"

# Pipeline options - ENHANCED
SKIP_PLOTS=false
VERBOSE=false
SAMPLE_TYPE="auto"
INPUT_BAM_DIR="/data/SMN/BAM"
ENABLE_ML=false
BOOTSTRAP_SAMPLES=1000
CONSOLIDATED_REPORT=true

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
PURPLE='\033[0;35m'
NC='\033[0m' # No Color

# Function to print colored output
print_status() {
    local color=$1
    local message=$2
    echo -e "${color}[$(date '+%Y-%m-%d %H:%M:%S')] ${message}${NC}"
}

print_error() {
    print_status "$RED" "ERROR: $1"
}

print_warning() {
    print_status "$YELLOW" "WARNING: $1"
}

print_info() {
    print_status "$BLUE" "INFO: $1"
}

print_success() {
    print_status "$GREEN" "SUCCESS: $1"
}

print_ml_info() {
    print_status "$PURPLE" "ML: $1"
}

# Enhanced dependency checking
check_dependencies() {
    print_info "Checking dependencies for enhanced pipeline..."
    
    local missing_tools=()
    local missing_python_packages=()
    
    # Check required command-line tools
    for tool in samtools python3; do
        if ! command -v "$tool" &> /dev/null; then
            missing_tools+=("$tool")
        fi
    done
    
    # Check Python packages - ENHANCED for ML
    python3 -c "
import sys
missing = []
required_packages = [
    'pandas', 'numpy', 'matplotlib', 'seaborn', 'scipy',
    'sklearn', 'joblib'  # Added for ML functionality
]

for package in required_packages:
    try:
        __import__(package)
    except ImportError:
        missing.append(package)

if missing:
    print('Missing Python packages:', ', '.join(missing))
    sys.exit(1)
else:
    print('All Python packages found')
" 2>/dev/null || missing_python_packages+=("ML packages (sklearn, joblib)")
    
    if [ ${#missing_tools[@]} -ne 0 ]; then
        print_error "Missing required tools: ${missing_tools[*]}"
        exit 1
    fi
    
    if [ ${#missing_python_packages[@]} -ne 0 ]; then
        print_warning "Missing packages: ${missing_python_packages[*]}"
        print_warning "ML features will be disabled. Install with: pip install scikit-learn joblib"
        ENABLE_ML=false
    fi
    
    print_success "Dependencies validated"
}

# Enhanced configuration validation
validate_config() {
    print_info "Validating enhanced configuration for exon 8 analysis..."
    
    # Check configuration files
    for file in "$BED_FILE" "$SNP_FILE"; do
        if [ ! -f "$file" ]; then
            print_error "Configuration file not found: $file"
            exit 1
        fi
    done
    
    # Validate BED file contains only exon 8
    if ! grep -q "exon8" "$BED_FILE"; then
        print_error "BED file does not contain exon 8 entries"
        exit 1
    fi
    
    # Check for exon 7 (should not be present in updated pipeline)
    if grep -q "exon7" "$BED_FILE"; then
        print_warning "BED file contains exon 7 entries - pipeline optimized for exon 8 only"
    fi
    
    # Validate input BAM directory
    if [ ! -d "$INPUT_BAM_DIR" ]; then
        print_error "Input BAM directory not found: $INPUT_BAM_DIR"
        exit 1
    fi
    
    # Check for BAM files
    local bam_count=$(find "$INPUT_BAM_DIR" -name "*.bam" | wc -l)
    if [ "$bam_count" -eq 0 ]; then
        print_error "No BAM files found in directory: $INPUT_BAM_DIR"
        exit 1
    fi
    
    # Validate minimum sample requirements for ML
    if [ "$ENABLE_ML" = true ] && [ "$bam_count" -lt 10 ]; then
        print_warning "Fewer than 10 samples - ML features may be unreliable"
    fi
    
    print_success "Enhanced configuration validated ($bam_count BAM files found)"
    print_info "Pipeline optimized for SMN1/SMN2 exon 8 analysis"
}

# Enhanced directory setup
setup_directories() {
    print_info "Setting up enhanced output directories..."
    
    mkdir -p "$RESULTS_DIR"/{depth,allele_counts,normalized,cnv_calls,reports,logs}
    mkdir -p "$RESULTS_DIR/normalized"/{plots,models}
    mkdir -p "$RESULTS_DIR/cnv_calls"/plots
    
    # Create analysis metadata
    cat > "$RESULTS_DIR/analysis_metadata.txt" << EOF
SMN CNV Analysis Metadata
========================
Pipeline Version: 2.0 Enhanced (Exon 8 Focus)
Analysis Date: $(date)
Input Directory: $INPUT_BAM_DIR
Sample Type: $SAMPLE_TYPE
ML Enabled: $ENABLE_ML
Bootstrap Samples: $BOOTSTRAP_SAMPLES
Configuration: Exon 8 only (SMN1/SMN2)
EOF
    
    print_success "Enhanced output directories created"
}

# Enhanced depth extraction (unchanged but with better logging)
run_depth_extraction() {
    print_info "Step 1: Extracting read depth per exon (EXON 8 FOCUS)..."
    
    local output_dir="$RESULTS_DIR/depth"
    local log_file="$LOG_DIR/depth_extraction.log"
    
    if ! bash "$BIN_DIR/extract_depth.sh" "$INPUT_BAM_DIR" "$BED_FILE" "$output_dir" "$SAMPLE_TYPE" 2>&1 | tee "$log_file"; then
        print_error "Depth extraction failed. Check log: $log_file"
        exit 1
    fi
    
    local depth_files=$(find "$output_dir" -name "*_depth.txt" | wc -l)
    if [ "$depth_files" -eq 0 ]; then
        print_error "No depth files were created"
        exit 1
    fi
    
    print_success "Depth extraction completed ($depth_files files created)"
    print_info "Focus: SMN1_exon8 and SMN2_exon8 regions only"
}

# Enhanced coverage calculation
run_coverage_calculation() {
    print_info "Step 2: Calculating average coverage per exon (ENHANCED)..."
    
    local input_dir="$RESULTS_DIR/depth"
    local output_file="$RESULTS_DIR/depth/coverage_summary.txt"
    local log_file="$LOG_DIR/coverage_calculation.log"
    
    if ! python3 "$BIN_DIR/calculate_coverage.py" "$input_dir" "$BED_FILE" "$output_file" 2>&1 | tee "$log_file"; then
        print_error "Coverage calculation failed. Check log: $log_file"
        exit 1
    fi
    
    if [ ! -f "$output_file" ]; then
        print_error "Coverage summary file was not created"
        exit 1
    fi
    
    # Enhanced coverage QC
    python3 -c "
import pandas as pd
df = pd.read_csv('$output_file', sep='\t')
print(f'Coverage summary: {len(df)} sample-exon combinations')
print(f'Mean coverage: {df[\"avg_coverage\"].mean():.1f}x')
print(f'Samples with low coverage (<20x): {(df[\"avg_coverage\"] < 20).sum()}')
"
    
    print_success "Enhanced coverage calculation completed"
}

# Enhanced allele counting with exon 8 focus
run_allele_counting() {
    print_info "Step 3: Performing allele-specific counting (EXON 8 SNPS)..."
    
    local output_dir="$RESULTS_DIR/allele_counts"
    local log_file="$LOG_DIR/allele_counting.log"
    
    local cmd="python3 $BIN_DIR/generate_report.py $cn_file $allele_file $output_dir"
    if [ "$SKIP_PLOTS" = true ]; then
        cmd="$cmd --format html"
    fi
    
    if ! eval "$cmd" 2>&1 | tee "$log_file"; then
        print_error "Enhanced report generation failed. Check log: $log_file"
        exit 1
    fi
    
    local report_count=$(find "$output_dir" -name "*_report.html" | wc -l)
    print_success "Enhanced report generation completed ($report_count reports created)"
}

# NEW: Generate consolidated MultiQC-style report
run_consolidated_report() {
    if [ "$CONSOLIDATED_REPORT" = false ]; then
        return 0
    fi
    
    print_info "Step 7: Generating consolidated MultiQC-style report..."
    
    local output_prefix="$RESULTS_DIR/SMN_CNV_Analysis"
    local log_file="$LOG_DIR/consolidated_report.log"
    
    if ! python3 "$BIN_DIR/generate_consolidated_report.py" "$RESULTS_DIR" "$output_prefix" --format html,txt 2>&1 | tee "$log_file"; then
        print_warning "Consolidated report generation failed. Check log: $log_file"
        return 1
    fi
    
    print_success "Consolidated MultiQC-style report generated"
    print_info "HTML report: ${output_prefix}_consolidated_report.html"
    print_info "Text summary: ${output_prefix}_consolidated_report.txt"
}

# Enhanced pipeline summary
create_enhanced_summary() {
    print_info "Creating enhanced pipeline summary..."
    
    local summary_file="$RESULTS_DIR/pipeline_summary.txt"
    local sample_count=$(find "$INPUT_BAM_DIR" -name "*.bam" | wc -l)
    
    cat > "$summary_file" << EOF
SMN CNV Detection Pipeline - Enhanced Summary (v2.0)
===================================================

PIPELINE CONFIGURATION:
- Version: 2.0 Enhanced (Exon 8 Focus)
- Analysis Date: $(date)
- Pipeline Directory: $PIPELINE_DIR
- Configuration Directory: $CONFIG_DIR
- Results Directory: $RESULTS_DIR
- Input BAM Directory: $INPUT_BAM_DIR

TARGET REGIONS:
- Primary Focus: SMN1/SMN2 Exon 8 Only
- BED File: $BED_FILE
- SNP Configuration: $SNP_FILE
- Rationale: Consistent lack of reads in exon 7 across datasets

SAMPLE INFORMATION:
- Total BAM Files: $sample_count
- Sample Type Strategy: $SAMPLE_TYPE
- Reference Samples Required: ≥3 (recommended ≥5)

ENHANCED FEATURES:
- ML-Enhanced Analysis: $ENABLE_ML
- Bootstrap Confidence: $BOOTSTRAP_SAMPLES samples
- Robust Outlier Detection: Multiple methods (IQR, Isolation Forest, Modified Z-score)
- Probabilistic CNV Calling: Consensus-based with uncertainty quantification
- Advanced Normalization: Robust scaling with MAD and trimmed statistics

ANALYSIS THRESHOLDS (Z-score based):
- Homozygous Deletion (CN=0): ≤-2.5
- Heterozygous Deletion (CN=1): -2.5 to -1.5  
- Normal Copy Number (CN=2): -1.5 to +1.5
- Duplication (CN=3): +1.5 to +2.5
- High Amplification (CN=4+): >+2.5

CLINICAL INTERPRETATION:
- SMN1 CN=0: Likely SMA affected (homozygous deletion)
- SMN1 CN=1: SMA carrier (heterozygous deletion)  
- SMN1 CN=2: Normal copy number
- SMN1 CN≥3: Gene duplication (may be protective)

QUALITY CONTROL RECOMMENDATIONS:
- Minimum consensus score: 0.7
- Minimum confidence level: medium or high
- Maximum outlier probability: 0.1
- Minimum bootstrap confidence: 0.5 (if enabled)
- Coverage requirement: ≥20x across exon 8

OUTPUT STRUCTURE:
- Depth Files: $RESULTS_DIR/depth/
- Coverage Summary: $RESULTS_DIR/depth/coverage_summary.txt
- Allele Counts: $RESULTS_DIR/allele_counts/allele_counts.txt
- Z-scores: $RESULTS_DIR/normalized/z_scores.txt
- Copy Numbers: $RESULTS_DIR/cnv_calls/copy_numbers.txt
- Gene-level Results: $RESULTS_DIR/cnv_calls/copy_numbers_gene_level.txt
- Individual Reports: $RESULTS_DIR/reports/
- Consolidated Report: $RESULTS_DIR/SMN_CNV_Analysis_consolidated_report.html
- Log Files: $RESULTS_DIR/logs/

IMPORTANT NOTES:
- This pipeline focuses exclusively on exon 8 for improved accuracy
- All CNV calls should be validated through orthogonal methods
- Clinical interpretation requires correlation with patient phenotype
- Consider family studies for carrier samples when appropriate
- ML models are trained on reference samples and may need retraining for different populations

EOF

    # Add runtime statistics if available
    if [ -f "$LOG_DIR/depth_extraction.log" ]; then
        echo "" >> "$summary_file"
        echo "RUNTIME STATISTICS:" >> "$summary_file"
        local start_time=$(head -1 "$LOG_DIR/depth_extraction.log" | grep -o '\[.*\]' | head -1)
        local end_time=$(date '+[%Y-%m-%d %H:%M:%S]')
        echo "- Start Time: $start_time" >> "$summary_file"
        echo "- End Time: $end_time" >> "$summary_file"
    fi
    
    print_success "Enhanced pipeline summary created: $summary_file"
}

# Enhanced usage function
show_enhanced_usage() {
    cat << EOF
SMN CNV Detection Pipeline - Enhanced v2.0 (Exon 8 Focus)

Usage: $0 <input_bam_dir> [OPTIONS]

REQUIRED:
    input_bam_dir       Directory containing BAM files to analyze

OPTIONS:
    --config DIR        Configuration directory (default: $CONFIG_DIR)
    --results DIR       Results directory (default: $RESULTS_DIR)
    --sample-type TYPE  Sample type: reference, test, or auto (default: auto)
    --enable-ml         Enable machine learning enhanced analysis
    --bootstrap N       Number of bootstrap samples (default: 1000, 0 to disable)
    --skip-plots        Skip generating plots to speed up analysis
    --no-consolidated   Skip consolidated MultiQC-style report
    --verbose           Enable verbose output
    --help              Show this help message

ENHANCED FEATURES (v2.0):
    ✨ Exon 8 Focus: Optimized for SMN1/SMN2 exon 8 only (improved accuracy)
    🤖 ML Enhancement: Random Forest, Gaussian Mixture Models, Isolation Forest
    📊 Bootstrap Confidence: Statistical confidence intervals for CNV calls
    🎯 Robust Statistics: Multiple outlier detection methods, trimmed means
    🔬 Probabilistic Calling: Consensus-based CNV calling with uncertainty
    📋 Consolidated Reports: MultiQC-style HTML and text summary reports
    🔍 Advanced QC: Enhanced quality control metrics and recommendations

SAMPLE TYPE AUTO-DETECTION:
    When --sample-type is 'auto' (default):
    - Files with 'ref', 'control', or 'normal' → reference samples
    - All other files → test samples
    - Override with --sample-type reference/test for all samples

MACHINE LEARNING FEATURES:
    --enable-ml activates:
    - Gaussian Mixture Models for copy number clustering
    - Random Forest classifiers for enhanced prediction
    - Isolation Forest for anomaly detection
    - Feature importance analysis
    - Cross-validation scoring

BOOTSTRAP ANALYSIS:
    --bootstrap N enables:
    - Statistical confidence intervals
    - Threshold optimization using reference samples
    - Uncertainty quantification for individual calls
    - Improved reliability assessment

QUALITY REQUIREMENTS:
    - Minimum 3 reference samples (5+ recommended)
    - Coverage ≥20x across exon 8 regions
    - BAM files must be indexed (.bai files)
    - Python packages: pandas, numpy, matplotlib, seaborn, scipy
    - ML packages (optional): scikit-learn, joblib

EXAMPLES:
    # Basic enhanced analysis
    $0 /path/to/bam/files/
    
    # Full ML-enhanced analysis
    $0 /path/to/bam/files/ --enable-ml --bootstrap 2000
    
    # All samples are references (for building reference database)
    $0 /path/to/bam/files/ --sample-type reference --enable-ml
    
    # Fast analysis without plots
    $0 /path/to/bam/files/ --skip-plots --bootstrap 0
    
    # Custom output location
    $0 /path/to/bam/files/ --results /custom/output/ --enable-ml

OUTPUT HIGHLIGHTS:
    📊 Consolidated Report: SMN_CNV_Analysis_consolidated_report.html
    📋 Text Summary: SMN_CNV_Analysis_consolidated_report.txt
    📈 Individual Reports: reports/[SAMPLE_ID]/[SAMPLE_ID]_report.html
    🎯 Gene-level Results: cnv_calls/copy_numbers_gene_level.txt

CLINICAL NOTES:
    - SMN1 CN≤0.5: Potential SMA affected - requires immediate attention
    - SMN1 CN=1: SMA carrier - consider genetic counseling  
    - SMN2 CN variations may modify SMA severity
    - All positive results should be confirmed by orthogonal methods
EOF
}

# Parse enhanced command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --help)
            show_enhanced_usage
            exit 0
            ;;
        --config)
            CONFIG_DIR="$2"
            BED_FILE="$CONFIG_DIR/smn_exons.bed"
            SNP_FILE="$CONFIG_DIR/discriminating_snps.txt"
            shift 2
            ;;
        --results)
            RESULTS_DIR="$2"
            LOG_DIR="$RESULTS_DIR/logs"
            shift 2
            ;;
        --sample-type)
            SAMPLE_TYPE="$2"
            if [[ ! "$SAMPLE_TYPE" =~ ^(reference|test|auto)$ ]]; then
                print_error "Invalid sample type: $SAMPLE_TYPE. Must be 'reference', 'test', or 'auto'"
                exit 1
            fi
            shift 2
            ;;
        --enable-ml)
            ENABLE_ML=true
            print_ml_info "Machine Learning features enabled"
            shift
            ;;
        --bootstrap)
            BOOTSTRAP_SAMPLES="$2"
            if ! [[ "$BOOTSTRAP_SAMPLES" =~ ^[0-9]+$ ]]; then
                print_error "Bootstrap samples must be a positive integer"
                exit 1
            fi
            shift 2
            ;;
        --skip-plots)
            SKIP_PLOTS=true
            shift
            ;;
        --no-consolidated)
            CONSOLIDATED_REPORT=false
            shift
            ;;
        --verbose)
            VERBOSE=true
            set -x
            shift
            ;;
        -*)
            print_error "Unknown option: $1"
            show_enhanced_usage
            exit 1
            ;;
        *)
            if [ -z "$INPUT_BAM_DIR" ]; then
                INPUT_BAM_DIR="$1"
            else
                print_error "Multiple input directories specified"
                show_enhanced_usage
                exit 1
            fi
            shift
            ;;
    esac
done

# Validate input directory
if [ -z "$INPUT_BAM_DIR" ]; then
    print_error "Input BAM directory is required"
    show_enhanced_usage
    exit 1
fi

# Main enhanced pipeline execution
main() {
    print_info "Starting Enhanced SMN CNV Detection Pipeline v2.0"
    print_info "🎯 FOCUS: Exon 8 analysis for improved accuracy"
    print_info "Input BAM directory: $INPUT_BAM_DIR"
    print_info "Configuration directory: $CONFIG_DIR"  
    print_info "Results directory: $RESULTS_DIR"
    print_info "Sample type strategy: $SAMPLE_TYPE"
    
    if [ "$ENABLE_ML" = true ]; then
        print_ml_info "Machine Learning enhancement: ENABLED"
    fi
    
    if [ "$BOOTSTRAP_SAMPLES" -gt 0 ]; then
        print_info "Bootstrap confidence analysis: $BOOTSTRAP_SAMPLES samples"
    fi
    
    # Pre-flight checks
    check_dependencies
    validate_config
    setup_directories
    
    # Execute enhanced pipeline steps
    local start_time=$(date +%s)
    
    run_depth_extraction
    run_coverage_calculation  
    run_allele_counting
    run_enhanced_normalization
    run_enhanced_cnv_estimation
    run_enhanced_report_generation
    run_consolidated_report
    
    # Create enhanced summary
    create_enhanced_summary
    
    local end_time=$(date +%s)
    local duration=$((end_time - start_time))
    
    print_success "Enhanced SMN CNV Pipeline completed successfully!"
    print_info "Total runtime: ${duration} seconds"
    print_info "Results available in: $RESULTS_DIR"
    
    # Show enhanced results summary
    if [ -f "$RESULTS_DIR/cnv_calls/copy_numbers_gene_level.txt" ]; then
        print_info "📊 Enhanced Results Summary:"
        python3 -c "
import pandas as pd
try:
    # Gene-level results
    df = pd.read_csv('$RESULTS_DIR/cnv_calls/copy_numbers_gene_level.txt', sep='\t')
    samples = df['sample_id'].unique()
    print(f'  📋 Analyzed samples: {len(samples)}')
    
    # SMN1 clinical significance
    smn1_data = df[df['gene'] == 'SMN1']
    if not smn1_data.empty and 'clinical_significance' in smn1_data.columns:
        clin_summary = smn1_data['clinical_significance'].value_counts()
        print('  🧬 SMN1 Clinical Findings:')
        for category, count in clin_summary.items():
            emoji = '🚨' if 'affected' in category.lower() else '⚠️' if 'carrier' in category.lower() else '✅'
            print(f'    {emoji} {category}: {count} samples')
    
    # Copy number distribution
    if not smn1_data.empty:
        cn_dist = smn1_data['estimated_copy_number'].round().value_counts().sort_index()
        print('  📈 SMN1 Copy Number Distribution:')
        for cn, count in cn_dist.items():
            print(f'    CN={cn}: {count} samples')
    
    # Quality metrics
    if 'confidence_score' in df.columns:
        high_conf = (df['confidence_score'] > 0.8).sum()
        total = len(df)
        print(f'  🎯 High confidence calls: {high_conf}/{total} ({high_conf/total*100:.1f}%)')
    
except Exception as e:
    print(f'  ⚠️ Could not generate summary: {e}')
"
    fi
    
    # Enhanced output guidance
    echo ""
    print_info "🗂️ Key Output Files:"
    print_info "  📊 Consolidated Report: $RESULTS_DIR/SMN_CNV_Analysis_consolidated_report.html"
    print_info "  📋 Text Summary: $RESULTS_DIR/SMN_CNV_Analysis_consolidated_report.txt"
    print_info "  🎯 Gene Results: $RESULTS_DIR/cnv_calls/copy_numbers_gene_level.txt"
    print_info "  📈 Individual Reports: $RESULTS_DIR/reports/"
    
    echo ""
    print_info "🔍 Quality Control Recommendations:"
    print_info "  • Review consolidated report for overall cohort analysis"
    print_info "  • Validate any SMN1 CN≤1 samples with orthogonal methods"
    print_info "  • Check confidence scores and bootstrap intervals"
    print_info "  • Consider clinical correlation for carrier samples"
    
    if [ "$ENABLE_ML" = true ]; then
        print_ml_info "🤖 ML models saved in: $RESULTS_DIR/normalized/models/"
    fi
}

# Run enhanced main function
main "$@"d="python3 $BIN_DIR/allele_count.py $INPUT_BAM_DIR $SNP_FILE $output_dir"
    if [ "$SAMPLE_TYPE" != "auto" ]; then
        cmd="$cmd --sample-type $SAMPLE_TYPE"
    fi
    
    if ! eval "$cmd" 2>&1 | tee "$log_file"; then
        print_error "Allele counting failed. Check log: $log_file"
        exit 1
    fi
    
    local allele_file="$output_dir/allele_counts.txt"
    if [ ! -f "$allele_file" ]; then
        print_error "Allele counts file was not created"
        exit 1
    fi
    
    print_success "Enhanced allele counting completed (exon 8 discriminating SNPs)"
}

# Enhanced normalization with ML features
run_enhanced_normalization() {
    print_info "Step 4: Enhanced normalization with ML features..."
    
    local coverage_file="$RESULTS_DIR/depth/coverage_summary.txt"
    local sample_info_file="$RESULTS_DIR/allele_counts/sample_info.txt"
    local allele_file="$RESULTS_DIR/allele_counts/allele_counts.txt"
    local output_file="$RESULTS_DIR/normalized/z_scores.txt"
    local log_file="$LOG_DIR/normalization.log"
    
    local cmd="python3 $BIN_DIR/normalize_coverage.py $coverage_file $sample_info_file $output_file"
    
    if [ "$ENABLE_ML" = true ]; then
        cmd="$cmd --enable-ml --allele-file $allele_file"
        print_ml_info "ML features enabled for enhanced normalization"
    fi
    
    if ! eval "$cmd" 2>&1 | tee "$log_file"; then
        print_error "Enhanced normalization failed. Check log: $log_file"
        exit 1
    fi
    
    if [ ! -f "$output_file" ]; then
        print_error "Z-scores file was not created"
        exit 1
    fi
    
    print_success "Enhanced normalization completed with robust statistics"
}

# Enhanced CNV estimation with ML and bootstrap
run_enhanced_cnv_estimation() {
    print_info "Step 5: Enhanced CNV estimation with ML and bootstrap confidence..."
    
    local z_scores_file="$RESULTS_DIR/normalized/z_scores.txt"
    local output_file="$RESULTS_DIR/cnv_calls/copy_numbers.txt"
    local log_file="$LOG_DIR/copy_number_estimation.log"
    
    local cmd="python3 $BIN_DIR/estimate_copy_number.py $z_scores_file $output_file"
    cmd="$cmd --bootstrap-samples $BOOTSTRAP_SAMPLES"
    
    if [ "$ENABLE_ML" = true ]; then
        cmd="$cmd --enable-ml"
        print_ml_info "ML-enhanced CNV calling enabled"
    fi
    
    if [ "$SKIP_PLOTS" = true ]; then
        cmd="$cmd --no-plots"
    fi
    
    if ! eval "$cmd" 2>&1 | tee "$log_file"; then
        print_error "Enhanced CNV estimation failed. Check log: $log_file"
        exit 1
    fi
    
    if [ ! -f "$output_file" ]; then
        print_error "Copy numbers file was not created"
        exit 1
    fi
    
    print_success "Enhanced CNV estimation completed with bootstrap confidence intervals"
}

# Enhanced report generation
run_enhanced_report_generation() {
    print_info "Step 6: Generating enhanced per-sample reports..."
    
    local cn_file="$RESULTS_DIR/cnv_calls/copy_numbers.txt"
    local allele_file="$RESULTS_DIR/allele_counts/allele_counts.txt"
    local output_dir="$RESULTS_DIR/reports"
    local log_file="$LOG_DIR/report_generation.log"
    
    local cm
