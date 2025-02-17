#!/bin/bash


# Activate Conda env
conda env create -f krasia.yml # This need fixing but no time now
conda activate krasia
# Function to parse YAML
parse_yaml() {
    local yaml_file="$1"
    sed -n -e "s|^\s*\([^:#]\+\)\s*:\s*\(.*\)\s*$|\1=\"\2\"|p" "$yaml_file"
}

# Check if YAML file is provided
if [ $# -ne 1 ]; then
    echo "Usage: $0 <config.yaml>"
    exit 1
fi

# Parse YAML file
yaml_file="$1"
eval $(parse_yaml "$yaml_file")

# Validate required paths
if [ -z "$ml_scripts_dir" ] || [ -z "$output_dir" ] || [ -z "$vcf_input" ] || [ -z "$ml_algorithm" ]; then
    echo "Error: Missing required fields in the YAML file."
    exit 1
fi

# Define transform_gt.sh path (always in ML directory)
transform_gt_script="$ml_scripts_dir/transform_gt.sh"

# Ensure required directories exist
mkdir -p "$output_dir/vcf_files"
mkdir -p "$output_dir/transformed_data"

# **Compress and Index Input VCF**
if [[ "$vcf_input" != *.gz ]]; then
    echo "Compressing VCF file: $vcf_input"
    bgzip -c "$vcf_input" > "$vcf_input.gz"
    vcf_input="$vcf_input.gz"
fi

echo "Indexing VCF file: $vcf_input"
bcftools index "$vcf_input"

# Save original sample names before merging
original_samples_file="$output_dir/vcf_files/original_samples.txt"
bcftools query -l "$vcf_input" > "$original_samples_file"

# Handle imputation if enabled
if [ "$enable_imputation" == "true" ]; then
    if [ -z "$beagle_dir" ] || [ -z "$wanted_positions_panel" ]; then
        echo "Error: beagle_dir and wanted_positions_panel are required when imputation is enabled."
        exit 1
    fi

    # Compress and index wanted positions panel
    if [[ "$wanted_positions_panel" != *.gz ]]; then
        echo "Compressing wanted positions panel..."
        bgzip -c "$wanted_positions_panel" > "$wanted_positions_panel.gz"
        wanted_positions_panel="$wanted_positions_panel.gz"
    fi

    echo "Indexing wanted positions panel..."
    bcftools index "$wanted_positions_panel"

    echo "Merging VCF with wanted positions panel..."
    merged_vcf="$output_dir/vcf_files/merged.vcf.gz"
    bcftools merge "$vcf_input" "$wanted_positions_panel" -Oz -o "$merged_vcf" --force-samples
    bcftools index "$merged_vcf"

    echo "Running Beagle for imputation..."
    imputed_vcf="$output_dir/vcf_files/imputed"
    java -jar "$beagle_dir" gt="$merged_vcf" out="$imputed_vcf"
    bcftools index "${imputed_vcf}.vcf.gz"

    vcf_input="${imputed_vcf}.vcf.gz"
fi

# **Filter out only original samples from imputed VCF**
filtered_vcf="$output_dir/vcf_files/original_samples_filtered.vcf.gz"
bcftools view -S "$original_samples_file" "$vcf_input" -Oz -o "$filtered_vcf"
bcftools index "$filtered_vcf"

# **Run ML scripts based on selected algorithm**
if [ "$ml_algorithm" == "algorithm1" ] || [ "$ml_algorithm" == "both" ]; then
    if [ -z "$bed_file_algorithm1" ]; then
        echo "Error: bed_file_algorithm1 is required for algorithm1."
        exit 1
    fi

    # Filter VCF for algorithm1
    filtered_vcf_algo1="$output_dir/vcf_files/final_filtered_algo1.vcf.gz"
    bcftools view -R "$bed_file_algorithm1" "$filtered_vcf" -Oz -o "$filtered_vcf_algo1"
    bcftools index "$filtered_vcf_algo1"

    # Transform genotype data for algorithm1
    transformed_file_algo1="$output_dir/transformed_data/transformed_algo1.txt"
    bash "$transform_gt_script" "$filtered_vcf_algo1" "$transformed_file_algo1"

    # Run ML script for algorithm1
    ml_script1="$ml_scripts_dir/ml_script1.py"
    python3 "$ml_script1" "$transformed_file_algo1" "$output_dir/ml_predictions_algo1.txt"
fi

if [ "$ml_algorithm" == "algorithm2" ] || [ "$ml_algorithm" == "both" ]; then
    if [ -z "$bed_file_algorithm2" ]; then
        echo "Error: bed_file_algorithm2 is required for algorithm2."
        exit 1
    fi

    # Filter VCF for algorithm2
    filtered_vcf_algo2="$output_dir/vcf_files/final_filtered_algo2.vcf.gz"
    bcftools view -R "$bed_file_algorithm2" "$filtered_vcf" -Oz -o "$filtered_vcf_algo2"
    bcftools index "$filtered_vcf_algo2"

    # Transform genotype data for algorithm2
    transformed_file_algo2="$output_dir/transformed_data/transformed_algo2.txt"
    bash "$transform_gt_script" "$filtered_vcf_algo2" "$transformed_file_algo2"

    # Run ML script for algorithm2
    ml_script2="$ml_scripts_dir/ml_script2.py"
    python3 "$ml_scripts_dir/ml_script2.py" "$ml_scripts_dir" "$transformed_file_algo2" "$bed_file_algorithm2" "$output_dir/ml_predictions_algo2.txt"

fi

echo "Processing completed successfully."
