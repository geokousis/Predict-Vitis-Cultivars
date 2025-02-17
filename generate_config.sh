#!/bin/bash

# Check if Zenity is installed
if ! command -v zenity &> /dev/null; then
    echo "Zenity is not installed. Please install it using: sudo apt install zenity"
    exit 1
fi

# Welcome message
zenity --info --title="Configuration Setup" --text="Welcome! This tool will help you create a config.yaml file for your pipeline."

# Ask for output directory
output_dir=$(zenity --file-selection --directory --title="Select Output Directory")
if [[ -z "$output_dir" ]]; then
    zenity --error --text="Output directory is required!"
    exit 1
fi

# Ask for ML scripts directory (transform_gt.sh is inside this directory)
ml_scripts_dir=$(zenity --file-selection --directory --title="Select ML Scripts Directory")
if [[ -z "$ml_scripts_dir" ]]; then
    zenity --error --text="ML scripts directory is required!"
    exit 1
fi

# Ask for VCF input file
vcf_input=$(zenity --file-selection --title="Select VCF Input File")
if [[ -z "$vcf_input" ]]; then
    zenity --error --text="VCF input file is required!"
    exit 1
fi

# Ask user whether to enable imputation
enable_imputation=$(zenity --list --title="Enable Imputation?" --radiolist --column="Select" --column="Option" FALSE "true" TRUE "false")
if [[ -z "$enable_imputation" ]]; then
    zenity --error --text="Please select an option for imputation."
    exit 1
fi

# If imputation is enabled, ask for Beagle JAR file and wanted positions panel
if [[ "$enable_imputation" == "true" ]]; then
    beagle_file=$(zenity --file-selection --title="Select Beagle JAR File (beagle.jar)")
    if [[ -z "$beagle_file" ]]; then
        zenity --error --text="Beagle JAR file is required!"
        exit 1
    fi

    wanted_positions_panel=$(zenity --file-selection --title="Select Wanted Positions Panel (VCF File)")
    if [[ -z "$wanted_positions_panel" ]]; then
        zenity --error --text="Wanted positions panel is required!"
        exit 1
    fi
fi

# Ask user to select ML algorithm
ml_algorithm=$(zenity --list --title="Select ML Algorithm" --radiolist --column="Select" --column="Algorithm" FALSE "algorithm1" FALSE "algorithm2" TRUE "both")
if [[ -z "$ml_algorithm" ]]; then
    zenity --error --text="Please select an ML algorithm."
    exit 1
fi

# Ask for BED files based on ML algorithm selection
if [[ "$ml_algorithm" == "algorithm1" || "$ml_algorithm" == "both" ]]; then
    bed_file_algorithm1=$(zenity --file-selection --title="Select BED File for Algorithm 1")
    if [[ -z "$bed_file_algorithm1" ]]; then
        zenity --error --text="BED file for Algorithm 1 is required!"
        exit 1
    fi
fi

if [[ "$ml_algorithm" == "algorithm2" || "$ml_algorithm" == "both" ]]; then
    bed_file_algorithm2=$(zenity --file-selection --title="Select BED File for Algorithm 2")
    if [[ -z "$bed_file_algorithm2" ]]; then
        zenity --error --text="BED file for Algorithm 2 is required!"
        exit 1
    fi
fi

# Create config.yaml file
config_file="config.yaml"
echo "Creating $config_file..."

cat > "$config_file" <<EOL
ml_scripts_dir: "$ml_scripts_dir"
output_dir: "$output_dir"
vcf_input: "$vcf_input"
enable_imputation: $enable_imputation
EOL

if [[ "$enable_imputation" == "true" ]]; then
cat >> "$config_file" <<EOL
beagle_file: "$beagle_file"
wanted_positions_panel: "$wanted_positions_panel"
EOL
fi

cat >> "$config_file" <<EOL
ml_algorithm: "$ml_algorithm"
EOL

if [[ "$ml_algorithm" == "algorithm1" || "$ml_algorithm" == "both" ]]; then
cat >> "$config_file" <<EOL
bed_file_algorithm1: "$bed_file_algorithm1"
EOL
fi

if [[ "$ml_algorithm" == "algorithm2" || "$ml_algorithm" == "both" ]]; then
cat >> "$config_file" <<EOL
bed_file_algorithm2: "$bed_file_algorithm2"
EOL
fi

zenity --info --title="Success!" --text="Configuration file '$config_file' has been created successfully!"

# Ask if the user wants to start the pipeline
run_pipeline=$(zenity --question --title="Run Pipeline?" --text="Do you want to start the pipeline now?" --ok-label="Yes" --cancel-label="No")

if [[ $? -eq 0 ]]; then
    if [[ ! -f "run.sh" ]]; then
        zenity --error --text="Error: run.sh not found! Please ensure it exists in the current directory."
        exit 1
    fi
    zenity --info --text="Starting the pipeline..."
    bash run.sh "$config_file"
else
    zenity --info --text="You can start the pipeline later by running: ./run.sh config.yaml"
fi

echo "Configuration complete."
