#!/bin/bash

# Usage: transform_gt.sh input_vcf output_file
# This script extracts genotype data from a VCF file (gzipped or not),
# transforms the genotypes into numeric format, and writes the results.

input_vcf="$1"
output="$2"

if [[ -z "$input_vcf" || -z "$output" ]]; then
    echo "Usage: $0 input_vcf output_file"
    echo "Transforms genotype data without removing duplicates."
    exit 1
fi

echo "Processing $input_vcf..."

# Function to output VCF contents regardless of compression
query_cmd() {
    if [[ "$input_vcf" == *.gz ]]; then
        bcftools view "$input_vcf"
    else
        cat "$input_vcf"
    fi
}

# Extract data and transform genotypes without duplicate filtering
query_cmd | bcftools query -f '%CHROM\t%POS[\t%GT]\n' | \
awk 'BEGIN { FS=OFS="\t" }
function transform(gt) {
    if (gt == "0/0" || gt == "0|0") return "0"
    else if (gt == "1/1" || gt == "1|1") return "2"
    else if (gt == "0/1" || gt == "1/0" || gt == "0|1" || gt == "1|0") return "1"
    else return "NA"
}
{
    # Transform genotype fields
    for (i = 3; i <= NF; i++) {
        $i = transform($i)
    }
    # Print every line (including duplicates)
    print
}' > "$output"

# Display line counts
total_lines=$(query_cmd | bcftools query -f '%CHROM\t%POS[\t%GT]\n' | wc -l)
output_lines=$(wc -l < "$output")

echo "Total lines processed: $total_lines"
echo "Lines in output: $output_lines"
