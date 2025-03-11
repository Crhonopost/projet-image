#!/bin/bash

# Output file to collect all PSNR values
output_file="all_psnr_values.psnr"

# Clear the output file if it exists
> "$output_file"

# Function to process a single .psnr file
process_psnr_file() {
    psnr_file=$1
    psnr_value=$(head -n 1 "$psnr_file")
    echo "$psnr_file : $psnr_value" >> "$output_file"
}

# Find all .psnr files and process them
for psnr_file in $(find . -name '*.psnr'); do
    process_psnr_file "$psnr_file"
done