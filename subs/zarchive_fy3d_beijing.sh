#!/bin/bash

# Base directory containing the 2020, 2021, 2022 folders
source_dir="descend"
# Base directory where the files will be moved
target_base_dir="/home/jihenghu/fy03/FY3D/descend"

# Loop over the years
for year in 2020 2021 2022; do
    # Find all .HDF files in the current year's directory
    find "$source_dir/$year" -type f -name "*.HDF" | while read filepath; do
        # Extract the subdirectory structure after the year, e.g., '20200116/FY3D_MWRIA_GBAL_L1_20200116_1142_010KM_MS.HDF'
        # echo $filepath 
        
        sub_dir=$(echo "$filepath" | cut -d "_" -f5)
     
        # Construct the target directory
        target_dir=$target_base_dir/$year/$sub_dir
        # echo $target_dir
        # exit
        # Create the target directory if it doesn't exist
        mkdir -p "$target_dir"
        
        # Move the file
        mv -v "$filepath" "$target_dir/"
    done
done

echo "Files have been moved successfully!"
