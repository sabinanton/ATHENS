#!/bin/bash

# Function to display help message if incorrect usage
usage() {
    echo "Usage: $0 <directory>"
    echo "Provide the directory from which subfolders will be removed if they exist."
    exit 1
}

# Check if the correct number of arguments is passed
if [ "$#" -ne 1 ]; then
    echo "Error: Incorrect number of arguments."
    usage
fi

# The directory from which to remove subfolders
TARGET_DIRECTORY=$1

# Check if the provided argument is a directory
if [ ! -d "$TARGET_DIRECTORY" ]; then
    echo "Error: Provided path ($TARGET_DIRECTORY) is not a valid directory."
    exit 1
fi

# Subfolders to remove
SUBFOLDERS=("Batch_Results" "Batch_Plots" "Single_Results" "Single_Plots")

echo "Starting the script for directory: $TARGET_DIRECTORY"

# Loop through each subfolder and remove it if it exists
for SUBFOLDER in "${SUBFOLDERS[@]}"
do
    FULL_PATH="$TARGET_DIRECTORY/$SUBFOLDER"
    if [ -d "$FULL_PATH" ]; then
        rm -rf "$FULL_PATH"
    fi
done

echo "Folder cleaned."
