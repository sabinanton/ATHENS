#!/bin/bash

# Script to run Plotter.py with specified arguments

# Check for correct number of arguments
if [ "$#" -ne 1 ]; then
    echo "Usage: ./postprocess.sh <folder_name>"
    exit 1
fi

# Assign arguments to variables for clarity
FOLDER_NAME="$1"

# Path to the Python script. Update this if the script is in a different location.
PLOTTER_SCRIPT="./src/PostProcess.py"

# Check if the Plotter.py script exists
if [ ! -f "$PLOTTER_SCRIPT" ]; then
    echo "Error: $PLOTTER_SCRIPT does not exist."
    exit 1
fi

# Run the Python script with the arguments
python3 "$PLOTTER_SCRIPT" "$FOLDER_NAME"
