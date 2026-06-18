#!/bin/bash

# Script to run Plotter.py with specified arguments

# Check for correct number of arguments
if [ "$#" -ne 5 ]; then
    echo "Usage: ./run_plotter.sh <folder_name> <is_batch> <latitude> <longitude> <zoom>"
    exit 1
fi

# Assign arguments to variables for clarity
FOLDER_NAME="$1"
BATCH="$2"  # This is expected to be "--batch", you can validate this if needed
LATITUDE="$3"
LONGITUDE="$4"
ZOOM="$5"

# Path to the Python script. Update this if the script is in a different location.
PLOTTER_SCRIPT="./src/Plotter.py"

# Check if the Plotter.py script exists
if [ ! -f "$PLOTTER_SCRIPT" ]; then
    echo "Error: $PLOTTER_SCRIPT does not exist."
    exit 1
fi

# Run the Python script with the arguments
python3 "$PLOTTER_SCRIPT" "$FOLDER_NAME" "$BATCH" "$LATITUDE" "$LONGITUDE" "$ZOOM"
