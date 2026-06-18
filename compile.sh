#!/bin/bash

# Script to compile Run_Parallel.cpp using mpicxx

# Name of the source file
SOURCE_FILE="src/Run_parallel.cpp"

# Name of the output executable
EXECUTABLE="athens"

# Compilation options
COMPILE_OPTIONS=""

# Check if the source file exists
if [ ! -f "$SOURCE_FILE" ]; then
    echo "Error: Source file $SOURCE_FILE does not exist."
    exit 1
fi

# Command to compile the code
mpicxx $COMPILE_OPTIONS $SOURCE_FILE -o $EXECUTABLE

# Check if the compilation was successful
if [ $? -eq 0 ]; then
    echo "Compilation successful. Executable created: $EXECUTABLE"
else
    echo "Compilation failed."
    exit 1
fi
