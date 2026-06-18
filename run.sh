#!/bin/bash

# Function to display help message
usage() {
    echo ""
    echo "Usage: $0 -project=<name_of_project> -np=<number_of_processes> -id=<id_of_project> [-batch]"
    echo "  -project     Name of the project."
    echo "  -np          Number of processes to use for mpirun."
    echo "  -id          ID of the simulation run."
    echo "  -batch       Optional. If provided, runs the executable with the -batch option."
    exit 1
}

# Initialize variables
PROJECT_NAME=""
NUM_PROCESSES=""
SIM_ID=""
BATCH_MODE=false

# Check if no arguments were provided
if [ "$#" -eq 0 ]; then
    usage
fi

# Parse command-line arguments
for i in "$@"
do
case $i in
    -project=*)
    PROJECT_NAME="${i#*=}"
    shift # past argument=value
    ;;
    -np=*)
    NUM_PROCESSES="${i#*=}"
    shift # past argument=value
    ;;
    -id=*)
    SIM_ID="${i#*=}"
    shift # past argument=value
    ;;
    -batch)
    BATCH_MODE=true
    shift # past argument with no value
    ;;
    *)
    # unknown option
    usage
    ;;
esac
done

# Validate necessary arguments
if [ -z "$PROJECT_NAME" ] || [ -z "$NUM_PROCESSES" ] || [ -z "$SIM_ID" ]; then
    usage
fi

# Define the path to the executable.
# This example assumes that the executable is in the current directory.
EXECUTABLE="./athens"

# Check if the executable exists
if [ ! -f "$EXECUTABLE" ]; then
    echo "Error: Executable $EXECUTABLE does not exist."
    exit 1
fi

# Check if mpirun is installed
if ! command -v mpirun &> /dev/null
then
    echo "Error: mpirun could not be found. Please install MPI."
    exit 1
fi

# Construct the command with the necessary options
COMMAND="mpirun -np $NUM_PROCESSES $EXECUTABLE --project=$PROJECT_NAME --id=$SIM_ID"

if $BATCH_MODE ; then
  COMMAND+=" --batch"
fi

# Run the command
eval $COMMAND

# Check the exit status of the mpirun command and respond accordingly
if [ "$?" -ne 0 ]; then
    echo "Error: Execution with mpirun failed."
    exit 1
else
    echo "Parallel execution successful."
fi