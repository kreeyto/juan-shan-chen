#!/bin/bash

if [ "$#" -lt 2 ]; then
    echo "Usage: ./run_post.sh <ID> <velocity_set>"
    exit 1
fi

ID="$1"
VELOCITY_SET="$2"

if command -v python3 &> /dev/null; then
    PYTHON_CMD="python3"
elif command -v python &> /dev/null; then
    PYTHON_CMD="python"
else
    echo "Python não encontrado no sistema."
    exit 1
fi

$PYTHON_CMD process_steps.py "$ID" "$VELOCITY_SET"
