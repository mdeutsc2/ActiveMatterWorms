#!/bin/bash

set -e

SIMULATION="../amatter.x"
ARGS=""
dirname="data_$(date +%d_%m_%Y)"

if [ -d "$dirname" ]; then
    count=2
    while [ -d "$dirname-$count" ]; do
        count=$((count+1))
    done
    dirname="$dirname-$count"
fi

if [ -z "$ARGS" ] && [ $# -gt 0 ]; then
    ARGS="$@"
fi

mkdir -p "$dirname"
cd "$dirname"
echo "Running with arguments: $ARGS"
time $SIMULATION $ARGS
