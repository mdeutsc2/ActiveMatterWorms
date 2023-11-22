#!/bin/bash

# recursive script runner


script_dir=$(dirname "$(realpath "$0")")
echo $script_dir
for dir in data_*/; do
    if [ -d "$dir" ]; then
        date_part=$(echo "$dir" | sed 's/data_\([0-9]*_[0-9]*_[0-9]*-[0-9]*\)\//\1/')
        arg="amatter_${date_part}.xyzv"
        cwd=$(pwd)
        file_path="$cwd""/""$dir""$arg"
        echo $file_path
        if ! [ -f $file_path ]; then
            echo "File does not exist."
        else
            (cd "$dir" && OVITO_GUI_MODE=1 python3 "${script_dir}/ovito_render.py" "--verbose" "$arg")
        fi
    fi
done