#!/bin/bash

set -e
cd ./data_osc/
cd ./data_01_01_2024/
python3 ../../scripts/merge_xyz.py amatter_osc_01_01_2024.xyzv
OVITO_GUI_MODE=1 python3 ../../scripts/ovito_render.py --reencode amatter_osc_01_01_2024.xyzv
cd ../data_01_01_2024-2/
python3 ../../scripts/merge_xyz.py amatter_osc_01_01_2024-2.xyzv
OVITO_GUI_MODE=1 python3 ../../scripts/ovito_render.py --reencode amatter_osc_01_01_2024-2.xyzv
cd ../data_01_01_2024-3/
python3 ../../scripts/merge_xyz.py amatter_osc_01_01_2024-3.xyzv
OVITO_GUI_MODE=1 python3 ../../scripts/ovito_render.py --reencode amatter_osc_01_01_2024-3.xyzv
cd ../data_01_01_2024-4/
python3 ../../scripts/merge_xyz.py amatter_osc_01_01_2024-4.xyzv
OVITO_GUI_MODE=1 python3 ../../scripts/ovito_render.py --reencode amatter_osc_01_01_2024-4.xyzv
cd ../data_01_01_2024-5/
python3 ../../scripts/merge_xyz.py amatter_osc_01_01_2024-5.xyzv
OVITO_GUI_MODE=1 python3 ../../scripts/ovito_render.py --reencode amatter_osc_01_01_2024-5.xyzv
cd ../data_01_01_2024-6/
python3 ../../scripts/merge_xyz.py amatter_osc_01_01_2024-6.xyzv
OVITO_GUI_MODE=1 python3 ../../scripts/ovito_render.py --reencode amatter_osc_01_01_2024-6.xyzv