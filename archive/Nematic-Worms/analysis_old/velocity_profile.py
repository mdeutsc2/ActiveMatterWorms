import argparse
from multiprocessing import Pool, cpu_count
from tqdm import tqdm
import itertools
import numpy as np

from afa.data import FilamentFrameLoader, Frame



class VelocityProfileByAxisAnalysis:

    def __init__(self, axis, filename, bins=10):
        self._axis = axis
        self._bins = bins
        self._filename = filename

    def update(self, frame):




    def analyze(self):
        # print("PositionByAxis Analyzing...")
        # for n, (label, positions) in enumerate(self._positions_by_charge.items()):
        #     print(f"{len(positions)} defects of label {label}")
        #     lab = self._type_labels[n] if self._type_labels is not None else label
        #     plt.hist(positions, bins=self._bins, density=True, label=f"{lab}", orientation="horizontal", alpha=0.6)
        # plt.xlabel("Density")
        # plt.ylabel("Position")
        # plt.title("Defect Positions")
        # plt.legend()
        # plt.savefig(self._filename)



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute and analyze the velocities of filament XYZ files.")
    parser.add_argument("-i", "--input", type=str, required=True)
    parser.add_argument("-o", "--output", type=str, default="velocity_profile.png")
    parser.add_argument('-a', '--profile-axis', type=int, choices=[0, 1, 2], default=0, help="Axis to compute velocity profile over.")
    parser.add_argument('-n', '--nbins', type=int, default=10, help="Number of bins to analyze using.")
    parser.add_argument("-s", "--skip", type=int, default=0)
    parser.add_argument("-np", "--num_processes", type=int, default=cpu_count() + 1)
    args = parser.parse_args()

    calculations = [
        VelocityProfileByAxisAnalysis(args.profile-axis, args.output, args.nbins)
    ]

    loader = FilamentFrameLoader(args.input, args.filament_length, args.limit, args.skip)

