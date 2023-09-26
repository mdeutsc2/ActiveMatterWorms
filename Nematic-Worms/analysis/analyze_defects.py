import argparse
from multiprocessing import cpu_count

import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm


def read_input(filename, skip_count):
    frames = []
    with open(filename, "r") as f:
        count = 0

        while True:
            line = f.readline()
            if line in ["", None]:
                break

            count += 1
            num_sites = int(line)
            comment = f.readline()

            R = np.zeros((num_sites, 3))
            L = np.empty(num_sites, dtype=str)
            Q = np.zeros(num_sites)
            for i in range(num_sites):
                data = f.readline().split()
                L[i] = data[0]
                R[i] = [float(x) for x in data[1:4]]
                Q[i] = float(data[-1])

            if count <= skip_count:
                continue

            frames.append({
                "index": count,
                "num_sites": num_sites,
                "labels": L,
                "positions": R,
                "charge": Q,
                "contains_defects": comment != "No defects detected",
            })
    return frames


class PositionByAxisAnalysis:

    def __init__(self, axis, filename, bins=10, type_labels=None):
        self._axis = axis
        self._bins = bins
        self._filename = filename
        self._positions_by_charge = {}
        self._type_labels = type_labels

    def update(self, frame):
        R = frame["positions"]
        Q = frame["charge"]
        L = frame["labels"]

        if len(Q) != len(R) != len(L):
            raise ValueError("R, Q, L are required to be the same length")

        if R.shape[1] <= self._axis:
            raise ValueError(f"Axis {self._axis} cannot be used on positions of shape {R.shape}")

        for i in range(len(R)):
            self._positions_by_charge.setdefault(L[i], []).append(R[i, self._axis])

    def analyze(self):
        print("PositionByAxis Analyzing...")
        for n, (label, positions) in enumerate(self._positions_by_charge.items()):
            print(f"{len(positions)} defects of label {label}")
            lab = self._type_labels[n] if self._type_labels is not None else label
            plt.hist(positions, bins=self._bins, density=True, label=f"{lab}", orientation="horizontal", alpha=0.6)
        plt.xlabel("Density")
        plt.ylabel("Position")
        plt.title("Defect Positions")
        plt.legend()
        plt.savefig(self._filename)


class NumberOfDefectsAnalysis:
    def __init__(self, filename, type_labels=None):
        self._filename = filename
        self._type_labels = type_labels
        self._counts = {}

    def update(self, frame):
        R = frame["positions"]
        Q = frame["charge"]
        L = frame["labels"]

        if len(Q) != len(R) != len(L):
            raise ValueError("R, Q, L are required to be the same length")

        _counts = {}
        for i in range(len(R)):
            _counts.setdefault(L[i], 0)
            _counts[L[i]] += 1

        for label,counts in _counts.items():
            self._counts.setdefault(label, []).append(counts)

    def analyze(self):
        print("NumberOfDefects Analyzing...")
        for n, (label, counts) in enumerate(self._counts.items()):
            plt.plot(counts, label=f"{self._type_labels[n]}")
        plt.xlabel("Time Steps")
        plt.ylabel("Defect Counts")
        plt.legend()
        plt.savefig(self._filename)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analyze XYZ sequence of defect positions.")
    parser.add_argument("-i", "--input", type=str, required=True, help="A director field sequence XYZV file")
    parser.add_argument("-o", "--output", type=str, default="defects.png")
    parser.add_argument('--position-by-axis', type=int, choices=[0, 1, 2], default=0, help="Axis to analyze positions by.")
    parser.add_argument('--nbins', type=int, default=10, help="Number of bins to analyze using.")
    parser.add_argument('--labels', type=str, default=["+1/2", "-1/2"], help="Labels for defects of type A, B, etc")
    parser.add_argument("-s", "--skip", type=int, default=0)
    parser.add_argument("-np", "--num_processes", type=int, default=cpu_count() + 1)
    args = parser.parse_args()

    calculations = [
        # NumberOfDefectsAnalysis(args.output, args.labels),
        PositionByAxisAnalysis(args.position_by_axis, args.output, args.nbins, args.labels)
    ]

    print(f"Reading {args.input} ...")
    for frame in tqdm(read_input(args.input, args.skip), desc=" Progress"):
        if frame["contains_defects"]:
            for calculation in calculations:
                calculation.update(frame)

    for calculation in calculations:
        calculation.analyze()