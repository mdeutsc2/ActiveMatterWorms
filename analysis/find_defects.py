import argparse
from multiprocessing import cpu_count
import math

from tqdm import tqdm
import numpy as np
from sklearn.cluster import AffinityPropagation, MeanShift, OPTICS, DBSCAN
from sklearn.exceptions import ConvergenceWarning


def read_vector_field(filename):
    frames = []
    with open(filename, "r") as f:
        count = 0
        while True:
            line = f.readline()
            if line in ["", None]:
                break

            # TODO infer format of line based on number of elements
            count += 1
            num_sites = int(line)
            comment = f.readline()
            R = np.zeros((num_sites, 3))
            V = np.zeros((num_sites, 3))
            for i in range(num_sites):
                data = [float(x) for x in f.readline().split()[1:]]
                R[i] = data[:3]
                V[i] = data[3:6]

            origin = R.min(axis=0)
            frames.append({
                "index": count,
                "num_sites": num_sites,
                "positions": R,
                "director": V,
                "origin": origin,
                "box_size": R.max(axis=0) - origin,
            })
    return frames


def compute_theta(frame, cellsize):
    num_sites = frame["num_sites"]
    R = frame["positions"]
    V = frame["director"]
    origin = frame["origin"]
    box_size = frame["box_size"]
    frame["cellsize"] = cellsize
    box_shape = tuple(np.ceil((box_size + 0.0001) / cellsize).astype(np.int64))
    theta = np.zeros(box_shape[:2])
    for i in range(num_sites):
        r, v = R[i], V[i]
        cell = np.floor((r - origin) / cellsize).astype(np.int64)
        theta[cell[0], cell[1]] = math.atan2(v[1], v[0])
    frame["theta"] = theta
    return frame


def theta_diff(t1, t2):
    dt = t2 - t1
    if dt >= math.pi / 2.:
        dt -= math.pi
    elif dt <= -math.pi / 2.:
        dt += math.pi
    return dt

def burgers_circuit(i, j, theta, loop):
    circuit_angle = 0
    di, dj = loop[0]
    t1 = theta[i + di, j + dj]
    for di,dj in loop[1:]:
        t2 = theta[i + di, j + dj]
        circuit_angle += theta_diff(t1, t2)
        t1 = t2
    return circuit_angle

DELTA = math.pi / 16


def contains_defect(angle):
    return angle > math.pi - DELTA or angle < -math.pi + DELTA


LOOP_A = [
    [0, 0],
    [1, 0],
    [1, 1],
    [0, 1],
    [0, 0],
]

LOOP_B = [
    [1, 0],
    [1, 1],
    [0, 1],
    [-1, 1],
    [-1, 0],
    [-1, -1],
    [0, -1],
    [1, -1],
    [1, 0],
]

LOOP_C = [
    [2, 0],
    [1, 1],
    [0, 2],
    [-1, 1],
    [-2, 0],
    [-1, -1],
    [0, -2],
    [1, -1],
    [2, 0],
]

INTERIOR_C = [
    [0, 0],
    [1, 0],
    [0, 1],
    [-1, 0],
    [0, -1],
]


def _cluster(counts):
    labels = np.zeros_like(counts)
    N,M = counts.shape

    # i,j grid positions of where defects are.
    I,J = np.where(counts > 0)

    # Initialize each defect w/ unique index.
    id = 1
    for i,j in zip(I,J):
        labels[i, j] = id
        id += 1

    # Iterate over each defect and merge w/ neighbor if it has lower index.
    def _replace_with_lowest_neighbor_id():
        swaps = 0
        for i,j in zip(I,J):
            for di,dj in LOOP_B[1:]:
                ni = i + di
                nj = j + dj
                # No periodic boundarys are implied here.
                # TODO: Implement clustering over PBC
                if ni >= 0 and ni < N and nj >= 0 and nj < M and labels[ni,nj] > 0 and labels[ni,nj] < labels[i,j]:
                    labels[i,j] = labels[ni,nj]
                    swaps += 1
        return swaps > 0

    # Do this until no more defects are renumbered.
    while _replace_with_lowest_neighbor_id():
        pass

    # Provide the x,y position of each cluster.
    clusters = []
    for cluster_id in np.unique(labels[labels > 0]):
        x,y = np.where(labels == cluster_id)
        clusters.append([np.mean(x), np.mean(y)])
    return clusters


def find_defects(frame):
    theta = frame["theta"]
    cellsize = frame["cellsize"]
    origin = frame["origin"]
    counts_pos = np.zeros_like(theta)
    counts_neg = np.zeros_like(theta)
    charge = {}


    for i in range(1, theta.shape[0] - 1):
        for j in range(1, theta.shape[1] - 1):
            # angle_circuit_a = burgers_circuit(i, j, theta, LOOP_A)
            # if contains_defect(angle_circuit_a):
            #     for di, dj in LOOP_A:
            #         counts[i + di, j + dj] += 1
            #         charge.setdefault((i + di, j + dj), []).append(angle_circuit_a / (2*math.pi))

            angle_circuit_b = burgers_circuit(i, j, theta, LOOP_B)
            if contains_defect(angle_circuit_b):
                q = angle_circuit_b / (2*math.pi)
                charge.setdefault((i, j), []).append(q)
                if q > 0:
                    counts_pos[i, j] += 1
                elif q < 0:
                    counts_neg[i, j] += 1

            # angle_circuit_c = burgers_circuit(i, j, theta, LOOP_C)
            # if contains_defect(angle_circuit_c):
            #     for di, dj in INTERIOR_C:
            #         counts[i + di, j + dj] += 1
            #         charge.setdefault((i + di, j + dj), []).append(angle_circuit_c / (2*math.pi))

    charge = {k: np.mean(v) for k,v in charge.items()}

    defects = []
    if counts_pos.sum() > 0:
        for x,y in _cluster(counts_pos):
            defects.append(("A", x*cellsize + origin[0], y*cellsize + origin[1], 0, 0.5))

    if counts_neg.sum() > 0:
        for x,y in _cluster(counts_neg):
            defects.append(("B", x*cellsize + origin[0], y*cellsize + origin[1], 0, -0.5))
    return defects


def write_defects(defects, file_handle):
    if len(defects) > 0:
        file_handle.write(f'{len(defects)}\n')
        file_handle.write(f'Defects detected\n')
        for c,x,y,z,q in defects:
            file_handle.write("%c %f %f %f %f\n" % (c, x, y, z, q))
    else:
        origin = frame["origin"]
        file_handle.write(f'2\n')
        file_handle.write(f'No defects detected\n')
        file_handle.write(f'A {origin[0]} {origin[1]} 0 0.5\n')
        file_handle.write(f'B {origin[0]} {origin[1]} 0 -0.5\n')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analyze XYZ files of active filament simulation.")
    parser.add_argument("-i", "--input", type=str, required=True, help="A director field sequence XYZV file")
    parser.add_argument("-o", "--output", type=str, default="defects.xyz")
    parser.add_argument("-np", "--num_processes", type=int, default=cpu_count() + 1)
    parser.add_argument("-c", "--cellsize", type=float, required=True) # 2 in make_directory_field.py
    args = parser.parse_args()

    print(f"Reading {args.input} ...")
    with open(args.output, "w+") as out:
        for frame in tqdm(read_vector_field(args.input), desc=" Progress"):
            compute_theta(frame, args.cellsize)
            defects = find_defects(frame)
            write_defects(defects, out)
