import argparse
from collections import defaultdict
from multiprocessing import Pool, cpu_count
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm
import itertools
import numpy as np

from afa.data import FilamentFrameLoader, Frame


def compute_relative_orientations(frame: Frame, n_samples: int, box: list):
    frame.compute_filament_tangents()

    i, j = np.random.randint(0, len(frame), size=[2, n_samples])
    t_i = frame.t[i]
    t_j = frame.t[j]
    t_dot_t = np.sum(t_i * t_j, axis=1)
    n_dot_n = np.abs(t_dot_t)

    Q = 1.5 * np.einsum('...j,...k->...jk', t_i, t_i)
    Q[:, 0, 0] -= 0.5
    Q[:, 1, 1] -= 0.5
    Q[:, 2, 2] -= 0.5

    eigenvalues, eigenvectors = np.linalg.eig(Q.mean(axis=0))
    S = eigenvalues.max()
    P = np.linalg.norm(frame.t.mean(axis=0))

    # Periodic box displacements
    dr = frame.r[i] - frame.r[j]
    dr[:, 0][dr[:, 0] >  box[0]/2] -= box[0]
    dr[:, 0][dr[:, 0] < -box[0]/2] += box[0]
    dr[:, 1][dr[:, 1] >  box[1]/2] -= box[1]
    dr[:, 1][dr[:, 1] < -box[1]/2] += box[1]

    return np.linalg.norm(dr, axis=1), n_dot_n, t_dot_t, S, P



if __name__ == "__main__":

    # Plot order parameters
    # alpha = np.array([0, 0.15, 0.3, 0.45, 0.6])
    # S = [0.384, 0.388, 0.447, 0.462, 0.473]
    # dS = [0.045, 0.048, 0.076, 0.063, 0.103]
    # P = [0.019, 0.347, 0.513, 0.585, 0.561]
    # dP = [0.007, 0.038, 0.106, 0.054, 0.148]

    # plt.errorbar(alpha - 0.0025, S, yerr=dS, label="Nematic", marker="o", linewidth=2)
    # plt.errorbar(alpha + 0.0025, P, yerr=dP, label="Polar", marker="s", linewidth=2)
    # plt.xlabel("Activity")
    # plt.ylabel("Order Parameter")
    # plt.legend(loc="lower right")
    # plt.show()
    # exit()

    # Plot persistence lengths
    # alpha = np.array([0, 0.15, 0.3, 0.45, 0.6])
    # l_s = [10, 13, 22, 24, 29]
    # l_p = [0.5, 6.2, 16.75, 27.5, 31.78]

    # plt.plot(alpha, l_s, label="Nematic", marker="o", linewidth=2)
    # plt.plot(alpha, l_p, label="Polar", marker="s", linewidth=2)
    # plt.xlabel("Activity")
    # plt.ylabel("Correlation Length")
    # plt.legend(loc="lower right")
    # plt.show()
    # exit()


    parser = argparse.ArgumentParser(description="Compute orientational correlation functions from XYZ files of active filament simulation.")
    parser.add_argument("-i", "--inputs", nargs='+', type=str, required=True)
    parser.add_argument("-n", "--filament_length", type=int, required=True)
    parser.add_argument("--bin_size", type=float, default=1)
    parser.add_argument("-b", "--box", nargs='+', type=float, required=True)
    parser.add_argument("-s", "--skip", type=int, default=0)
    parser.add_argument("-m", "--samples", type=float, default=10000, help="For each frame, sample m pairs of orientations to update running function estimate")
    parser.add_argument("-l", "--limit", type=int, default=1000000)
    parser.add_argument("-np", "--num_processes", type=int, default=cpu_count() + 1)
    args = parser.parse_args()

    for input_file in args.inputs:

        loader = FilamentFrameLoader(input_file, args.filament_length, args.limit, args.skip, exclude_n_end_of_frame=4)

        X_nematic = defaultdict(list)
        X_polar = defaultdict(list)
        S = []
        P = []

        # https://pubs-rsc-org.proxy.library.kent.edu/en/content/articlelanding/2016/sm/c6sm00812g
        def nematic_distribution(X):
            data = np.array([(i*args.bin_size, 2*np.mean(values) - 1) for i,values in X.items()])
            order = np.argsort(data[:,0])
            return data[order,0], data[order,1]

        def polar_distribution(X):
            data = np.array([(i*args.bin_size, np.mean(values)) for i,values in X.items()])
            order = np.argsort(data[:,0])
            return data[order,0], data[order,1]

        # with open(args.output, "w+") as out:
        with Pool(processes=args.num_processes) as pool:
            for frame_batch in tqdm(loader.iter_batches(args.num_processes), desc=f" Progress"):
                process_args = [(frame, args.samples, args.box) for frame in frame_batch]

                for R, nn_values, tt_values, s, p in pool.starmap(compute_relative_orientations, process_args):

                    S.append(s)
                    P.append(p)

                    # convert to bins
                    for i in range(len(R)):
                        if R[i] < args.box[0] / 2:
                            X_nematic[R[i] // args.bin_size].append(nn_values[i])
                            X_polar[R[i] // args.bin_size].append(tt_values[i])

            print("Nematic Order (S) =", round(np.mean(S), 3), "+\-", round(np.std(S), 3))
            print("Polar Order   (P) =", round(np.mean(P), 3), "+\-", round(np.std(P), 3))

            # x, y = nematic_distribution(X_nematic)
            # plt.plot(x, y, label=input_file + " Nematic", linewidth=2)
            x, y = polar_distribution(X_polar)
            plt.plot(x, y, label=input_file + " Polar", linewidth=2)

    plt.legend()
    plt.show()


