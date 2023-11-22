import argparse
from multiprocessing import Pool, cpu_count
from tqdm import tqdm
import itertools
import numpy as np

from afa.data import FilamentFrameLoader, Frame


def compute_director(particle_frame: Frame, cellsize):
    particle_frame.compute_filament_tangents()

    origin = particle_frame.r.min(axis=0)
    box_max = particle_frame.r.max(axis=0)
    box_size = box_max - origin
    box_shape = tuple(np.ceil((box_size + 0.0001) / cellsize).astype(np.int64))
    n_cells = np.prod(box_shape)

    N = np.zeros(box_shape, dtype=int)
    Q = np.zeros(box_shape + (3, 3))
    S = np.zeros_like(N, dtype=float)
    R = np.zeros(box_shape + (3,))
    director = np.zeros(box_shape + (3,))

    # For each particle, update the Q for the cell it is in.
    for i in range(len(particle_frame)):
        r = particle_frame.get(i)
        t = particle_frame.get_t(i)
        cell = np.floor((r - origin) / cellsize).astype(np.int64)
        N[cell[0], cell[1], cell[2]] += 1

        for a in range(3):
            for b in range(3):
                Q[cell[0], cell[1], cell[2], a, b] += 3*t[a]*t[b] - int(a == b)

    # Compute the position of each cell and compute the director
    # of all cells containing particles
    for cell in itertools.product(*[range(s) for s in box_shape]):
        R[cell] = origin + (np.array(cell) * cellsize) + cellsize / 2

        if N[cell] > 0:
            Q[cell] /= 2 * N[cell] # Complete the Q calculation
            eigenvalues, eigenvectors = np.linalg.eig(Q[cell])
            S[cell] = eigenvalues.max()
            director[cell] = eigenvectors[eigenvalues.argmax()]

    return np.hstack((R.reshape(n_cells, 3), director.reshape(n_cells, 3), S.reshape(n_cells, 1), N.reshape(n_cells, 1)))




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analyze XYZ files of active filament simulation.")
    parser.add_argument("-i", "--input", type=str, required=True)
    parser.add_argument("-n", "--filament_length", type=int, required=True)
    parser.add_argument("-o", "--output", type=str, default="defects.xyz")
    parser.add_argument("-s", "--skip", type=int, default=0)
    parser.add_argument("-c", "--cellsize", type=float, default=2)
    parser.add_argument("-l", "--limit", type=int, default=1000000)
    parser.add_argument("-np", "--num_processes", type=int, default=cpu_count() + 1)
    args = parser.parse_args()

    loader = FilamentFrameLoader(args.input, args.filament_length, args.limit, args.skip)

    with open(args.output, "w+") as out:
        with Pool(processes=args.num_processes) as pool:
            for frame_batch in tqdm(loader.iter_batches(args.num_processes), desc=f" Progress"):
                process_args = [(frame, args.cellsize) for frame in frame_batch]
                results = pool.starmap(compute_director, process_args)
                for director_field in results:
                    out.write(f"# {len(director_field)}\n")
                    out.write("# type,x,y,z,nx,ny,nz,s,N\n")
                    for x in director_field:
                        x,y,z,vx,vy,vz,s,n = tuple(x)
                        out.write("A %f %f %f %f %f %f %f %f\n" % (x,y,z,vx,-vy,vz,s,n))
