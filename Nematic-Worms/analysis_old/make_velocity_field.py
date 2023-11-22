import argparse
from multiprocessing import Pool, cpu_count
from tqdm import tqdm
import itertools
import numpy as np

from afa.data import FilamentFrameLoader, Frame


def compute_velocity_field(particle_frame: Frame, cellsize):
    origin = particle_frame.r.min(axis=0)
    box_max = particle_frame.r.max(axis=0)
    box_size = box_max - origin
    box_shape = tuple(np.ceil((box_size + 0.0001) / cellsize).astype(np.int))
    n_cells = np.prod(box_shape)

    N = np.zeros(box_shape, dtype=int)
    R = np.zeros(box_shape + (3,))
    V = np.zeros(box_shape + (3,))

    # For each particle, sum velocity of particles in for the cell
    for i in range(len(particle_frame)):
        r = particle_frame.get(i)
        v = particle_frame.dr[i]

        cell = np.floor((r - origin) / cellsize).astype(np.int)
        N[cell[0], cell[1], cell[2]] += 1
        V[cell[0], cell[1], cell[2]] += v

    # Give a position to each cell and find maximum velocity in system to normalize
    max_magnitude = 0
    for cell in itertools.product(*[range(s) for s in box_shape]):
        R[cell] = origin + (np.array(cell) * cellsize) + cellsize / 2
        v_mag = np.linalg.norm(V[cell])
        if v_mag > max_magnitude:
            max_magnitude = v_mag
    # V /= max_magnitude

    return np.hstack((R.reshape(n_cells, 3), V.reshape(n_cells, 3)))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create a velocity profile of nematic filament XYZ files.")
    parser.add_argument("-i", "--input", type=str, required=True)
    parser.add_argument("-n", "--filament_length", type=int, required=True)
    parser.add_argument("-o", "--output", type=str, default="velocity_field.xyzv")
    parser.add_argument("-c", "--cellsize", type=float, default=2)
    parser.add_argument("-l", "--limit", type=int, default=1000000)
    parser.add_argument("-s", "--skip", type=int, default=0)
    parser.add_argument("-np", "--num_processes", type=int, default=cpu_count() + 1)
    args = parser.parse_args()

    loader = FilamentFrameLoader(args.input, args.filament_length, args.limit, args.skip, compute_displacements=True)

    with open(args.output, "w+") as out:
        with Pool(processes=args.num_processes) as pool:
            for frame_batch in tqdm(loader.iter_batches(args.num_processes), desc=f" Progress"):
                process_args = [(frame, args.cellsize) for frame in frame_batch]
                results = pool.starmap(compute_velocity_field, process_args)
                for velocity_field in results:
                    out.write(f"{len(velocity_field)}\n")
                    out.write("type,x,y,z,vx,vy,vz\n")
                    for x in velocity_field:
                        x,y,z,vx,vy,vz = tuple(x)
                        out.write("A %f %f %f %f %f %f\n" % (x,y,z,vx,-vy,vz))