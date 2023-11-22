import os
import argparse
from multiprocessing import Pool, Manager

import numpy as np
from tqdm import tqdm

from afa.data import FilamentFrameLoader, DefectCollection, Frame
from afa.utils import *


def search_for_defects(particle_frame: Frame, params):
    pid = os.getpid()
    radius = params["radius"]
    n_points = params["n_points"]
    defect_collection = DefectCollection()

    for test_point in particle_frame.sample_points(n_points):
        particles_indices_near_point = particle_frame.particles_near(test_point, radius)
        particles_near_point = particle_frame.get(particles_indices_near_point)
        particle_bin_assignments = bins_by_displacement_angle(particles_near_point, test_point)
        bin_ids = np.unique(particle_bin_assignments)
        director_angles = np.zeros((max(bin_ids+1),))

        for bin_id in bin_ids:
            particles_accessor = particles_indices_near_point[particle_bin_assignments == bin_id]
            directory_angle = compute_director_angle(particle_frame.get_t(particles_accessor))
            director_angles[bin_id] = directory_angle

        attempt_value = winding_number(director_angles, np.pi)
        defect_collection.try_append(test_point, attempt_value)
    return defect_collection


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analyze XYZ files of active filament simulation.")
    parser.add_argument("-i", "--input", type=str, required=True)
    parser.add_argument("-o", "--output", type=str, default="defects.xyz")
    parser.add_argument("-s", "--skip", type=int, default=0)
    parser.add_argument("-l", "--limit", type=int, default=1000)
    parser.add_argument("-n", "--filament_length", type=int, required=True)
    parser.add_argument("-a", "--attempts", type=int, default=100)
    parser.add_argument("-r", "--compute_radius", type=float, default=10)
    parser.add_argument("-np", "--num_processes", type=int, default=4)
    args = parser.parse_args()

    defect_collections = []

    with Manager() as manager:
        params = manager.dict({ "n_points": args.attempts, "radius": args.compute_radius })
        with Pool(processes=args.num_processes) as pool:
            for particle_frame in tqdm(FilamentFrameLoader(args.input, args.filament_length, args.limit, args.skip), desc=f" Progress"):
                particle_frame.build_neighbors_list()
                particle_frame.compute_filament_tangents()
                process_args = [(particle_frame, params)] * args.num_processes
                defect_collection = DefectCollection.merge(pool.starmap(search_for_defects, process_args))
                defect_collections.append(defect_collection)

    with open(args.output, "w") as f:
        [dc.to_xyz(f, zshift=1.0) for dc in defect_collections]