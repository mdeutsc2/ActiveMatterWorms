import numpy as np


def winding_number(loop_values: np.ndarray, increment: float):
    values = np.array(loop_values) / increment
    return round(np.sum(np.diff(values)))


BIN_EDGES = {
    n_bins: [(i - n_bins // 2) * (2 * np.pi / n_bins) for i in range(n_bins)]
    for n_bins in range(4, 16)
}

def bins_by_displacement_angle(points: np.ndarray, center: np.ndarray, n_bins: int = 8):
    d = points - center
    angle = np.arctan2(d[:,1], d[:,0])
    return np.digitize(angle, BIN_EDGES[n_bins]) - 1


def strided(a, L):
    # Store the shape and strides info
    shp = a.shape
    s  = a.strides

    # Compute length of output array along the first axis
    nd0 = shp[0]-L+1

    # Setup shape and strides for use with np.lib.stride_tricks.as_strided
    # and get (n+1) dim output array
    shp_in = (nd0,L)+shp[1:]
    strd_in = (s[0],) + s
    return np.lib.stride_tricks.as_strided(a, shape=shp_in, strides=strd_in)


def compute_tangents(particle_positions: np.ndarray, filament_length: int):
    tangents = np.zeros_like(particle_positions)
    for (im1, i, ip1) in strided(np.arange(len(particle_positions)), 3):
        index_in_filament = i  % filament_length

        if index_in_filament == 0: # is head
            im1 = i

        if index_in_filament == (filament_length - 1): # is tail
            ip1 = i

        dr = particle_positions[ip1] - particle_positions[im1]
        tangents[i] = dr / np.linalg.norm(dr)

    return tangents


def compute_director_angle(particle_tangents: np.ndarray):
    r_avg = np.mean(particle_tangents, axis=0)
    return np.arctan2(r_avg[1], r_avg[0])
