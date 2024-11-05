import numpy as np


def load_next_frame(fp):
    pcount, _ = _read_xyz_header(fp)
    return _read_frame_body(fp, pcount)


def _read_xyz_header(fp):
    line = fp.readline()
    if line == '':
        raise EmptyLine

    particle_count = int(line)
    comment = fp.readline()
    return particle_count, comment


def _read_frame_body(fp, pcount):
    positions = np.empty((pcount, 3))
    labels = np.empty((pcount,), dtype=str)
    for i in range(pcount):
        labels[i], *data = fp.readline().split()
        positions[i] = data[:3]
    return positions, labels


class EmptyLine(Exception):
    pass