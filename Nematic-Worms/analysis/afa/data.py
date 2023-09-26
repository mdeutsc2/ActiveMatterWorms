import numpy as np
from scipy.spatial import KDTree
from tqdm import tqdm
import subprocess
from typing import Optional

from .io import load_next_frame, EmptyLine
from .utils import compute_tangents

from multiprocessing.managers import BaseManager


class Frame:
    def __init__(self,
        positions: np.ndarray,
        labels: np.ndarray,
        filament_length: int,
        displacements: Optional[np.ndarray] = None,
        **kwargs
    ):
        if len(positions) != len(labels):
            raise ValueError("Positions and Labels need to be same length.")
        self.r = positions
        self.t: np.ndarray = None
        self.dr = displacements
        self.labels = labels
        self.kd_tree: KDTree = None
        self.np: int = filament_length
        self.metadata = kwargs

    def __len__(self):
        return len(self.r)

    def get(self, indices):
        return self.r[indices]

    def get_t(self, indices):
        return self.t[indices]

    def build_neighbors_list(self):
        self.kd_tree = KDTree(self.r)

    def compute_filament_tangents(self):
        self.t = compute_tangents(self.r, self.np)

    def sample_points(self, n: int):
        sample_indices = np.random.randint(0, len(self), size=n)
        return self.r[sample_indices]

    def particles_near(self, point: np.ndarray, radius: float):
        return np.array(self.kd_tree.query_ball_point(point, radius))


class FilamentFrameLoader:

    def __init__(self,
        filename: str,
        filament_length: int,
        max_frames: int = 1000,
        skip: int = 0,
        compute_displacements=False,
        exclude_n_end_of_frame: int = 0,
    ):
        self.frames = []
        self.max_frames = max_frames
        self.skip = skip
        self.ppf = filament_length
        self.compute_displacements = compute_displacements
        self.exclude_last = exclude_n_end_of_frame
        self.load(filename)


    def __iter__(self) -> "FilamentFrameLoader":
        return iter(self.frames)

    def __len__(self):
        return len(self.frames)

    def __next__(self) -> Frame:
        return next(self.frames)

    def iter_batches(self, batch_size: int):
        for i in range(0, len(self), batch_size):
            yield self.frames[i:i+batch_size]

    def load(self, filename: str):
        with open(filename, "r") as f:
            r_old = None
            for frame_count in tqdm(range(self.skip + self.max_frames), desc=f" Reading {filename}"):
                try:
                    r, t = load_next_frame(f)
                    dr = None

                    if frame_count < self.skip:
                        continue

                    if self.exclude_last > 0:
                        r = r[:-self.exclude_last]
                        t = t[:-self.exclude_last]

                    if self.compute_displacements:
                        if r_old is None:
                            r_old = r
                        dr = r - r_old

                    self.frames.append(Frame(
                        positions=r,
                        displacements=dr,
                        labels=t,
                        filament_length=self.ppf,
                        frame_number=frame_count))
                except EmptyLine:
                    pass



class DefectCollection:

    def __init__(self):
        self.r = []
        self.v = []

    def try_append(self, point, value: float):
        if value > 0:
            self.r.append(point)
            self.v.append(value)

    def __len__(self):
        return len(self.r)

    def to_xyz(self, filebuffer, zshift: float = 0.0):
        filebuffer.write(f"{len(self.r)}\n")
        filebuffer.write(f"DefectCollection\n")
        for _r, _v in zip(self.r, self.v):
            _r[2] += zshift
            filebuffer.write("%s %f %f %f\n" % (_v, *_r))

    @staticmethod
    def merge(collections: list):
        merged = DefectCollection()
        for c in collections:
            merged.r.extend(c.r)
            merged.v.extend(c.v)
        return merged


# BaseManager.register('Frame', Frame)
