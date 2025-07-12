from scipy.spatial.transform import Rotation as R
import numpy as np

import numpy as np
from .transform import (
    tr2o3,
    from_o3,
)
from scipy.spatial.transform import Rotation as R

## cellの有無で変換をスイッチする

def identity():
    return np.identity(3)


def inversion(cell=None):
    rot = -np.identity(3)
    if cell is not None:
        rot = from_o3(cell, rot)
    return rot



def rotation_matrix(axis, angle_deg, cell=None):
    angle_rad = np.radians(angle_deg)
    c, s = np.cos(angle_rad), np.sin(angle_rad)
    if axis == 'x':
        rot = np.array([[1, 0, 0],
                        [0, c, -s],
                        [0, s, c]])
    elif axis == 'y':
        rot = np.array([[c, 0, s],
                         [0, 1, 0],
                         [-s, 0, c]])
    elif axis == 'z':
        rot = np.array([[c, -s, 0],
                         [s, c, 0],
                         [0, 0, 1]])
    if cell is not None:
        rot = from_o3(cell, rot)
    return rot

def improper_rotation_matrix(axis, angle_deg, cell=None):
    """ S_n = C_n followed by σh (mirror perpendicular to axis) """
    rot = reflection_matrix({'x': 'yz', 'y': 'xz', 'z': 'xy'}[axis]) @ rotation_matrix(axis, angle_deg)
    if cell is not None:
        rot = from_o3(cell, rot)
    return rot


def reflection_matrix(plane, cell=None):
    if plane == 'xy':
        rot = np.diag([1, 1, -1])
    elif plane == 'yz':
        rot = np.diag([-1, 1, 1])
    elif plane == 'xz':
        rot = np.diag([1, -1, 1])
    
    if cell is not None:
        rot = from_o3(cell, rot)
    return rot
