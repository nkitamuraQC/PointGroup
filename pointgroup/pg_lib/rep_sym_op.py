
from scipy.spatial.transform import Rotation as R
import numpy as np

import numpy as np
from .transform import (
    tr2o3,
    from_o3,
)
from scipy.spatial.transform import Rotation as R

## cellの有無で変換をスイッチする

def rotation_matrix_axis(axis_vec, angle_deg, cell=None):
    """
    任意の軸(axis_vec: 3次元ベクトル)まわりの回転行列を返す (Rodriguesの公式, scipy使用)
    angle_deg: 回転角（度）
    cell: 格子情報（必要に応じて渡す）
    return: 3x3回転行列
    """
    axis_vec = np.asarray(axis_vec)
    axis_vec = axis_vec / np.linalg.norm(axis_vec)
    rot = R.from_rotvec(np.deg2rad(angle_deg) * axis_vec).as_matrix()
    if cell is not None:
        rot = from_o3(cell, rot)
    return rot

def reflection_matrix_axis(axis_vec, cell):
    n = np.asarray(axis_vec)
    n = n / np.linalg.norm(n)
    rot = np.eye(3) - 2 * np.outer(n, n)
    if cell is not None:
        rot = from_o3(cell, rot)
    return rot

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
