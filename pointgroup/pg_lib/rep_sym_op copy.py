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

def inversion(cell):
    rot = -np.identity(3)
    rot = from_o3(cell, rot)
    return rot



def rotation_matrix(axis, angle_deg, cell):
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
    rot = from_o3(cell, rot)
    return rot

def improper_rotation_matrix(axis, angle_deg, cell):
    """ S_n = C_n followed by σh (mirror perpendicular to axis) """
    rot = reflection_matrix({'x': 'yz', 'y': 'xz', 'z': 'xy'}[axis]) @ rotation_matrix(axis, angle_deg)
    rot = from_o3(cell, rot)
    return rot


def reflection_matrix(plane, cell):
    if plane == 'xy':
        rot = np.diag([1, 1, -1])
    elif plane == 'yz':
        rot = np.diag([-1, 1, 1])
    elif plane == 'xz':
        rot = np.diag([1, -1, 1])

    rot = from_o3(cell, rot)
    return rot

def gen_sigmav2(cell):
    rot = rotation_matrix('z', 120) @ reflection_matrix('xz') @ rotation_matrix('z', -120)
    rot = from_o3(cell, rot)
    return rot

def gen_sigmav2_2(cell):
    rot = rotation_matrix('z', 60) @ reflection_matrix('xz') @ rotation_matrix('z', -60)
    rot = from_o3(cell, rot)
    return rot

def gen_sigmav3(cell):
    rot = rotation_matrix('z', 240) @ reflection_matrix('xz') @ rotation_matrix('z', -240)
    rot = from_o3(cell, rot)
    return rot

def gen_sigmav3_2(cell):
    rot = rotation_matrix('z', 120) @ reflection_matrix('xz') @ rotation_matrix('z', -120)
    rot = from_o3(cell, rot)
    return rot

def gen_c2_d1(cell):
    rot = rotation_matrix('z', 45) @ rotation_matrix('x', 180) @ rotation_matrix('z', -45)
    rot = from_o3(cell, rot)
    return rot


def gen_c2_d2(cell):
    rot = rotation_matrix('z', -45) @ rotation_matrix('x', 180) @ rotation_matrix('z', 45)
    rot = from_o3(cell, rot)
    return rot

def gen_s6(cell, axis, angle_rad):
    axis = np.array([0, 0, 1])
    rotvec = axis * angle_rad  # 回転ベクトル（ロドリゲスベクトル）
    
    # Rotation オブジェクトを生成
    rotation = R.from_rotvec(rotvec)
    
    # 回転行列を取得
    rot = rotation.as_matrix()
    rot = from_o3(cell, rot)
    return rot


def gen_sigmad1(cell):
    rot = rotation_matrix('z', 45) @ reflection_matrix('xz') @ rotation_matrix('z', -45)
    rot = from_o3(cell, rot)
    return rot

def gen_sigmad2(cell):
    rot = rotation_matrix('z', -45) @ reflection_matrix('xz') @ rotation_matrix('z', 45)
    rot = from_o3(cell, rot)
    return rot

def gen_sigmad1_2(cell):
    rot = rotation_matrix('z', 0) @ reflection_matrix('xz'),
    rot = from_o3(cell, rot)
    return rot

def gen_sigmad2_2(cell):
    rot = rotation_matrix('z', 120) @ reflection_matrix('xz') @ rotation_matrix('z', -120)
    rot = from_o3(cell, rot)
    return rot

def gen_sigmad2_3(cell):
    rot = rotation_matrix('z', 60) @ reflection_matrix('xz') @ rotation_matrix('z', -60)
    rot = from_o3(cell, rot)
    return rot

def gen_sigmad3(cell):
    rot = rotation_matrix('z', 240) @ reflection_matrix('xz') @ rotation_matrix('z', -240)
    rot = from_o3(cell, rot)
    return rot


def gen_c2_2(cell):
    rot = rotation_matrix('z', 120) @ rotation_matrix('x', 180) @ rotation_matrix('z', -120)
    rot = from_o3(cell, rot)
    return rot

def gen_c2_3(cell):
    rot = rotation_matrix('z', 240) @ rotation_matrix('x', 180) @ rotation_matrix('z', -240)
    rot = from_o3(cell, rot)
    return rot
