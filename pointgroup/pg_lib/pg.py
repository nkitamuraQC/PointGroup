### 主要な32点群のうち、代表的な点群ごとに対称操作の行列表現を辞書で整理します。
# 各点群に対応する操作は簡潔な3x3行列（直交変換）として与えます。

import numpy as np

# 基本操作の行列表現（軸: x=0, y=1, z=2）
E = np.identity(3)

def inversion():
    return -np.identity(3)

def improper_rotation_matrix(axis, angle_deg):
    """ S_n = C_n followed by σh (mirror perpendicular to axis) """
    return reflection_matrix({'x': 'yz', 'y': 'xz', 'z': 'xy'}[axis]) @ rotation_matrix(axis, angle_deg)


def rotation_matrix(axis, angle_deg):
    angle_rad = np.radians(angle_deg)
    c, s = np.cos(angle_rad), np.sin(angle_rad)
    if axis == 'x':
        return np.array([[1, 0, 0],
                         [0, c, -s],
                         [0, s, c]])
    elif axis == 'y':
        return np.array([[c, 0, s],
                         [0, 1, 0],
                         [-s, 0, c]])
    elif axis == 'z':
        return np.array([[c, -s, 0],
                         [s, c, 0],
                         [0, 0, 1]])

def reflection_matrix(plane):
    if plane == 'xy':
        return np.diag([1, 1, -1])
    elif plane == 'yz':
        return np.diag([-1, 1, 1])
    elif plane == 'xz':
        return np.diag([1, -1, 1])

# 代表的な点群の操作（例: C1, Ci, Cs, C2, C2v, D2, D2h, C3, C3v, D3, D3h, Td, Oh, Ih）
# 一部を例示し、残りは後で補完可能です

point_groups = {}

# C1 (identity only)
point_groups['C1'] = {'E': E}

# Ci (identity, inversion)
point_groups['Ci'] = {
    'E': E,
    'i': -np.identity(3)
}

# Cs (mirror plane, e.g. xz)
point_groups['Cs'] = {
    'E': E,
    'σ(xz)': reflection_matrix('xz')
}

# C2 (rotation 180° around z)
point_groups['C2'] = {
    'E': E,
    'C2(z)': rotation_matrix('z', 180)
}

# C2v (E, C2, σv(xz), σv'(yz))
point_groups['C2v'] = {
    'E': E,
    'C2(z)': rotation_matrix('z', 180),
    'σ(xz)': reflection_matrix('xz'),
    'σ(yz)': reflection_matrix('yz')
}

# D2 (E, C2(x), C2(y), C2(z))
point_groups['D2'] = {
    'E': E,
    'C2(x)': rotation_matrix('x', 180),
    'C2(y)': rotation_matrix('y', 180),
    'C2(z)': rotation_matrix('z', 180)
}

# C3 (E, C3, C3^2)
C3 = rotation_matrix('z', 120)
C3_2 = rotation_matrix('z', 240)
point_groups['C3'] = {
    'E': E,
    'C3(z)': C3,
    'C3^2(z)': C3_2
}

# C3v (E, 2C3, 3σv)
point_groups['C3v'] = {
    'E': E,
    'C3(z)': C3,
    'C3^2(z)': C3_2,
    'σv1': reflection_matrix('xz'),
    'σv2': rotation_matrix('z', 120) @ reflection_matrix('xz') @ rotation_matrix('z', -120),
    'σv3': rotation_matrix('z', 240) @ reflection_matrix('xz') @ rotation_matrix('z', -240)
}

# Td (tetrahedral)
# Approximating with identity and representative operations
point_groups['Td'] = {
    'E': E,
    'C3': rotation_matrix('z', 120),
    'C2': rotation_matrix('x', 180),
    'S4': rotation_matrix('z', 90),  # improper
    'σd': reflection_matrix('xz')
}

# Oh (octahedral)
point_groups['Oh'] = {
    'E': E,
    'C4(z)': rotation_matrix('z', 90),
    'C3': rotation_matrix('z', 120),
    'C2(x)': rotation_matrix('x', 180),
    'i': -np.identity(3)
}

# point_groups.keys()

# 残りの32点群の対称操作を追加して完成させる

# 補助関数（再掲）

# Monoclinic point groups
point_groups['C2h'] = {
    'E': E,
    'C2(z)': rotation_matrix('z', 180),
    'i': inversion(),
    'σh(xy)': reflection_matrix('xy')
}

# Orthorhombic
point_groups['D2h'] = {
    'E': E,
    'C2(x)': rotation_matrix('x', 180),
    'C2(y)': rotation_matrix('y', 180),
    'C2(z)': rotation_matrix('z', 180),
    'i': inversion(),
    'σ(xy)': reflection_matrix('xy'),
    'σ(yz)': reflection_matrix('yz'),
    'σ(xz)': reflection_matrix('xz')
}

# Tetragonal point groups
point_groups['C4'] = {
    'E': E,
    'C4(z)': rotation_matrix('z', 90),
    'C4^2(z)': rotation_matrix('z', 180),
    'C4^3(z)': rotation_matrix('z', 270)
}

point_groups['S4'] = {
    'E': E,
    'S4': improper_rotation_matrix('z', 90),
    'S4^2': rotation_matrix('z', 180),
    'S4^3': improper_rotation_matrix('z', 270)
}

point_groups['C4h'] = {
    'E': E,
    'C4(z)': rotation_matrix('z', 90),
    'C2(z)': rotation_matrix('z', 180),
    'C4^3(z)': rotation_matrix('z', 270),
    'i': inversion(),
    'S4': improper_rotation_matrix('z', 90),
    'S4^3': improper_rotation_matrix('z', 270),
    'σh(xy)': reflection_matrix('xy')
}

point_groups['D4'] = {
    'E': E,
    'C4(z)': rotation_matrix('z', 90),
    'C2(z)': rotation_matrix('z', 180),
    'C4^3(z)': rotation_matrix('z', 270),
    'C2(x)': rotation_matrix('x', 180),
    'C2(y)': rotation_matrix('y', 180),
    'C2(d1)': rotation_matrix('z', 45) @ rotation_matrix('x', 180) @ rotation_matrix('z', -45),
    'C2(d2)': rotation_matrix('z', -45) @ rotation_matrix('x', 180) @ rotation_matrix('z', 45)
}

point_groups['C4v'] = {
    'E': E,
    'C4(z)': rotation_matrix('z', 90),
    'C2(z)': rotation_matrix('z', 180),
    'C4^3(z)': rotation_matrix('z', 270),
    'σv(xz)': reflection_matrix('xz'),
    'σv(yz)': reflection_matrix('yz'),
    'σd1': rotation_matrix('z', 45) @ reflection_matrix('xz') @ rotation_matrix('z', -45),
    'σd2': rotation_matrix('z', -45) @ reflection_matrix('xz') @ rotation_matrix('z', 45)
}

point_groups['D2d'] = {
    'E': E,
    'C2(z)': rotation_matrix('z', 180),
    'C2(x)': rotation_matrix('x', 180),
    'C2(y)': rotation_matrix('y', 180),
    'S4': improper_rotation_matrix('z', 90),
    'S4^3': improper_rotation_matrix('z', 270),
    'σd1': rotation_matrix('z', 45) @ reflection_matrix('xz') @ rotation_matrix('z', -45),
    'σd2': rotation_matrix('z', -45) @ reflection_matrix('xz') @ rotation_matrix('z', 45)
}

point_groups['D4h'] = {
    'E': E,
    'C4(z)': rotation_matrix('z', 90),
    'C2(z)': rotation_matrix('z', 180),
    'C4^3(z)': rotation_matrix('z', 270),
    'C2(x)': rotation_matrix('x', 180),
    'C2(y)': rotation_matrix('y', 180),
    'i': inversion(),
    'S4': improper_rotation_matrix('z', 90),
    'S4^3': improper_rotation_matrix('z', 270),
    'σh': reflection_matrix('xy'),
    'σv(xz)': reflection_matrix('xz'),
    'σv(yz)': reflection_matrix('yz')
}

# Trigonal point groups
point_groups['S6'] = {
    'E': E,
    'S6': improper_rotation_matrix('z', 60),
    'S6^5': improper_rotation_matrix('z', 300),
    'C3': rotation_matrix('z', 120),
    'C3^2': rotation_matrix('z', 240),
    'i': inversion()
}

point_groups['D3'] = {
    'E': E,
    'C3': rotation_matrix('z', 120),
    'C3^2': rotation_matrix('z', 240),
    'C2(1)': rotation_matrix('x', 180),
    'C2(2)': rotation_matrix('z', 120) @ rotation_matrix('x', 180) @ rotation_matrix('z', -120),
    'C2(3)': rotation_matrix('z', 240) @ rotation_matrix('x', 180) @ rotation_matrix('z', -240)
}

point_groups['D3d'] = {
    'E': E,
    'C3': rotation_matrix('z', 120),
    'C3^2': rotation_matrix('z', 240),
    'C2': rotation_matrix('x', 180),
    'i': inversion(),
    'S6': improper_rotation_matrix('z', 60),
    'σd1': rotation_matrix('z', 0) @ reflection_matrix('xz'),
    'σd2': rotation_matrix('z', 120) @ reflection_matrix('xz') @ rotation_matrix('z', -120),
    'σd3': rotation_matrix('z', 240) @ reflection_matrix('xz') @ rotation_matrix('z', -240)
}

# Hexagonal point groups
point_groups['C6'] = {
    'E': E,
    'C6': rotation_matrix('z', 60),
    'C3': rotation_matrix('z', 120),
    'C2': rotation_matrix('z', 180),
    'C6^5': rotation_matrix('z', 300)
}

point_groups['C3h'] = {
    'E': E,
    'C3': rotation_matrix('z', 120),
    'C3^2': rotation_matrix('z', 240),
    'σh': reflection_matrix('xy'),
    'S3': improper_rotation_matrix('z', 120),
    'S3^2': improper_rotation_matrix('z', 240)
}

point_groups['C6h'] = {
    'E': E,
    'C6': rotation_matrix('z', 60),
    'C3': rotation_matrix('z', 120),
    'C2': rotation_matrix('z', 180),
    'C6^5': rotation_matrix('z', 300),
    'i': inversion(),
    'S3': improper_rotation_matrix('z', 120),
    'σh': reflection_matrix('xy')
}

point_groups['D6'] = {
    'E': E,
    'C6': rotation_matrix('z', 60),
    'C3': rotation_matrix('z', 120),
    'C2': rotation_matrix('z', 180),
    'C2(x)': rotation_matrix('x', 180),
    'C2(y)': rotation_matrix('y', 180)
}

point_groups['C6v'] = {
    'E': E,
    'C6': rotation_matrix('z', 60),
    'C3': rotation_matrix('z', 120),
    'C2': rotation_matrix('z', 180),
    'σv': reflection_matrix('xz'),
    'σv2': rotation_matrix('z', 60) @ reflection_matrix('xz') @ rotation_matrix('z', -60),
    'σv3': rotation_matrix('z', 120) @ reflection_matrix('xz') @ rotation_matrix('z', -120)
}

point_groups['D3h'] = {
    'E': E,
    'C3': rotation_matrix('z', 120),
    'C3^2': rotation_matrix('z', 240),
    'C2': rotation_matrix('x', 180),
    'σh': reflection_matrix('xy'),
    'σv1': reflection_matrix('xz'),
    'σv2': rotation_matrix('z', 120) @ reflection_matrix('xz') @ rotation_matrix('z', -120),
    'σv3': rotation_matrix('z', 240) @ reflection_matrix('xz') @ rotation_matrix('z', -240)
}

point_groups['D6h'] = {
    'E': E,
    'C6': rotation_matrix('z', 60),
    'C3': rotation_matrix('z', 120),
    'C2': rotation_matrix('z', 180),
    'C2(x)': rotation_matrix('x', 180),
    'i': inversion(),
    'σh': reflection_matrix('xy'),
    'σd1': reflection_matrix('xz'),
    'σd2': rotation_matrix('z', 60) @ reflection_matrix('xz') @ rotation_matrix('z', -60)
}

# Cubic (すでに Td, Oh 追加済)
point_groups['T'] = {
    'E': E,
    'C3': rotation_matrix('z', 120),
    'C2': rotation_matrix('x', 180)
}

point_groups['Th'] = {
    'E': E,
    'C3': rotation_matrix('z', 120),
    'C2': rotation_matrix('x', 180),
    'i': inversion()
}

point_groups['O'] = {
    'E': E,
    'C4(z)': rotation_matrix('z', 90),
    'C3(1)': rotation_matrix('z', 120),
    'C3(2)': rotation_matrix('z', 240),
    'C2(x)': rotation_matrix('x', 180),
    'C2(y)': rotation_matrix('y', 180),
    'C2(z)': rotation_matrix('z', 180)
}


# 完成した点群の個数を確認
len(point_groups)  # Should be 32 now


