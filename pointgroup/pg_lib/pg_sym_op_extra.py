

"""
結晶系ごとに点群をまとめて返す補助関数群
三斜晶, 単斜晶, 斜方晶, 正方晶, 三方晶, 六方晶, 立方晶
"""
from .rep_sym_op import ( 
    identity, 
    inversion, 
    rotation_matrix, 
    improper_rotation_matrix, 
    reflection_matrix,
)

import numpy as np

def get_rep(cell=None):
    """
    全結晶系の点群操作辞書をまとめて返す
    cell: 格子情報（必要に応じて渡す）
    return: dict {点群名: {操作名: 3x3行列}}
    """
    all_pg = {}
    for f in [make_triclinic, make_monoclinic, make_orthorhombic, make_tetragonal, make_trigonal, make_hexagonal, make_cubic]:
        all_pg.update(f(cell=cell))
    return all_pg

############################################################
############################################################

def make_triclinic(cell=None):
    """
    三斜晶系の点群（C1, Ci）の対称操作をO(3)行列で返す
    cell: 格子情報（未使用でもOK）
    return: dict {点群名: {操作名: 3x3行列}}
    """
    point_groups = {}
    # C1 (identity only)
    point_groups['C1'] = {'E  E': identity()}
    # Ci (identity, inversion)
    point_groups['Ci'] = {
        'E  E': identity(),
        'i  I': inversion(cell=cell)
    }
    return point_groups

def make_monoclinic(cell=None):
    """
    単斜晶系の点群（C2, Cs, C2h）の対称操作をO(3)行列で返す
    cell: 格子情報（必要に応じて渡す）
    return: dict {点群名: {操作名: 3x3行列}}
    """
    point_groups = {}
    # C2 (2回回転)
    point_groups['C2'] = {
        'E  E': identity(),
        'C2  C2(y)': rotation_matrix('y', 180, cell=cell)
    }
    # Cs (鏡映)
    point_groups['Cs'] = {
        'E  E': identity(),
        'sigmah  sigma(xz)': reflection_matrix('xz', cell=cell)
    }
    # C2h (2回回転＋反転)
    point_groups['C2h'] = {
        'E  E': identity(),
        'C2  C2(y)': rotation_matrix('y', 180, cell=cell),
        'i  I': inversion(cell),
        'sigmah  sigmah(xy)': reflection_matrix('xy', cell=cell)
    }
    return point_groups

def make_orthorhombic(cell=None):
    """
    斜方晶系の点群（D2, C2v, D2h）の対称操作をO(3)行列で返す
    cell: 格子情報（必要に応じて渡す）
    return: dict {点群名: {操作名: 3x3行列}}
    """
    point_groups = {}
    # D2 (3つの2回回転)
    point_groups['D2'] = {
        'E  E': identity(),
        'C2(x)  C2(x)': rotation_matrix('x', 180, cell=cell),
        'C2(y)  C2(y)': rotation_matrix('y', 180, cell=cell),
        'C2(z)  C2(z)': rotation_matrix('z', 180, cell=cell)
    }
    # C2v (2回回転＋2鏡映)
    point_groups['C2v'] = {
        'E  E': identity(),
        'C2  C2(z)': rotation_matrix('z', 180, cell=cell),
        'sigmav(xz)  sigma(xz)': reflection_matrix('xz', cell=cell),
        'sigmav(yz)  sigma(yz)': reflection_matrix('yz', cell=cell)
    }
    # D2h (3つの2回回転＋反転＋3鏡映)
    point_groups['D2h'] = {
        'E  E': identity(),
        'C2(x)  C2(x)': rotation_matrix('x', 180, cell=cell),
        'C2(y)  C2(y)': rotation_matrix('y', 180, cell=cell),
        'C2(z)  C2(z)': rotation_matrix('z', 180, cell=cell),
        'i  I': inversion(cell),
        'sigmaxy  sigma(xy)': reflection_matrix('xy', cell=cell),
        'sigmayz  sigma(yz)': reflection_matrix('yz', cell=cell),
        'sigmaxz  sigma(xz)': reflection_matrix('xz', cell=cell)
    }
    return point_groups

def make_tetragonal(cell=None):
    """
    正方晶系の点群（C4, S4, C4h, D4, C4v, D2d, D4h）の対称操作をO(3)行列で返す
    cell: 格子情報（必要に応じて渡す）
    return: dict {点群名: {操作名: 3x3行列}}
    """
    point_groups = {}
    # C4 (4回回転)
    point_groups['C4'] = {
        'E  E': identity(),
        'C4  C4(z)': rotation_matrix('z', 90, cell=cell),
        'C2  C2(z)': rotation_matrix('z', 180, cell=cell),
        'C4_3  C4^3(z)': rotation_matrix('z', 270, cell=cell)
    }
    # S4 (4回不正回転)
    point_groups['S4'] = {
        'E  E': identity(),
        'S4  S4': improper_rotation_matrix('z', 90, cell=cell),
        'C2  C2(z)': rotation_matrix('z', 180, cell=cell),
        'S4_3  S4^3': improper_rotation_matrix('z', 270, cell=cell)
    }
    # C4h (4回回転＋水平鏡映)
    point_groups['C4h'] = {
        'E  E': identity(),
        'C4  C4(z)': rotation_matrix('z', 90, cell=cell),
        'C2  C2(z)': rotation_matrix('z', 180, cell=cell),
        'C4_3  C4^3(z)': rotation_matrix('z', 270, cell=cell),
        'i  I': inversion(cell),
        'S4  S4': improper_rotation_matrix('z', 90, cell=cell),
        'S4_3  S4^3': improper_rotation_matrix('z', 270, cell=cell),
        'sigmah  sigmah(xy)': reflection_matrix('xy', cell=cell)
    }
    # D4 (4回回転＋2回回転×2＋2回回転×2)
    # d1: x=y, d2: x=-y
    sqrt2_inv = 1/np.sqrt(2)
    d1 = [sqrt2_inv, sqrt2_inv, 0]
    d2 = [sqrt2_inv, -sqrt2_inv, 0]
    point_groups['D4'] = {
        'E  E': identity(),
        'C4  C4(z)': rotation_matrix('z', 90, cell=cell),
        'C2  C2(z)': rotation_matrix('z', 180, cell=cell),
        'C4_3  C4^3(z)': rotation_matrix('z', 270, cell=cell),
        'C2x  C2(x)': rotation_matrix('x', 180, cell=cell),
        'C2y  C2(y)': rotation_matrix('y', 180, cell=cell),
        'C2d1  C2(d1)': rotation_matrix(d1, 180, cell=cell),
        'C2d2  C2(d2)': rotation_matrix(d2, 180, cell=cell)
    }
    # C4v (4回回転＋2回回転＋2鏡映×2)
    # sigmad1: x=y, sigmad2: x=-y
    point_groups['C4v'] = {
        'E  E': identity(),
        'C4  C4(z)': rotation_matrix('z', 90, cell),
        'C2  C2(z)': rotation_matrix('z', 180, cell),
        'C4_3  C4^3(z)': rotation_matrix('z', 270, cell),
        'sigmavxz  sigmav(xz)': reflection_matrix('xz', cell),
        'sigmavyz  sigmav(yz)': reflection_matrix('yz', cell),
        'sigmad1  sigmad1': reflection_matrix([sqrt2_inv, sqrt2_inv, 0], cell),
        'sigmad2  sigmad2': reflection_matrix([sqrt2_inv, -sqrt2_inv, 0], cell)
    }
    # D2d (2回回転×3＋S4×2＋鏡映×2)
    point_groups['D2d'] = {
        'E  E': identity(),
        'C2z  C2(z)': rotation_matrix('z', 180, cell),
        'C2x  C2(x)': rotation_matrix('x', 180, cell),
        'C2y  C2(y)': rotation_matrix('y', 180, cell),
        'S4  S4': improper_rotation_matrix('z', 90, cell),
        'S4_3  S4^3': improper_rotation_matrix('z', 270, cell),
        'sigmad1  sigmad1': reflection_matrix([sqrt2_inv, sqrt2_inv, 0], cell),
        'sigmad2  sigmad2': reflection_matrix([sqrt2_inv, -sqrt2_inv, 0], cell)
    }
    # D4h (4回回転＋2回回転×2＋反転＋S4×2＋鏡映×3)
    point_groups['D4h'] = {
        'E  E': identity(),
        'C4  C4(z)': rotation_matrix('z', 90, cell),
        'C2  C2(z)': rotation_matrix('z', 180, cell),
        'C4_3  C4^3(z)': rotation_matrix('z', 270, cell),
        'C2x  C2(x)': rotation_matrix('x', 180, cell),
        'C2y  C2(y)': rotation_matrix('y', 180, cell),
        'i  I': inversion(cell),
        'S4  S4': improper_rotation_matrix('z', 90, cell),
        'S4_3  S4^3': improper_rotation_matrix('z', 270, cell),
        'sigmah  sigmah(xy)': reflection_matrix('xy', cell),
        'sigmavxz  sigmav(xz)': reflection_matrix('xz', cell),
        'sigmavyz  sigmav(yz)': reflection_matrix('yz', cell)
    }
    return point_groups

def make_trigonal(cell=None):
    """
    三方晶系の点群（C3, S6, D3, C3v, D3d）の対称操作をO(3)行列で返す
    cell: 格子情報（必要に応じて渡す）
    return: dict {点群名: {操作名: 3x3行列}}
    """
    point_groups = {}
    # C3 (3回回転)
    point_groups['C3'] = {
        'E  E': identity(),
        'C3  C3(z)': rotation_matrix('z', 120, cell=cell),
        'C3_2  C3^2(z)': rotation_matrix('z', 240, cell=cell)
    }
    # S6 (6回不正回転)
    point_groups['S6'] = {
        'E  E': identity(),
        'S6  S6': improper_rotation_matrix('z', 60, cell=cell),
        'S6_5  S6^5': improper_rotation_matrix('z', 300, cell=cell),
        'C3  C3(z)': rotation_matrix('z', 120, cell=cell),
        'C3_2  C3^2(z)': rotation_matrix('z', 240, cell=cell),
        'i  I': inversion(cell)
    }
    # D3 (3回回転＋2回回転×3)
    # C2(2): y軸120度回転したx軸, C2(3): y軸240度回転したx軸
    axis2 = np.dot(rotation_matrix('z', 120, cell=None), [1, 0, 0])
    axis3 = np.dot(rotation_matrix('z', 240, cell=None), [1, 0, 0])
    point_groups['D3'] = {
        'E  E': identity(),
        'C3  C3(z)': rotation_matrix('z', 120, cell=cell),
        'C3_2  C3^2(z)': rotation_matrix('z', 240, cell=cell),
        'C2_1  C2(1)': rotation_matrix('x', 180, cell=cell),
        'C2_2  C2(2)': rotation_matrix(axis2, 180, cell=cell),
        'C2_3  C2(3)': rotation_matrix(axis3, 180, cell=cell)
    }
    # C3v (3回回転＋鏡映×3)
    # sigmav2, sigmav3: z軸を含み120度, 240度回転した面
    sigmav2_axis = np.dot(rotation_matrix('z', 120, cell=None), [1,0,0])
    sigmav3_axis = np.dot(rotation_matrix('z', 240, cell=None), [1,0,0])
    point_groups['C3v'] = {
        'E  E': identity(),
        'C3  C3(z)': rotation_matrix('z', 120, cell=cell),
        'C3_2  C3^2(z)': rotation_matrix('z', 240, cell=cell),
        'sigmav1  sigmav1': reflection_matrix('xz', cell=cell),
        'sigmav2  sigmav2': reflection_matrix(sigmav2_axis, cell=cell),
        'sigmav3  sigmav3': reflection_matrix(sigmav3_axis, cell=cell)
    }
    # D3d (3回回転＋2回回転×3＋反転＋鏡映×3)
    # sigmad1: xy面, sigmad2,3: xy面をz軸中心に120,240度回転
    sigmad2_axis = np.dot(rotation_matrix('z', 120, cell=None), [0,0,1])
    sigmad3_axis = np.dot(rotation_matrix('z', 240, cell=None), [0,0,1])
    point_groups['D3d'] = {
        'E  E': identity(),
        'C3  C3(z)': rotation_matrix('z', 120, cell=cell),
        'C3_2  C3^2(z)': rotation_matrix('z', 240, cell=cell),
        'C2  C2': rotation_matrix('x', 180, cell=cell),
        'i  I': inversion(cell),
        'S6  S6': improper_rotation_matrix('z', 60, cell=cell),
        'sigmad1  sigmad1': reflection_matrix('xy', cell=cell),
        'sigmad2  sigmad2': reflection_matrix(sigmad2_axis, cell=cell),
        'sigmad3  sigmad3': reflection_matrix(sigmad3_axis, cell=cell)
    }
    return point_groups

def make_hexagonal(cell=None):
    """
    六方晶系の点群（C6, C3h, C6h, D6, C6v, D3h, D6h）の対称操作をO(3)行列で返す
    cell: 格子情報（必要に応じて渡す）
    return: dict {点群名: {操作名: 3x3行列}}
    """
    point_groups = {}
    # C6 (6回回転)
    point_groups['C6'] = {
        'E  E': identity(),
        'C6  C6': rotation_matrix('z', 60, cell),
        'C3  C3': rotation_matrix('z', 120, cell),
        'C2  C2': rotation_matrix('z', 180, cell),
        'C6_5  C6^5': rotation_matrix('z', 300, cell)
    }
    # C3h (3回回転＋水平鏡映)
    point_groups['C3h'] = {
        'E  E': identity(),
        'C3  C3': rotation_matrix('z', 120, cell),
        'C3_2  C3^2': rotation_matrix('z', 240, cell),
        'sigmah  sigmah': reflection_matrix('xy', cell),
        'S3  S3': improper_rotation_matrix('z', 120, cell),
        'S3_2  S3^2': improper_rotation_matrix('z', 240, cell)
    }
    # C6h (6回回転＋水平鏡映＋反転)
    point_groups['C6h'] = {
        'E  E': identity(),
        'C6  C6': rotation_matrix('z', 60, cell),
        'C3  C3': rotation_matrix('z', 120, cell),
        'C2  C2': rotation_matrix('z', 180, cell),
        'C6_5  C6^5': rotation_matrix('z', 300, cell),
        'i  I': inversion(cell),
        'S3  S3': improper_rotation_matrix('z', 120, cell),
        'sigmah  sigmah': reflection_matrix('xy', cell)
    }
    # D6 (6回回転＋2回回転×3)
    point_groups['D6'] = {
        'E  E': identity(),
        'C6  C6': rotation_matrix('z', 60, cell),
        'C3  C3': rotation_matrix('z', 120, cell),
        'C2  C2': rotation_matrix('z', 180, cell),
        'C2x  C2(x)': rotation_matrix('x', 180, cell),
        'C2y  C2(y)': rotation_matrix('y', 180, cell)
    }
    # C6v (6回回転＋鏡映×6)
    # sigmav2, sigmav3: xz面をz軸中心に120,240度回転
    sigmav2_axis = np.dot(rotation_matrix('z', 120, cell=None), [1,0,0])
    sigmav3_axis = np.dot(rotation_matrix('z', 240, cell=None), [1,0,0])
    point_groups['C6v'] = {
        'E  E': identity(),
        'C6  C6': rotation_matrix('z', 60, cell),
        'C3  C3': rotation_matrix('z', 120, cell),
        'C2  C2': rotation_matrix('z', 180, cell),
        'sigmav  sigmav': reflection_matrix('xz', cell),
        'sigmav2  sigmav2': reflection_matrix(sigmav2_axis, cell),
        'sigmav3  sigmav3': reflection_matrix(sigmav3_axis, cell)
    }
    # D3h (3回回転＋2回回転×3＋水平鏡映＋鏡映×3)
    # sigmav2, sigmav3: xz面をz軸中心に120,240度回転
    sigmav2_axis = np.dot(rotation_matrix('z', 120, cell=None), [1,0,0])
    sigmav3_axis = np.dot(rotation_matrix('z', 240, cell=None), [1,0,0])
    point_groups['D3h'] = {
        'E  E': identity(),
        'C3  C3': rotation_matrix('z', 120, cell),
        'C3_2  C3^2': rotation_matrix('z', 240, cell),
        'C2  C2': rotation_matrix('x', 180, cell),
        'sigmah  sigmah': reflection_matrix('xy', cell),
        'sigmav1  sigmav1': reflection_matrix('xz', cell),
        'sigmav2  sigmav2': reflection_matrix(sigmav2_axis, cell),
        'sigmav3  sigmav3': reflection_matrix(sigmav3_axis, cell)
    }
    # D6h (6回回転＋2回回転×3＋反転＋鏡映×6)
    # sigmad2: xz面をz軸中心に60度回転
    sigmad2_axis = np.dot(rotation_matrix('z', 60, cell=None), [1,0,0])
    point_groups['D6h'] = {
        'E  E': identity(),
        'C6  C6': rotation_matrix('z', 60, cell),
        'C3  C3': rotation_matrix('z', 120, cell),
        'C2  C2': rotation_matrix('z', 180, cell),
        'C2x  C2(x)': rotation_matrix('x', 180, cell),
        'i  I': inversion(cell),
        'sigmah  sigmah': reflection_matrix('xy', cell),
        'sigmad1  sigmad1': reflection_matrix('xz', cell),
        'sigmad2  sigmad2': reflection_matrix(sigmad2_axis, cell)
    }
    return point_groups

def make_cubic(cell=None):
    """
    立方晶系の点群（T, Th, O, Td, Oh）の対称操作をO(3)行列で返す
    cell: 格子情報（必要に応じて渡す）
    return: dict {点群名: {操作名: 3x3行列}}
    """
    point_groups = {}
    # T (三重回転＋2回回転)
    point_groups['T'] = {
        'E  E': identity(),
        'C3  C3': rotation_matrix('z', 120, cell),
        'C2  C2': rotation_matrix('x', 180, cell)
    }
    # Th (T＋反転)
    point_groups['Th'] = {
        'E  E': identity(),
        'C3  C3': rotation_matrix('z', 120, cell),
        'C2  C2': rotation_matrix('x', 180, cell),
        'i  I': inversion(cell)
    }
    # O (四重回転＋三重回転＋2回回転)
    point_groups['O'] = {
        'E  E': identity(),
        'C4  C4(z)': rotation_matrix('z', 90, cell),
        'C3_1  C3(1)': rotation_matrix('z', 120, cell),
        'C3_2  C3(2)': rotation_matrix('z', 240, cell),
        'C2x  C2(x)': rotation_matrix('x', 180, cell),
        'C2y  C2(y)': rotation_matrix('y', 180, cell),
        'C2z  C2(z)': rotation_matrix('z', 180, cell)
    }
    # Td (T＋鏡映)
    point_groups['Td'] = {
        'E  E': identity(),
        'C3  C3': rotation_matrix('z', 120, cell),
        'C2  C2': rotation_matrix('x', 180, cell),
        'S4  S4': rotation_matrix('z', 90, cell),
        'sigmad  sigmad': reflection_matrix('xz', cell)
    }
    # Oh (O＋反転＋鏡映)
    point_groups['Oh'] = {
        'E  E': identity(),
        'C4  C4(z)': rotation_matrix('z', 90, cell),
        'C3  C3': rotation_matrix('z', 120, cell),
        'C2x  C2(x)': rotation_matrix('x', 180, cell),
        'i  I': inversion(cell)
    }
    return point_groups
