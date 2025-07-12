"""
結晶系ごとに点群をまとめて返す補助関数群
三斜晶, 単斜晶, 斜方晶, 正方晶, 三方晶, 六方晶, 立方晶
"""
from .rep_sym_op import *

def make_triclinic(cell=None):
    """
    三斜晶系の点群（C1, Ci）の対称操作をO(3)行列で返す
    cell: 格子情報（未使用でもOK）
    return: dict {点群名: {操作名: 3x3行列}}
    """
    point_groups = {}
    # C1 (identity only)
    point_groups['C1'] = {'E': identity()}
    # Ci (identity, inversion)
    point_groups['Ci'] = {
        'E': identity(),
        'i': inversion(cell=cell)
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
        'E': identity(),
        'C2(y)': rotation_matrix('y', 180, cell=cell)
    }
    # Cs (鏡映)
    point_groups['Cs'] = {
        'E': identity(),
        'sigma(xz)': reflection_matrix('xz', cell=cell)
    }
    # C2h (2回回転＋反転)
    point_groups['C2h'] = {
        'E': identity(),
        'C2(y)': rotation_matrix('y', 180, cell=cell),
        'i': inversion(cell),
        'sigmah(xy)': reflection_matrix('xy', cell=cell)
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
        'E': identity(),
        'C2(x)': rotation_matrix('x', 180, cell=cell),
        'C2(y)': rotation_matrix('y', 180, cell=cell),
        'C2(z)': rotation_matrix('z', 180, cell=cell)
    }
    # C2v (2回回転＋2鏡映)
    point_groups['C2v'] = {
        'E': identity(),
        'C2(z)': rotation_matrix('z', 180, cell=cell),
        'sigma(xz)': reflection_matrix('xz', cell=cell),
        'sigma(yz)': reflection_matrix('yz', cell=cell)
    }
    # D2h (3つの2回回転＋反転＋3鏡映)
    point_groups['D2h'] = {
        'E': identity(),
        'C2(x)': rotation_matrix('x', 180, cell=cell),
        'C2(y)': rotation_matrix('y', 180, cell=cell),
        'C2(z)': rotation_matrix('z', 180, cell=cell),
        'i': inversion(cell),
        'sigma(xy)': reflection_matrix('xy', cell=cell),
        'sigma(yz)': reflection_matrix('yz', cell=cell),
        'sigma(xz)': reflection_matrix('xz', cell=cell)
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
        'E': identity(),
        'C4(z)': rotation_matrix('z', 90, cell=cell),
        'C2(z)': rotation_matrix('z', 180, cell=cell),
        'C4^3(z)': rotation_matrix('z', 270, cell=cell)
    }
    # S4 (4回不正回転)
    point_groups['S4'] = {
        'E': identity(),
        'S4': improper_rotation_matrix('z', 90, cell=cell),
        'C2(z)': rotation_matrix('z', 180, cell=cell),
        'S4^3': improper_rotation_matrix('z', 270, cell=cell)
    }
    # C4h (4回回転＋水平鏡映)
    point_groups['C4h'] = {
        'E': identity(),
        'C4(z)': rotation_matrix('z', 90, cell=cell),
        'C2(z)': rotation_matrix('z', 180, cell=cell),
        'C4^3(z)': rotation_matrix('z', 270, cell=cell),
        'i': inversion(cell),
        'S4': improper_rotation_matrix('z', 90, cell=cell),
        'S4^3': improper_rotation_matrix('z', 270, cell=cell),
        'sigmah(xy)': reflection_matrix('xy', cell=cell)
    }
    # D4 (4回回転＋2回回転×2＋2回回転×2)
    point_groups['D4'] = {
        'E': identity(),
        'C4(z)': rotation_matrix('z', 90, cell=cell),
        'C2(z)': rotation_matrix('z', 180, cell=cell),
        'C4^3(z)': rotation_matrix('z', 270, cell=cell),
        'C2(x)': rotation_matrix('x', 180, cell=cell),
        'C2(y)': rotation_matrix('y', 180, cell=cell),
        'C2(d1)': gen_c2_d1(cell=cell),
        'C2(d2)': gen_c2_d2(cell=cell)
    }
    # C4v (4回回転＋2回回転＋2鏡映×2)
    point_groups['C4v'] = {
        'E': identity(),
        'C4(z)': rotation_matrix('z', 90, cell),
        'C2(z)': rotation_matrix('z', 180, cell),
        'C4^3(z)': rotation_matrix('z', 270, cell),
        'sigmav(xz)': reflection_matrix('xz', cell),
        'sigmav(yz)': reflection_matrix('yz', cell),
        'sigmad1': gen_sigmad1(cell),
        'sigmad2': gen_sigmad2(cell)
    }
    # D2d (2回回転×3＋S4×2＋鏡映×2)
    point_groups['D2d'] = {
        'E': identity(),
        'C2(z)': rotation_matrix('z', 180, cell),
        'C2(x)': rotation_matrix('x', 180, cell),
        'C2(y)': rotation_matrix('y', 180, cell),
        'S4': improper_rotation_matrix('z', 90, cell),
        'S4^3': improper_rotation_matrix('z', 270, cell),
        'sigmad1': gen_sigmad1(cell),
        'sigmad2': gen_sigmad2(cell)
    }
    # D4h (4回回転＋2回回転×2＋反転＋S4×2＋鏡映×3)
    point_groups['D4h'] = {
        'E': identity(),
        'C4(z)': rotation_matrix('z', 90, cell),
        'C2(z)': rotation_matrix('z', 180, cell),
        'C4^3(z)': rotation_matrix('z', 270, cell),
        'C2(x)': rotation_matrix('x', 180, cell),
        'C2(y)': rotation_matrix('y', 180, cell),
        'i': inversion(cell),
        'S4': improper_rotation_matrix('z', 90, cell),
        'S4^3': improper_rotation_matrix('z', 270, cell),
        'sigmah(xy)': reflection_matrix('xy', cell),
        'sigmav(xz)': reflection_matrix('xz', cell),
        'sigmav(yz)': reflection_matrix('yz', cell)
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
        'E': identity(),
        'C3(z)': rotation_matrix('z', 120, cell=cell),
        'C3^2(z)': rotation_matrix('z', 240, cell=cell)
    }
    # S6 (6回不正回転)
    point_groups['S6'] = {
        'E': identity(),
        'S6': improper_rotation_matrix('z', 60, cell=cell),
        'S6^5': improper_rotation_matrix('z', 300, cell=cell),
        'C3(z)': rotation_matrix('z', 120, cell=cell),
        'C3^2(z)': rotation_matrix('z', 240, cell=cell),
        'i': inversion(cell)
    }
    # D3 (3回回転＋2回回転×3)
    point_groups['D3'] = {
        'E': identity(),
        'C3(z)': rotation_matrix('z', 120, cell=cell),
        'C3^2(z)': rotation_matrix('z', 240, cell=cell),
        'C2(1)': rotation_matrix('x', 180, cell=cell),
        'C2(2)': gen_c2_2(cell=cell),
        'C2(3)': gen_c2_3(cell=cell)
    }
    # C3v (3回回転＋鏡映×3)
    point_groups['C3v'] = {
        'E': identity(),
        'C3(z)': rotation_matrix('z', 120, cell=cell),
        'C3^2(z)': rotation_matrix('z', 240, cell=cell),
        'sigmav1': reflection_matrix('xz', cell=cell),
        'sigmav2': gen_sigmav2(cell=cell),
        'sigmav3': gen_sigmav3(cell=cell)
    }
    # D3d (3回回転＋2回回転×3＋反転＋鏡映×3)
    point_groups['D3d'] = {
        'E': identity(),
        'C3(z)': rotation_matrix('z', 120, cell=cell),
        'C3^2(z)': rotation_matrix('z', 240, cell=cell),
        'C2': rotation_matrix('x', 180, cell=cell),
        'i': inversion(cell),
        'S6': improper_rotation_matrix('z', 60, cell=cell),
        'sigmad1': gen_sigmad1_2(cell=cell),
        'sigmad2': gen_sigmad2_2(cell=cell),
        'sigmad3': gen_sigmad3(cell=cell)
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
        'E': identity(),
        'C6': rotation_matrix('z', 60, cell),
        'C3': rotation_matrix('z', 120, cell),
        'C2': rotation_matrix('z', 180, cell),
        'C6^5': rotation_matrix('z', 300, cell)
    }
    # C3h (3回回転＋水平鏡映)
    point_groups['C3h'] = {
        'E': identity(),
        'C3': rotation_matrix('z', 120, cell),
        'C3^2': rotation_matrix('z', 240, cell),
        'sigmah': reflection_matrix('xy', cell),
        'S3': improper_rotation_matrix('z', 120, cell),
        'S3^2': improper_rotation_matrix('z', 240, cell)
    }
    # C6h (6回回転＋水平鏡映＋反転)
    point_groups['C6h'] = {
        'E': identity(),
        'C6': rotation_matrix('z', 60, cell),
        'C3': rotation_matrix('z', 120, cell),
        'C2': rotation_matrix('z', 180, cell),
        'C6^5': rotation_matrix('z', 300, cell),
        'i': inversion(cell),
        'S3': improper_rotation_matrix('z', 120, cell),
        'sigmah': reflection_matrix('xy', cell)
    }
    # D6 (6回回転＋2回回転×3)
    point_groups['D6'] = {
        'E': identity(),
        'C6': rotation_matrix('z', 60, cell),
        'C3': rotation_matrix('z', 120, cell),
        'C2': rotation_matrix('z', 180, cell),
        'C2(x)': rotation_matrix('x', 180, cell),
        'C2(y)': rotation_matrix('y', 180, cell)
    }
    # C6v (6回回転＋鏡映×6)
    point_groups['C6v'] = {
        'E': identity(),
        'C6': rotation_matrix('z', 60, cell),
        'C3': rotation_matrix('z', 120, cell),
        'C2': rotation_matrix('z', 180, cell),
        'sigmav': reflection_matrix('xz', cell),
        'sigmav2': gen_sigmav2_2(cell),
        'sigmav3': gen_sigmav3_2(cell)
    }
    # D3h (3回回転＋2回回転×3＋水平鏡映＋鏡映×3)
    point_groups['D3h'] = {
        'E': identity(),
        'C3': rotation_matrix('z', 120, cell),
        'C3^2': rotation_matrix('z', 240, cell),
        'C2': rotation_matrix('x', 180, cell),
        'sigmah': reflection_matrix('xy', cell),
        'sigmav1': reflection_matrix('xz', cell),
        'sigmav2': gen_sigmav2(cell),
        'sigmav3': gen_sigmav3(cell)
    }
    # D6h (6回回転＋2回回転×3＋反転＋鏡映×6)
    point_groups['D6h'] = {
        'E': identity(),
        'C6': rotation_matrix('z', 60, cell),
        'C3': rotation_matrix('z', 120, cell),
        'C2': rotation_matrix('z', 180, cell),
        'C2(x)': rotation_matrix('x', 180, cell),
        'i': inversion(cell),
        'sigmah': reflection_matrix('xy', cell),
        'sigmad1': reflection_matrix('xz', cell),
        'sigmad2': gen_sigmad2_3(cell)
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
        'E': identity(),
        'C3': rotation_matrix('z', 120, cell),
        'C2': rotation_matrix('x', 180, cell)
    }
    # Th (T＋反転)
    point_groups['Th'] = {
        'E': identity(),
        'C3': rotation_matrix('z', 120, cell),
        'C2': rotation_matrix('x', 180, cell),
        'i': inversion(cell)
    }
    # O (四重回転＋三重回転＋2回回転)
    point_groups['O'] = {
        'E': identity(),
        'C4(z)': rotation_matrix('z', 90, cell),
        'C3(1)': rotation_matrix('z', 120, cell),
        'C3(2)': rotation_matrix('z', 240, cell),
        'C2(x)': rotation_matrix('x', 180, cell),
        'C2(y)': rotation_matrix('y', 180, cell),
        'C2(z)': rotation_matrix('z', 180, cell)
    }
    # Td (T＋鏡映)
    point_groups['Td'] = {
        'E': identity(),
        'C3': rotation_matrix('z', 120, cell),
        'C2': rotation_matrix('x', 180, cell),
        'S4': rotation_matrix('z', 90, cell),
        'sigmad': reflection_matrix('xz', cell)
    }
    # Oh (O＋反転＋鏡映)
    point_groups['Oh'] = {
        'E': identity(),
        'C4(z)': rotation_matrix('z', 90, cell),
        'C3': rotation_matrix('z', 120, cell),
        'C2(x)': rotation_matrix('x', 180, cell),
        'i': inversion(cell)
    }
    return point_groups
