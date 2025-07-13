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
    rotation_matrix_axis,
    reflection_matrix_axis,
    improper_rotation_matrix_axis,
)

import numpy as np


def get_rep(cell=None):
    """
    全結晶系の点群操作辞書をまとめて返す
    cell: 格子情報（必要に応じて渡す）
    return: dict {点群名: {操作名: 3x3行列}}
    """
    all_pg = {}
    for f in [
        make_triclinic,
        make_monoclinic,
        make_orthorhombic,
        make_tetragonal,
        make_trigonal,
        make_hexagonal,
        make_cubic,
    ]:
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
    point_groups["C1"] = {"E  E": identity()}
    # Ci (identity, inversion)
    point_groups["Ci"] = {"E  E": identity(), "i  I": inversion(cell=cell)}
    return point_groups


def make_monoclinic(cell=None):
    """
    単斜晶系の点群（C2, Cs, C2h）の対称操作をO(3)行列で返す
    cell: 格子情報（必要に応じて渡す）
    return: dict {点群名: {操作名: 3x3行列}}
    """
    point_groups = {}
    # C2 (2回回転)
    point_groups["C2"] = {
        "E  E": identity(),
        "C2  C2(y)": rotation_matrix("y", 180, cell=cell),
    }
    # Cs (鏡映)
    point_groups["Cs"] = {
        "E  E": identity(),
        "sigmah  sigma(xz)": reflection_matrix("xz", cell=cell),
    }
    # C2h (2回回転＋反転)
    point_groups["C2h"] = {
        "E  E": identity(),
        "C2  C2(y)": rotation_matrix("y", 180, cell=cell),
        "i  I": inversion(cell),
        "sigmah  sigmah(xy)": reflection_matrix("xy", cell=cell),
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
    point_groups["D2"] = {
        "E  E": identity(),
        "C2(z)  C2(z)": rotation_matrix("z", 180, cell=cell),
        "C2(y)  C2(y)": rotation_matrix("y", 180, cell=cell),
        "C2(x)  C2(x)": rotation_matrix("x", 180, cell=cell),
    }
    # C2v (2回回転＋2鏡映)
    point_groups["C2v"] = {
        "E  E": identity(),
        "C2  C2(z)": rotation_matrix("z", 180, cell=cell),
        "sigmav(xz)  sigma(xz)": reflection_matrix("xz", cell=cell),
        "sigmav(yz)  sigma(yz)": reflection_matrix("yz", cell=cell),
    }
    # D2h (3つの2回回転＋反転＋3鏡映)
    point_groups["D2h"] = {
        "E  E": identity(),
        "C2(z)  C2(z)": rotation_matrix("z", 180, cell=cell),
        "C2(y)  C2(y)": rotation_matrix("y", 180, cell=cell),
        "C2(x)  C2(x)": rotation_matrix("x", 180, cell=cell),
        "i  I": inversion(cell),
        "sigmaxy  sigma(xy)": reflection_matrix("xy", cell=cell),
        "sigmaxz  sigma(xz)": reflection_matrix("xz", cell=cell),
        "sigmayz  sigma(yz)": reflection_matrix("yz", cell=cell),
    }
    return point_groups


def make_tetragonal(cell=None):
    """
    正方晶系の点群（C4, S4, C4h, D4, C4v, D2d, D4h）の対称操作をO(3)行列で返す
    cell: 格子情報（必要に応じて渡す）
    return: dict {点群名: {操作名: 3x3行列}}
    """
    point_groups = {}

    # 対角線方向のベクトル定義
    sqrt2_inv = 1 / np.sqrt(2)
    d1 = [sqrt2_inv, sqrt2_inv, 0]  # x=y方向
    d2 = [sqrt2_inv, -sqrt2_inv, 0]  # x=-y方向

    # C4 (4回回転) - 4個の操作
    point_groups["C4"] = {
        "E  E": identity(),
        "C4  C4(z)": rotation_matrix("z", 90, cell=cell),
        "C2  C2(z)": rotation_matrix("z", 180, cell=cell),
        "C4_3  C4^3(z)": rotation_matrix("z", 270, cell=cell),
    }

    # S4 (4回不正回転) - 4個の操作
    point_groups["S4"] = {
        "E  E": identity(),
        "S4  S4": improper_rotation_matrix("z", 90, cell=cell),
        "C2  C2(z)": rotation_matrix("z", 180, cell=cell),
        "S4_3  S4^3": improper_rotation_matrix("z", 270, cell=cell),
    }

    # C4h (4回回転＋水平鏡映) - 8個の操作
    point_groups["C4h"] = {
        "E  E": identity(),
        "C4  C4(z)": rotation_matrix("z", 90, cell=cell),
        "C2  C2(z)": rotation_matrix("z", 180, cell=cell),
        "C4_3  C4^3(z)": rotation_matrix("z", 270, cell=cell),
        "i  I": inversion(cell),
        "S4  S4": improper_rotation_matrix("z", 90, cell=cell),
        "sigmah  sigmah(xy)": reflection_matrix("xy", cell=cell),
        "S4_3  S4^3": improper_rotation_matrix("z", 270, cell=cell),
    }

    # D4 (4回回転＋2回回転×4) - 8個の操作
    point_groups["D4"] = {
        "E  E": identity(),
        "2C4  C4(z)": rotation_matrix("z", 90, cell=cell),
        "C2  C2(z)": rotation_matrix("z", 180, cell=cell),
        "2C4_3  C4^3(z)": rotation_matrix("z", 270, cell=cell),
        "2C2_1  C2(x)": rotation_matrix("x", 180, cell=cell),
        "2C2_2  C2(y)": rotation_matrix("y", 180, cell=cell),
        "C2d1  C2(d1)": rotation_matrix_axis(d1, 180, cell=cell),
        "C2d2  C2(d2)": rotation_matrix_axis(d2, 180, cell=cell),
    }

    # C4v (4回回転＋鏡映×4) - 8個の操作
    point_groups["C4v"] = {
        "E  E": identity(),
        "2C4  C4(z)": rotation_matrix("z", 90, cell=cell),
        "C2  C2(z)": rotation_matrix("z", 180, cell=cell),
        "2C4_3  C4^3(z)": rotation_matrix("z", 270, cell=cell),
        "2sigmav  sigmav(xz)": reflection_matrix("xz", cell=cell),
        "2sigmav  sigmav(yz)": reflection_matrix("yz", cell=cell),
        "2sigmad  sigmad1": reflection_matrix_axis(d1, cell=cell),
        "2sigmad  sigmad2": reflection_matrix_axis(d2, cell=cell),
    }

    # D2d (2回回転×3＋S4×2＋鏡映×2) - 8個の操作
    point_groups["D2d"] = {
        "E  E": identity(),
        "C2z  C2(z)": rotation_matrix("z", 180, cell=cell),
        "C2x  C2(x)": rotation_matrix("x", 180, cell=cell),
        "C2y  C2(y)": rotation_matrix("y", 180, cell=cell),
        "S4  S4": improper_rotation_matrix("z", 90, cell=cell),
        "S4_3  S4^3": improper_rotation_matrix("z", 270, cell=cell),
        "sigmad1  sigmad1": reflection_matrix_axis(d1, cell=cell),
        "sigmad2  sigmad2": reflection_matrix_axis(d2, cell=cell),
    }

    # D4h (4回回転＋2回回転×4＋反転＋S4×2＋鏡映×5) - 16個の操作
    point_groups["D4h"] = {
        "E  E": identity(),
        "C4  C4(z)": rotation_matrix("z", 90, cell=cell),
        "C2  C2(z)": rotation_matrix("z", 180, cell=cell),
        "C4_3  C4^3(z)": rotation_matrix("z", 270, cell=cell),
        "C2x  C2(x)": rotation_matrix("x", 180, cell=cell),
        "C2y  C2(y)": rotation_matrix("y", 180, cell=cell),
        "C2d1  C2(d1)": rotation_matrix_axis(d1, 180, cell=cell),
        "C2d2  C2(d2)": rotation_matrix_axis(d2, 180, cell=cell),
        "i  I": inversion(cell),
        "S4  S4": improper_rotation_matrix("z", 90, cell=cell),
        "S4_3  S4^3": improper_rotation_matrix("z", 270, cell=cell),
        "sigmah  sigmah(xy)": reflection_matrix("xy", cell=cell),
        "sigmavxz  sigmav(xz)": reflection_matrix("xz", cell=cell),
        "sigmavyz  sigmav(yz)": reflection_matrix("yz", cell=cell),
        "sigmad1  sigmad1": reflection_matrix_axis(d1, cell=cell),
        "sigmad2  sigmad2": reflection_matrix_axis(d2, cell=cell),
    }

    # より一般的な正方晶系の点群も追加
    # C2v (2回回転＋鏡映×2) - 正方格子でも使用される
    point_groups["C2v"] = {
        "E  E": identity(),
        "C2  C2(z)": rotation_matrix("z", 180, cell=cell),
        "sigmavxz  sigmav(xz)": reflection_matrix("xz", cell=cell),
        "sigmavyz  sigmav(yz)": reflection_matrix("yz", cell=cell),
    }

    # D2 (2回回転×3) - 正方格子でも使用される
    point_groups["D2"] = {
        "E  E": identity(),
        "C2x  C2(x)": rotation_matrix("x", 180, cell=cell),
        "C2y  C2(y)": rotation_matrix("y", 180, cell=cell),
        "C2z  C2(z)": rotation_matrix("z", 180, cell=cell),
    }

    # D2h (2回回転×3＋反転＋鏡映×3) - 正方格子でも使用される
    point_groups["D2h"] = {
        "E  E": identity(),
        "C2x  C2(x)": rotation_matrix("x", 180, cell=cell),
        "C2y  C2(y)": rotation_matrix("y", 180, cell=cell),
        "C2z  C2(z)": rotation_matrix("z", 180, cell=cell),
        "i  I": inversion(cell),
        "sigmaxy  sigma(xy)": reflection_matrix("xy", cell=cell),
        "sigmaxz  sigma(xz)": reflection_matrix("xz", cell=cell),
        "sigmayz  sigma(yz)": reflection_matrix("yz", cell=cell),
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
    point_groups["C3"] = {
        "E  E": identity(),
        "C3  C3(z)": rotation_matrix("z", 120, cell=cell),
        "C3_2  C3^2(z)": rotation_matrix("z", 240, cell=cell),
    }
    # S6 (6回不正回転)
    point_groups["S6"] = {
        "E  E": identity(),
        "S6  S6": improper_rotation_matrix("z", 60, cell=cell),
        "S6_5  S6^5": improper_rotation_matrix("z", 300, cell=cell),
        "C3  C3(z)": rotation_matrix("z", 120, cell=cell),
        "C3_2  C3^2(z)": rotation_matrix("z", 240, cell=cell),
        "i  I": inversion(cell),
    }
    # D3 (3回回転＋2回回転×3)
    # C2(2): y軸120度回転したx軸, C2(3): y軸240度回転したx軸
    axis2 = np.dot(rotation_matrix("z", 120, cell=None), [1, 0, 0])
    axis3 = np.dot(rotation_matrix("z", 240, cell=None), [1, 0, 0])
    point_groups["D3"] = {
        "E  E": identity(),
        "2C3  C3(z)": rotation_matrix("z", 120, cell=cell),
        "2C3  C3^2(z)": rotation_matrix("z", 240, cell=cell),
        "3C2  C2(1)": rotation_matrix("x", 180, cell=cell),
        "3C2  C2(2)": rotation_matrix_axis(axis2, 180, cell=cell),
        "3C2  C2(3)": rotation_matrix_axis(axis3, 180, cell=cell),
    }
    # C3v (3回回転＋鏡映×3)
    # sigmav2, sigmav3: z軸を含み120度, 240度回転した面
    sigmav2_axis = np.dot(rotation_matrix("z", 120, cell=None), [1, 0, 0])
    sigmav3_axis = np.dot(rotation_matrix("z", 240, cell=None), [1, 0, 0])
    point_groups["C3v"] = {
        "E  E": identity(),
        "2C3  C3(z)": rotation_matrix("z", 120, cell=cell),
        "2C3  C3^2(z)": rotation_matrix("z", 240, cell=cell),
        "3sigmav  sigmav1": reflection_matrix("xz", cell=cell),
        "3sigmav  sigmav2": reflection_matrix_axis(sigmav2_axis, cell=cell),
        "3sigmav  sigmav3": reflection_matrix_axis(sigmav3_axis, cell=cell),
    }
    # D3d (3回回転＋2回回転×3＋反転＋鏡映×3)
    # sigmad1: xy面, sigmad2,3: xy面をz軸中心に120,240度回転
    sigmad2_axis = np.dot(rotation_matrix("z", 120, cell=None), [0, 0, 1])
    sigmad3_axis = np.dot(rotation_matrix("z", 240, cell=None), [0, 0, 1])
    point_groups["D3d"] = {
        "E  E": identity(),
        "2C3  C3(z)": rotation_matrix("z", 120, cell=cell),
        "2C3  C3^2(z)": rotation_matrix("z", 240, cell=cell),
        "3C2  C2": rotation_matrix("x", 180, cell=cell),
        "i  I": inversion(cell),
        "2S6  S6": improper_rotation_matrix("z", 60, cell=cell),
        "3sigmad  sigmad1": reflection_matrix("xy", cell=cell),
        "3sigmad  sigmad2": reflection_matrix_axis(sigmad2_axis, cell=cell),
        "3sigmad  sigmad3": reflection_matrix_axis(sigmad3_axis, cell=cell),
    }
    return point_groups


def make_hexagonal(cell=None):
    """
    六方晶系の点群（C6, C3h, C6h, D6, C6v, D3h, D6h）の対称操作をO(3)行列で返す
    cell: 格子情報（必要に応じて渡す）
    return: dict {点群名: {操作名: 3x3行列}}
    """
    point_groups = {}

    # C6 (6回回転) - 完全版
    point_groups["C6"] = {
        "E  E": identity(),
        "C6  C6": rotation_matrix("z", 60, cell),
        "C6_2  C6^2": rotation_matrix("z", 120, cell),
        "C6_3  C6^3": rotation_matrix("z", 180, cell),
        "C6_4  C6^4": rotation_matrix("z", 240, cell),
        "C6_5  C6^5": rotation_matrix("z", 300, cell),
    }

    # C3h (3回回転＋水平鏡映) - 完全版
    point_groups["C3h"] = {
        "E  E": identity(),
        "C3  C3": rotation_matrix("z", 120, cell),
        "C3_2  C3^2": rotation_matrix("z", 240, cell),
        "sigmah  sigmah": reflection_matrix("xy", cell),
        "S3  S3": improper_rotation_matrix("z", 120, cell),
        "S3_2  S3^2": improper_rotation_matrix("z", 240, cell),
    }

    # C6h (6回回転＋水平鏡映＋反転) - 完全版
    point_groups["C6h"] = {
        "E  E": identity(),
        "C6  C6": rotation_matrix("z", 60, cell),
        "C6_2  C6^2": rotation_matrix("z", 120, cell),
        "C6_3  C6^3": rotation_matrix("z", 180, cell),
        "C6_4  C6^4": rotation_matrix("z", 240, cell),
        "C6_5  C6^5": rotation_matrix("z", 300, cell),
        "i  I": inversion(cell),
        "S3  S3": improper_rotation_matrix("z", 120, cell),
        "S3_5  S3^5": improper_rotation_matrix("z", 240, cell),
        "sigmah  sigmah": reflection_matrix("xy", cell),
        "S6  S6": improper_rotation_matrix("z", 60, cell),
        "S6_5  S6^5": improper_rotation_matrix("z", 300, cell),
    }

    # D6 (6回回転＋2回回転×6) - 完全版
    point_groups["D6"] = {
        "E  E": identity(),
        "C6  C6": rotation_matrix("z", 60, cell),
        "C6_2  C6^2": rotation_matrix("z", 120, cell),
        "C6_3  C6^3": rotation_matrix("z", 180, cell),
        "C6_4  C6^4": rotation_matrix("z", 240, cell),
        "C6_5  C6^5": rotation_matrix("z", 300, cell),
        "C2x  C2(x)": rotation_matrix("x", 180, cell),
        "C2y  C2(y)": rotation_matrix("y", 180, cell),
        "C2_1  C2(1)": rotation_matrix_axis([1, 1, 0], 180, cell),  # xy面内の対角線
        "C2_2  C2(2)": rotation_matrix_axis([1, -1, 0], 180, cell),  # xy面内の対角線
        "C2_3  C2(3)": rotation_matrix_axis(
            [np.sqrt(3) / 2, 1 / 2, 0], 180, cell
        ),  # 30度回転した軸
        "C2_4  C2(4)": rotation_matrix_axis(
            [np.sqrt(3) / 2, -1 / 2, 0], 180, cell
        ),  # -30度回転した軸
    }

    # C6v (6回回転＋鏡映×6) - 完全版
    # 鏡映面の法線ベクトルを計算
    sigmav1_axis = [1, 0, 0]  # xz面
    sigmav2_axis = np.dot(rotation_matrix("z", 60, cell=None), [1, 0, 0])  # 60度回転
    sigmav3_axis = np.dot(rotation_matrix("z", 120, cell=None), [1, 0, 0])  # 120度回転
    sigmav4_axis = np.dot(rotation_matrix("z", 180, cell=None), [1, 0, 0])  # 180度回転
    sigmav5_axis = np.dot(rotation_matrix("z", 240, cell=None), [1, 0, 0])  # 240度回転
    sigmav6_axis = np.dot(rotation_matrix("z", 300, cell=None), [1, 0, 0])  # 300度回転

    point_groups["C6v"] = {
        "E  E": identity(),
        "C6  C6": rotation_matrix("z", 60, cell),
        "C6_2  C6^2": rotation_matrix("z", 120, cell),
        "C6_3  C6^3": rotation_matrix("z", 180, cell),
        "C6_4  C6^4": rotation_matrix("z", 240, cell),
        "C6_5  C6^5": rotation_matrix("z", 300, cell),
        "sigmav1  sigmav1": reflection_matrix_axis(sigmav1_axis, cell),
        "sigmav2  sigmav2": reflection_matrix_axis(sigmav2_axis, cell),
        "sigmav3  sigmav3": reflection_matrix_axis(sigmav3_axis, cell),
        "sigmav4  sigmav4": reflection_matrix_axis(sigmav4_axis, cell),
        "sigmav5  sigmav5": reflection_matrix_axis(sigmav5_axis, cell),
        "sigmav6  sigmav6": reflection_matrix_axis(sigmav6_axis, cell),
    }

    # D6h (6回回転＋2回回転×6＋反転＋鏡映×7)
    # z軸: C6, C3, C2, S3, S6, σh
    # 3C2': x, y, -x-y軸
    # 3C2'': x', y', z'軸 (x'=x+y, y'=y-x, z'=z)
    # 3σd: C2'軸を含む面, 3σv: C2''軸を含む面
    sqrt3 = np.sqrt(3)
    # C2'軸
    c2p_axes = [
        [1, 0, 0],
        [-0.5, sqrt3/2, 0],
        [-0.5, -sqrt3/2, 0],
    ]
    # C2''軸
    c2pp_axes = [
        [0, 1, 0],
        [sqrt3/2, 0.5, 0],
        [-sqrt3/2, 0.5, 0],
    ]
    # σd面法線（C2'軸とz軸の中間）
    sigmad_axes = [
        np.array([1, 0, 0]) + np.array([0, 0, 1]),
        np.array([-0.5, sqrt3/2, 0]) + np.array([0, 0, 1]),
        np.array([-0.5, -sqrt3/2, 0]) + np.array([0, 0, 1]),
    ]
    # σv面法線（C2''軸とz軸の中間）
    sigmav_axes = [
        np.array([0, 1, 0]) + np.array([0, 0, 1]),
        np.array([sqrt3/2, 0.5, 0]) + np.array([0, 0, 1]),
        np.array([-sqrt3/2, 0.5, 0]) + np.array([0, 0, 1]),
    ]
    point_groups["D6h"] = {
        "E  E": identity(),
        "2C6  C6(z)": rotation_matrix("z", 60, cell),
        "2C6  C6^5(z)": rotation_matrix("z", 300, cell),
        "2C3  C3(z)": rotation_matrix("z", 120, cell),
        "2C3  C3^2(z)": rotation_matrix("z", 240, cell),
        "C2  C2(z)": rotation_matrix("z", 180, cell),
        # 3C2'
        **{
            f"3C2_1  C2'({i+1})": rotation_matrix_axis(np.array(axis), 180, cell)
            for i, axis in enumerate(c2p_axes)
        },
        # 3C2''
        **{
            f"3C2_2  C2''({i+1})": rotation_matrix_axis(np.array(axis), 180, cell)
            for i, axis in enumerate(c2pp_axes)
        },
        "i  I": inversion(cell),
        "2S3  S3(z)": improper_rotation_matrix("z", 120, cell),
        "2S3  S3^2(z)": improper_rotation_matrix("z", 240, cell),
        "2S6  S6(z)": improper_rotation_matrix("z", 60, cell),
        "2S6  S6^5(z)": improper_rotation_matrix("z", 300, cell),
        "sigmah  sigmah(xy)": reflection_matrix("xy", cell),
        # 3σd
        **{
            f"3sigmad  sigmad({i+1})": reflection_matrix_axis(axis/np.linalg.norm(axis), cell)
            for i, axis in enumerate(sigmad_axes)
        },
        # 3σv
        **{
            f"3sigmav  sigmav({i+1})": reflection_matrix_axis(axis/np.linalg.norm(axis), cell)
            for i, axis in enumerate(sigmav_axes)
        },
    }

    # 他の重要な点群も追加
    # D3h (3回回転＋2回回転×3＋水平鏡映＋鏡映×3)
    point_groups["D3h"] = {
        "E  E": identity(),
        "2C3  C3": rotation_matrix("z", 120, cell),
        "2C3  C3^2": rotation_matrix("z", 240, cell),
        "3C2  C2(x)": rotation_matrix("x", 180, cell),
        "3C2  C2(1)": rotation_matrix_axis([np.sqrt(3) / 2, 1 / 2, 0], 180, cell),
        "3C2  C2(2)": rotation_matrix_axis([np.sqrt(3) / 2, -1 / 2, 0], 180, cell),
        "sigmah  sigmah": reflection_matrix("xy", cell),
        "2S3  S3": improper_rotation_matrix("z", 120, cell),
        "2S3  S3^5": improper_rotation_matrix("z", 240, cell),
        "3sigmav  sigmav1": reflection_matrix("xz", cell),
        "3sigmav  sigmav2": reflection_matrix_axis([np.sqrt(3) / 2, 1 / 2, 0], cell),
        "3sigmav  sigmav3": reflection_matrix_axis([np.sqrt(3) / 2, -1 / 2, 0], cell),
    }

    # D3d (3回回転＋2回回転×3＋反転＋対角鏡映×3)
    point_groups["D3d"] = {
        "E  E": identity(),
        "2C3  C3": rotation_matrix("z", 120, cell),
        "2C3  C3^2": rotation_matrix("z", 240, cell),
        "3C2  C2(x)": rotation_matrix("x", 180, cell),
        "3C2  C2(1)": rotation_matrix_axis([np.sqrt(3) / 2, 1 / 2, 0], 180, cell),
        "3C2  C2(2)": rotation_matrix_axis([np.sqrt(3) / 2, -1 / 2, 0], 180, cell),
        "i  I": inversion(cell),
        "2S6  S6": improper_rotation_matrix("z", 60, cell),
        "2S6  S6^5": improper_rotation_matrix("z", 300, cell),
        "3sigmad  sigmad1": reflection_matrix_axis([1, 1, 0], cell),
        "3sigmad  sigmad2": reflection_matrix_axis([np.sqrt(3) / 2, 1 / 2, 0], cell),
        "3sigmad  sigmad3": reflection_matrix_axis([np.sqrt(3) / 2, -1 / 2, 0], cell),
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
    # 4本のC3（三重軸, 立方体対角線）・3本のC2（xyz軸）
    c3_axes = [
        [1, 1, 1],
        [-1, 1, 1],
        [1, -1, 1],
        [1, 1, -1],
    ]
    c2_axes = [
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1],
    ]
    sqrt3 = np.sqrt(3)
    # 3C2: xyz軸, 8C3: 4本の三重軸ごとに±120度
    point_groups["T"] = {
        "E  E": identity(),
        **{
            f"8C3  C3({i+1})": rotation_matrix_axis(np.array(axis) / sqrt3, 120, cell)
            for i, axis in enumerate(c3_axes)
        },
        **{
            f"8C3  C3^2({i+1})": rotation_matrix_axis(np.array(axis) / sqrt3, 240, cell)
            for i, axis in enumerate(c3_axes)
        },
        "3C2  C2(x)": rotation_matrix("x", 180, cell),
        "3C2  C2(y)": rotation_matrix("y", 180, cell),
        "3C2  C2(z)": rotation_matrix("z", 180, cell),
    }
    # Th (T＋反転)
    point_groups["Th"] = {
        "E  E": identity(),
        **{
            f"8C3  C3({i+1})": rotation_matrix_axis(np.array(axis) / sqrt3, 120, cell)
            for i, axis in enumerate(c3_axes)
        },
        **{
            f"8C3  C3^2({i+1})": rotation_matrix_axis(np.array(axis) / sqrt3, 240, cell)
            for i, axis in enumerate(c3_axes)
        },
        "3C2  C2(x)": rotation_matrix("x", 180, cell),
        "3C2  C2(y)": rotation_matrix("y", 180, cell),
        "3C2  C2(z)": rotation_matrix("z", 180, cell),
        "i  I": inversion(cell),
    }
    # O (四重回転＋三重回転＋2回回転)
    # 4本のC3（三重軸）・3本のC2（xyz軸）・3本のC4（xyz軸, ±90°）
    point_groups["O"] = {
        "E  E": identity(),
        **{
            f"8C3  C3({i+1})": rotation_matrix_axis(np.array(axis) / sqrt3, 120, cell)
            for i, axis in enumerate(c3_axes)
        },
        **{
            f"8C3  C3^2({i+1})": rotation_matrix_axis(np.array(axis) / sqrt3, 240, cell)
            for i, axis in enumerate(c3_axes)
        },
        "6C2  C2(x)": rotation_matrix("x", 180, cell),
        "6C2  C2(y)": rotation_matrix("y", 180, cell),
        "6C2  C2(z)": rotation_matrix("z", 180, cell),
        "6C4  C4(x)": rotation_matrix("x", 90, cell),
        "6C4  C4^3(x)": rotation_matrix("x", 270, cell),
        "6C4  C4(y)": rotation_matrix("y", 90, cell),
        "6C4  C4^3(y)": rotation_matrix("y", 270, cell),
        "6C4  C4(z)": rotation_matrix("z", 90, cell),
        "6C4  C4^3(z)": rotation_matrix("z", 270, cell),
    }
    # Td (T＋鏡映)
    # 8C3: 4本の三重軸(立方体の対角線)ごとに±120度
    # 3C2: xyz軸
    # 6S4: xyz軸ごとに±90度の不正回転
    # 6sigmad: 立方体の対角面
    sqrt2 = np.sqrt(2)
    sqrt3 = np.sqrt(3)
    c3_axes = [
        [1, 1, 1],
        [-1, 1, 1],
        [1, -1, 1],
        [1, 1, -1],
    ]
    c2_axes = [
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1],
    ]
    sigmad_axes = [
        [1, 1, 0],
        [1, -1, 0],
        [1, 0, 1],
        [1, 0, -1],
        [0, 1, 1],
        [0, 1, -1],
    ]
    point_groups["Td"] = {
        "E  E": identity(),
        # 8C3
        **{
            f"8C3  C3({i+1})": rotation_matrix_axis(np.array(axis) / sqrt3, 120, cell)
            for i, axis in enumerate(c3_axes)
        },
        **{
            f"8C3  C3^2({i+1})": rotation_matrix_axis(np.array(axis) / sqrt3, 240, cell)
            for i, axis in enumerate(c3_axes)
        },
        # 3C2
        "3C2  C2(x)": rotation_matrix("x", 180, cell),
        "3C2  C2(y)": rotation_matrix("y", 180, cell),
        "3C2  C2(z)": rotation_matrix("z", 180, cell),
        # 6S4
        "6S4  S4(x)": improper_rotation_matrix("x", 90, cell),
        "6S4  S4^3(x)": improper_rotation_matrix("x", 270, cell),
        "6S4  S4(y)": improper_rotation_matrix("y", 90, cell),
        "6S4  S4^3(y)": improper_rotation_matrix("y", 270, cell),
        "6S4  S4(z)": improper_rotation_matrix("z", 90, cell),
        "6S4  S4^3(z)": improper_rotation_matrix("z", 270, cell),
        # 6sigmad
        **{
            f"6sigmad  sigmad({i+1})": reflection_matrix_axis(
                np.array(axis) / np.linalg.norm(axis), cell
            )
            for i, axis in enumerate(sigmad_axes)
        },
    }
    # Oh (O＋反転＋鏡映)
    # 立方体の全対称操作を網羅的に記述
    # 8C3: 4本の三重軸(立方体の対角線)ごとに±120度
    # 6C2: 3本の二重軸(x,y,z)ごとに180度
    # 6C4: 3本の四重軸(x,y,z)ごとに±90度
    # 3C2: 立方体の面心を結ぶ軸(立方体の辺の中心を通る軸)ごとに180度
    # 6S4: 3本の四重軸(x,y,z)ごとに±90度の不正回転
    # 8S6: 4本の三重軸ごとに±60度の不正回転
    # 3sigmah: xy, yz, zx
    # 6sigmad: 立方体の対角面
    sqrt2 = np.sqrt(2)
    sqrt3 = np.sqrt(3)
    # 三重軸(立方体の対角線)
    c3_axes = [
        [1, 1, 1],
        [-1, 1, 1],
        [1, -1, 1],
        [1, 1, -1],
    ]
    # 面心を結ぶC2軸
    c2_axes = [[0, 1, 1], [0, -1, 1], [1, 0, 1], [-1, 0, 1], [1, 1, 0], [1, -1, 0]]
    # 対角面の法線
    sigmad_axes = [
        [1, 1, 0],
        [1, -1, 0],
        [1, 0, 1],
        [1, 0, -1],
        [0, 1, 1],
        [0, 1, -1],
    ]
    # S6軸(三重軸と同じ)
    s6_axes = c3_axes
    # S4軸(四重軸と同じ)
    # s4_axes = ... 未使用のため削除
    # 四重軸
    # c4_axes, c2_main_axes, sigmah_planes は未使用のため削除
    point_groups["Oh"] = {
        # 単位元
        "E  E": identity(),
        # 8C3
        **{
            f"8C3  C3({i+1})": rotation_matrix_axis(np.array(axis) / sqrt3, 120, cell)
            for i, axis in enumerate(c3_axes)
        },
        **{
            f"8C3  C3^2({i+1})": rotation_matrix_axis(np.array(axis) / sqrt3, 240, cell)
            for i, axis in enumerate(c3_axes)
        },
        # 6C2 (xyz)
        "6C2  C2(x)": rotation_matrix("x", 180, cell),
        "6C2  C2(y)": rotation_matrix("y", 180, cell),
        "6C2  C2(z)": rotation_matrix("z", 180, cell),
        # 6C4 (xyz)
        "6C4  C4(x)": rotation_matrix("x", 90, cell),
        "6C4  C4^3(x)": rotation_matrix("x", 270, cell),
        "6C4  C4(y)": rotation_matrix("y", 90, cell),
        "6C4  C4^3(y)": rotation_matrix("y", 270, cell),
        "6C4  C4(z)": rotation_matrix("z", 90, cell),
        "6C4  C4^3(z)": rotation_matrix("z", 270, cell),
        # 3C2 (面心)
        **{
            f"3C2  C2(f{i+1})": rotation_matrix_axis(np.array(axis) / sqrt2, 180, cell)
            for i, axis in enumerate(c2_axes)
        },
        # 反転
        "i  I": inversion(cell),
        # 6S4 (xyz)
        "6S4  S4(x)": improper_rotation_matrix("x", 90, cell),
        "6S4  S4^3(x)": improper_rotation_matrix("x", 270, cell),
        "6S4  S4(y)": improper_rotation_matrix("y", 90, cell),
        "6S4  S4^3(y)": improper_rotation_matrix("y", 270, cell),
        "6S4  S4(z)": improper_rotation_matrix("z", 90, cell),
        "6S4  S4^3(z)": improper_rotation_matrix("z", 270, cell),
        # 8S6 (三重軸)
        **{
            f"8S6  S6({i+1})": improper_rotation_matrix_axis(
                axis / np.linalg.norm(axis), 60, cell
            )
            for i, axis in enumerate(s6_axes)
        },
        **{
            f"8S6  S6^5({i+1})": improper_rotation_matrix_axis(
                axis / np.linalg.norm(axis), 300, cell
            )
            for i, axis in enumerate(s6_axes)
        },
        # 3sigmah
        "3sigmah  sigmah(xy)": reflection_matrix("xy", cell),
        "3sigmah  sigmah(yz)": reflection_matrix("yz", cell),
        "3sigmah  sigmah(xz)": reflection_matrix("xz", cell),
        # 6sigmad (対角面)
        **{
            f"6sigmad  sigmad({i+1})": reflection_matrix_axis(
                np.array(axis) / np.linalg.norm(axis), cell
            )
            for i, axis in enumerate(sigmad_axes)
        },
    }
    return point_groups
