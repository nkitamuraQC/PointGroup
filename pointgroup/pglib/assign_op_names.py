

import numpy as np


def assign_op_name(rot, tr=None, tol=1e-6, pg_name=None, crystal_system=None):
    """
    spglibの回転行列(rot: 3x3 int)と並進(tr: 3要素)から操作名(C2, C3, S4, sigma, i, E, ...)を結晶系ごとに判定して返す。
    crystal_system: "triclinic", "monoclinic", "orthorhombic", "tetragonal", "trigonal", "hexagonal", "cubic" など
    """
    if crystal_system == "triclinic":
        return _assign_op_name_triclinic(rot, tr, tol, pg_name)
    elif crystal_system == "monoclinic":
        return _assign_op_name_monoclinic(rot, tr, tol, pg_name)
    elif crystal_system == "orthorhombic":
        return _assign_op_name_orthorhombic(rot, tr, tol, pg_name)
    elif crystal_system == "tetragonal":
        return _assign_op_name_tetragonal(rot, tr, tol, pg_name)
    elif crystal_system == "trigonal":
        return _assign_op_name_trigonal(rot, tr, tol, pg_name)
    elif crystal_system == "hexagonal":
        return _assign_op_name_hexagonal(rot, tr, tol, pg_name)
    elif crystal_system == "cubic":
        return _assign_op_name_cubic(rot, tr, tol, pg_name)
    else:
        return _assign_op_name_general(rot, tr, tol, pg_name)

# --- 各結晶系ごとの判定関数 ---
def _assign_op_name_triclinic(rot, tr, tol, pg_name):
    # C1, Ciのみ
    if np.allclose(rot, np.eye(3), atol=tol):
        return 'E'
    if np.allclose(rot, -np.eye(3), atol=tol):
        return 'i'
    return 'unknown'

def _assign_op_name_monoclinic(rot, tr, tol, pg_name):
    # y軸中心のC2, sigma, E, i
    if np.allclose(rot, np.eye(3), atol=tol):
        return 'E'
    if np.allclose(rot, -np.eye(3), atol=tol):
        return 'i'
    eigvals, eigvecs = np.linalg.eig(rot)
    eigvals = np.round(eigvals.real, 6)
    # C2(y)
    if np.count_nonzero(np.abs(eigvals - 1) < tol) == 1:
        idx = np.argmin(np.abs(eigvals - 1))
        axis = eigvecs[:, idx].real
        axis_name = _rotation_axis_name(axis, pg_name)
        n = _find_rotation_order(rot, tol=tol)
        if n == 2 and axis_name == "y":
            return "C2(y)"
    # sigma(xz)
    if np.count_nonzero(np.abs(eigvals - 1) < tol) == 2 and np.count_nonzero(np.abs(eigvals + 1) < tol) == 1:
        normal = eigvecs[:, np.argmin(np.abs(eigvals + 1))].real
        plane = _mirror_plane_name(normal, pg_name)
        if plane == "sigmav(xz)":
            return plane
    return 'unknown'

def _assign_op_name_orthorhombic(rot, tr, tol, pg_name):
    # x, y, z軸C2, 各鏡映, E, i
    if np.allclose(rot, np.eye(3), atol=tol):
        return 'E'
    if np.allclose(rot, -np.eye(3), atol=tol):
        return 'i'
    eigvals, eigvecs = np.linalg.eig(rot)
    eigvals = np.round(eigvals.real, 6)
    # C2(x), C2(y), C2(z)
    if np.count_nonzero(np.abs(eigvals - 1) < tol) == 1:
        idx = np.argmin(np.abs(eigvals - 1))
        axis = eigvecs[:, idx].real
        axis_name = _rotation_axis_name(axis, pg_name)
        n = _find_rotation_order(rot, tol=tol)
        if n == 2 and axis_name in ("x", "y", "z"):
            return f"C2({axis_name})"
    # 各鏡映
    if np.count_nonzero(np.abs(eigvals - 1) < tol) == 2 and np.count_nonzero(np.abs(eigvals + 1) < tol) == 1:
        normal = eigvecs[:, np.argmin(np.abs(eigvals + 1))].real
        plane = _mirror_plane_name(normal, pg_name)
        if plane in ("sigmav(xz)", "sigmav(yz)", "sigmah"):
            return plane
    return 'unknown'

def _assign_op_name_tetragonal(rot, tr, tol, pg_name):
    # z軸C4, C2, S4, sigma, E, i
    if np.allclose(rot, np.eye(3), atol=tol):
        return 'E'
    if np.allclose(rot, -np.eye(3), atol=tol):
        return 'i'
    eigvals, eigvecs = np.linalg.eig(rot)
    eigvals = np.round(eigvals.real, 6)
    # C4(z), C2(z)
    if np.count_nonzero(np.abs(eigvals - 1) < tol) == 1:
        idx = np.argmin(np.abs(eigvals - 1))
        axis = eigvecs[:, idx].real
        axis_name = _rotation_axis_name(axis, pg_name)
        n = _find_rotation_order(rot, tol=tol)
        if n == 4 and axis_name == "z":
            return "C4"
        if n == 2 and axis_name == "z":
            return "C2(z)"
    # S4(z)
    if np.linalg.det(rot) < 0:
        n = _find_rotation_order(rot, tol=tol)
        if n == 4:
            return "S4"
    # 鏡映
    if np.count_nonzero(np.abs(eigvals - 1) < tol) == 2 and np.count_nonzero(np.abs(eigvals + 1) < tol) == 1:
        normal = eigvecs[:, np.argmin(np.abs(eigvals + 1))].real
        plane = _mirror_plane_name(normal, pg_name)
        if plane in ("sigmah", "sigmav(xz)", "sigmav(yz)", "sigmad"):
            return plane
    return 'unknown'

def _assign_op_name_trigonal(rot, tr, tol, pg_name):
    # z軸C3, S6, sigma, E, i
    if np.allclose(rot, np.eye(3), atol=tol):
        return 'E'
    if np.allclose(rot, -np.eye(3), atol=tol):
        return 'i'
    eigvals, eigvecs = np.linalg.eig(rot)
    eigvals = np.round(eigvals.real, 6)
    # C3(z)
    if np.count_nonzero(np.abs(eigvals - 1) < tol) == 1:
        idx = np.argmin(np.abs(eigvals - 1))
        axis = eigvecs[:, idx].real
        axis_name = _rotation_axis_name(axis, pg_name)
        n = _find_rotation_order(rot, tol=tol)
        if n == 3 and axis_name == "z":
            return "C3"
    # S6(z)
    if np.linalg.det(rot) < 0:
        n = _find_rotation_order(rot, tol=tol)
        if n == 6:
            return "S6"
    # 鏡映
    if np.count_nonzero(np.abs(eigvals - 1) < tol) == 2 and np.count_nonzero(np.abs(eigvals + 1) < tol) == 1:
        normal = eigvecs[:, np.argmin(np.abs(eigvals + 1))].real
        plane = _mirror_plane_name(normal, pg_name)
        if plane in ("sigmah", "sigmav(xz)", "sigmav(yz)", "sigmad"):
            return plane
    return 'unknown'

def _assign_op_name_hexagonal(rot, tr, tol, pg_name):
    # z軸C6, C3, C2, S6, sigma, E, i
    if np.allclose(rot, np.eye(3), atol=tol):
        return 'E'
    if np.allclose(rot, -np.eye(3), atol=tol):
        return 'i'
    eigvals, eigvecs = np.linalg.eig(rot)
    eigvals = np.round(eigvals.real, 6)
    # C6(z), C3(z), C2(z)
    if np.count_nonzero(np.abs(eigvals - 1) < tol) == 1:
        idx = np.argmin(np.abs(eigvals - 1))
        axis = eigvecs[:, idx].real
        axis_name = _rotation_axis_name(axis, pg_name)
        n = _find_rotation_order(rot, tol=tol)
        if n == 6 and axis_name == "z":
            return "C6"
        if n == 3 and axis_name == "z":
            return "C3"
        if n == 2 and axis_name == "z":
            return "C2(z)"
    # S6(z)
    if np.linalg.det(rot) < 0:
        n = _find_rotation_order(rot, tol=tol)
        if n == 6:
            return "S6"
    # 鏡映
    if np.count_nonzero(np.abs(eigvals - 1) < tol) == 2 and np.count_nonzero(np.abs(eigvals + 1) < tol) == 1:
        normal = eigvecs[:, np.argmin(np.abs(eigvals + 1))].real
        plane = _mirror_plane_name(normal, pg_name)
        if plane in ("sigmah", "sigmav(xz)", "sigmav(yz)", "sigmad"):
            return plane
    return 'unknown'

def _assign_op_name_cubic(rot, tr, tol, pg_name):
    # C3, C4, C2, S4, S6, sigma, E, i など多様
    if np.allclose(rot, np.eye(3), atol=tol):
        return 'E'
    if np.allclose(rot, -np.eye(3), atol=tol):
        return 'i'
    eigvals, eigvecs = np.linalg.eig(rot)
    eigvals = np.round(eigvals.real, 6)
    # C4(z), C3([111]), C2(x), C2(y), C2(z)
    if np.count_nonzero(np.abs(eigvals - 1) < tol) == 1:
        idx = np.argmin(np.abs(eigvals - 1))
        axis = eigvecs[:, idx].real
        axis_name = _rotation_axis_name(axis, pg_name)
        n = _find_rotation_order(rot, tol=tol)
        if n == 4 and axis_name == "z":
            return "C4"
        if n == 3 and axis_name == "d":
            return "C3"
        if n == 2 and axis_name in ("x", "y", "z"):
            return f"C2({axis_name})"
    # S4(z), S6(d) など
    if np.linalg.det(rot) < 0:
        # S4(z)
        n = _find_rotation_order(rot, tol=tol)
        if n == 4:
            return "S4"
        # S6([111])
        if n == 6:
            # S6軸が[111]方向か判定
            if np.count_nonzero(np.abs(eigvals - 1) < tol) == 1:
                idx = np.argmin(np.abs(eigvals - 1))
                axis = eigvecs[:, idx].real
                axis_name = _rotation_axis_name(axis, pg_name)
                if axis_name == "d":
                    return "S6"
    # 鏡映
    if np.count_nonzero(np.abs(eigvals - 1) < tol) == 2 and np.count_nonzero(np.abs(eigvals + 1) < tol) == 1:
        normal = eigvecs[:, np.argmin(np.abs(eigvals + 1))].real
        plane = _mirror_plane_name(normal, pg_name)
        if plane in ("sigmah", "sigmav(xz)", "sigmav(yz)", "sigmad"):
            return plane
    return 'unknown'

def _assign_op_name_general(rot, tr, tol, pg_name):
    # 一般的な判定（従来のロジック）
    if np.allclose(rot, np.eye(3), atol=tol):
        if tr is not None and not np.allclose(tr, 0, atol=tol):
            return f't({tuple(np.round(tr, 6))})'
        return 'E'
    if np.allclose(rot, -np.eye(3), atol=tol):
        return 'i'
    eigvals, eigvecs = np.linalg.eig(rot)
    eigvals = np.round(eigvals.real, 6)
    if np.count_nonzero(np.abs(eigvals - 1) < tol) == 2 and np.count_nonzero(np.abs(eigvals + 1) < tol) == 1:
        normal = eigvecs[:, np.argmin(np.abs(eigvals + 1))].real
        plane = _mirror_plane_name(normal, pg_name)
        return plane
    if np.count_nonzero(np.abs(eigvals - 1) < tol) == 1:
        idx = np.argmin(np.abs(eigvals - 1))
        axis = eigvecs[:, idx].real
        axis_name = _rotation_axis_name(axis, pg_name)
        n = _find_rotation_order(rot, tol=tol)
        if n is not None and n > 1:
            if n == 2 and axis_name in ("x", "y", "z"):
                return f"C2({axis_name})"
            if n == 3 and axis_name == "z":
                return "C3"
            if n == 4 and axis_name == "z":
                return "C4"
            if n == 6 and axis_name == "z":
                return "C6"
            return f"C{n}({axis_name})"
    if np.linalg.det(rot) < 0:
        if np.count_nonzero(np.abs(eigvals - 1) < tol) == 1:
            idx = np.argmin(np.abs(eigvals - 1))
            axis = eigvecs[:, idx].real
            axis_name = _rotation_axis_name(axis, pg_name)
        else:
            axis_name = '?'
        n = _find_rotation_order(rot, tol=tol)
        if n is not None and n > 1:
            if n == 4 and axis_name == "z":
                return "S4"
            if n == 6 and axis_name == "z":
                return "S6"
            return f"S{n}({axis_name})"
        return f"S操作({axis_name})"
    return 'unknown'

def _mirror_plane_name(normal, pg_name=None):
    # 鏡映面の法線ベクトルからpg_ch_extraの命名規則に合わせた名前を返す
    normal = np.real(normal)
    normal = normal / (np.max(np.abs(normal)) if np.max(np.abs(normal)) > 0 else 1)
    rounded = np.rint(normal * 3).astype(int)
    g = np.gcd.reduce(np.abs(rounded))
    if g > 0:
        rounded = rounded // g
    # xy, xz, yz, d, h, v
    if np.allclose(rounded, [0, 0, 1]):
        return "sigmah"
    if np.allclose(rounded, [1, 0, 0]):
        return "sigmav(yz)"
    if np.allclose(rounded, [0, 1, 0]):
        return "sigmav(xz)"
    # 対角面
    if np.allclose(np.abs(rounded), [1, 1, 0]):
        return "sigmad"
    # その他
    return f"sigma({rounded.tolist()})"

def _rotation_axis_name(axis, pg_name=None):
    # 回転軸ベクトルからpg_ch_extraの命名規則に合わせた軸名を返す
    axis = np.real(axis)
    axis = axis / (np.max(np.abs(axis)) if np.max(np.abs(axis)) > 0 else 1)
    rounded = np.rint(axis * 3).astype(int)
    g = np.gcd.reduce(np.abs(rounded))
    if g > 0:
        rounded = rounded // g
    if np.allclose(rounded, [0, 0, 1]):
        return "z"
    if np.allclose(rounded, [1, 0, 0]):
        return "x"
    if np.allclose(rounded, [0, 1, 0]):
        return "y"
    if np.allclose(np.abs(rounded), [1, 1, 0]):
        return "xy"
    if np.allclose(np.abs(rounded), [1, 0, 1]):
        return "xz"
    if np.allclose(np.abs(rounded), [0, 1, 1]):
        return "yz"
    if np.allclose(np.abs(rounded), [1, 1, 1]):
        return "d"
    return str(rounded.tolist())

def _axis_to_str(axis):
    # 近い整数比でxyz軸名や[111]型を返す
    axis = np.real(axis)
    axis = axis / (np.max(np.abs(axis)) if np.max(np.abs(axis)) > 0 else 1)
    rounded = np.rint(axis * 3).astype(int)  # 例: [1,1,1]や[1,0,0]など
    g = np.gcd.reduce(np.abs(rounded))
    if g > 0:
        rounded = rounded // g
    labels = ['x', 'y', 'z']
    if np.count_nonzero(rounded) == 1:
        idx = np.argmax(np.abs(rounded))
        return labels[idx]
    return '[' + ''.join(str(i) for i in rounded if i != 0) + ']'

def _find_rotation_order(rot, tol=1e-6):
    # rot^n = I となる最小のn (1 < n <= 6) を返す。なければNone。
    for n in range(2, 7):
        if np.allclose(np.linalg.matrix_power(rot, n), np.eye(3), atol=tol):
            return n
    return None