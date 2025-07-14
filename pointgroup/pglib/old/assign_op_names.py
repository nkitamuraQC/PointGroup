import numpy as np

def assign_op_name(rot, pg_name, tr=None, tol=1e-6):
    """
    点群名pg_nameごとに、その点群に属する操作の表現行列rotに名前をつけて返す。
    例: pg_name="D2h"ならC2(x), C2(y), C2(z), sigma(xy), sigma(xz), sigma(yz), i, Eなど。
    """
    pg_name_l = pg_name.lower() if pg_name else None
    # 点群ごとに分岐
    if pg_name_l in ("c1", "1"):
        return _assign_op_name_c1(rot, tr, tol)
    if pg_name_l in ("ci", "-1"):
        return _assign_op_name_ci(rot, tr, tol)
    if pg_name_l in ("c2", "2"):
        return _assign_op_name_c2(rot, tr, tol)
    if pg_name_l in ("cs", "m"):
        return _assign_op_name_cs(rot, tr, tol)
    if pg_name_l in ("c2h", "2/m"):
        return _assign_op_name_c2h(rot, tr, tol)
    if pg_name_l in ("d2", ):
        return _assign_op_name_d2(rot, tr, tol)
    if pg_name_l in ("c2v", ):
        return _assign_op_name_c2v(rot, tr, tol)
    if pg_name_l in ("d2h", ):
        return _assign_op_name_d2h(rot, tr, tol)
    if pg_name_l in ("c4", ):
        return _assign_op_name_c4(rot, tr, tol)
    if pg_name_l in ("s4", ):
        return _assign_op_name_s4(rot, tr, tol)
    if pg_name_l in ("c4h", ):
        return _assign_op_name_c4h(rot, tr, tol)
    if pg_name_l in ("d4", ):
        return _assign_op_name_d4(rot, tr, tol)
    if pg_name_l in ("c4v", ):
        return _assign_op_name_c4v(rot, tr, tol)
    if pg_name_l in ("d2d", ):
        return _assign_op_name_d2d(rot, tr, tol)
    if pg_name_l in ("d4h", ):
        return _assign_op_name_d4h(rot, tr, tol)
    if pg_name_l in ("c3", ):
        return _assign_op_name_c3(rot, tr, tol)
    if pg_name_l in ("c3i", ):
        return _assign_op_name_c3i(rot, tr, tol)
    if pg_name_l in ("d3", ):
        return _assign_op_name_d3(rot, tr, tol)
    if pg_name_l in ("c3v", ):
        return _assign_op_name_c3v(rot, tr, tol)
    if pg_name_l in ("d3d", ):
        return _assign_op_name_d3d(rot, tr, tol)
    if pg_name_l in ("c6", ):
        return _assign_op_name_c6(rot, tr, tol)
    if pg_name_l in ("c3h", ):
        return _assign_op_name_c3h(rot, tr, tol)
    if pg_name_l in ("c6h", ):
        return _assign_op_name_c6h(rot, tr, tol)
    if pg_name_l in ("d6", ):
        return _assign_op_name_d6(rot, tr, tol)
    if pg_name_l in ("c6v", ):
        return _assign_op_name_c6v(rot, tr, tol)
    if pg_name_l in ("d3h", ):
        return _assign_op_name_d3h(rot, tr, tol)
    if pg_name_l in ("d6h", ):
        return _assign_op_name_d6h(rot, tr, tol)
    if pg_name_l in ("t", ):
        return _assign_op_name_t(rot, tr, tol)
    if pg_name_l in ("th", ):
        return _assign_op_name_th(rot, tr, tol)
    if pg_name_l in ("o", ):
        return _assign_op_name_o(rot, tr, tol)
    if pg_name_l in ("td", ):
        return _assign_op_name_td(rot, tr, tol)
    if pg_name_l in ("oh", ):
        return _assign_op_name_oh(rot, tr, tol)
    # 未知の点群は一般判定
    return _assign_op_name_general(rot, tr, tol, pg_name)

def _assign_op_name_c1(rot, tr, tol):
    if np.allclose(rot, np.eye(3), atol=tol):
        return "E"
    return "unknown"

def _assign_op_name_ci(rot, tr, tol):
    if np.allclose(rot, np.eye(3), atol=tol):
        return "E"
    if np.allclose(rot, -np.eye(3), atol=tol):
        return "i"
    return "unknown"

def _assign_op_name_c2(rot, tr, tol):
    if np.allclose(rot, np.eye(3), atol=tol):
        return "E"
    # C2: det=1, trace=-1, 固有値1, -1, -1
    if np.isclose(np.linalg.det(rot), 1, atol=tol) and np.isclose(np.trace(rot), -1, atol=tol):
        return "C2"
    return "unknown"

def _assign_op_name_cs(rot, tr, tol):
    if np.allclose(rot, np.eye(3), atol=tol):
        return "E"
    # sigmah: det=-1, trace=1, 固有値1,1,-1
    if np.isclose(np.linalg.det(rot), -1, atol=tol) and np.isclose(np.trace(rot), 1, atol=tol):
        return "sigmah"
    return "unknown"

def _assign_op_name_c2h(rot, tr, tol):
    if np.allclose(rot, np.eye(3), atol=tol):
        return "E"
    if np.allclose(rot, -np.eye(3), atol=tol):
        return "i"
    # C2: det=1, trace=-1
    if np.isclose(np.linalg.det(rot), 1, atol=tol) and np.isclose(np.trace(rot), -1, atol=tol):
        return "C2"
    # sigmah: det=-1, trace=1
    if np.isclose(np.linalg.det(rot), -1, atol=tol) and np.isclose(np.trace(rot), 1, atol=tol):
        return "sigmah"
    return "unknown"

def _assign_op_name_d2(rot, tr, tol):
    if np.allclose(rot, np.eye(3), atol=tol):
        return "E"
    # C2(z), C2(y), C2(x)
    axes = {"C2(z)": [0, 0, 1], "C2(y)": [0, 1, 0], "C2(x)": [1, 0, 0]}
    for name, axis in axes.items():
        mat = _rotation_matrix(axis, np.pi)
        if np.allclose(rot, mat, atol=tol):
            return name
    return "unknown"

def _assign_op_name_c2v(rot, tr, tol):
    if np.allclose(rot, np.eye(3), atol=tol):
        return "E"
    # C2(z)
    mat_c2z = _rotation_matrix([0, 0, 1], np.pi)
    if np.allclose(rot, mat_c2z, atol=tol):
        return "C2(z)"
    # sigmav(xz), sigmav(yz)
    if np.allclose(rot, np.diag([1, -1, 1]), atol=tol):
        return "sigmav(xz)"
    if np.allclose(rot, np.diag([-1, 1, 1]), atol=tol):
        return "sigmav(yz)"
    return "unknown"

def _assign_op_name_d2h(rot, tr, tol):
    if np.allclose(rot, np.eye(3), atol=tol):
        return "E"
    if np.allclose(rot, -np.eye(3), atol=tol):
        return "i"
    # C2(z), C2(y), C2(x)
    axes = {"C2(z)": [0, 0, 1], "C2(y)": [0, 1, 0], "C2(x)": [1, 0, 0]}
    for name, axis in axes.items():
        mat = _rotation_matrix(axis, np.pi)
        if np.allclose(rot, mat, atol=tol):
            return name
    # sigmaxy, sigmaxz, sigmayz
    if np.allclose(rot, np.diag([1, 1, -1]), atol=tol):
        return "sigmaxy"
    if np.allclose(rot, np.diag([1, -1, 1]), atol=tol):
        return "sigmaxz"
    if np.allclose(rot, np.diag([-1, 1, 1]), atol=tol):
        return "sigmayz"
    return "unknown"

# ユーティリティ: 任意軸回転行列


def _assign_op_name_c4(rot, tr, tol):
    if np.allclose(rot, np.eye(3), atol=tol):
        return "E"
    mat_c4 = _rotation_matrix([0,0,1], np.pi/2)
    mat_c2 = _rotation_matrix([0,0,1], np.pi)
    mat_c4_3 = _rotation_matrix([0,0,1], 3*np.pi/2)
    if np.allclose(rot, mat_c4, atol=tol):
        return "C4"
    if np.allclose(rot, mat_c2, atol=tol):
        return "C2"
    if np.allclose(rot, mat_c4_3, atol=tol):
        return "C4_3"
    return "unknown"

def _assign_op_name_s4(rot, tr, tol):
    if np.allclose(rot, np.eye(3), atol=tol):
        return "E"
    mat_s4 = np.array([[0, 1, 0], [-1, 0, 0], [0, 0, -1]])
    if np.allclose(rot, mat_s4, atol=tol):
        return "S4"
    return "unknown"

def _assign_op_name_c4h(rot, tr, tol):
    if np.allclose(rot, np.eye(3), atol=tol):
        return "E"
    if np.allclose(rot, -np.eye(3), atol=tol):
        return "i"
    mat_c4 = _rotation_matrix([0,0,1], np.pi/2)
    mat_c2 = _rotation_matrix([0,0,1], np.pi)
    if np.allclose(rot, mat_c4, atol=tol):
        return "C4"
    if np.allclose(rot, mat_c2, atol=tol):
        return "C2"
    return "unknown"

def _assign_op_name_d4(rot, tr, tol):
    if np.allclose(rot, np.eye(3), atol=tol):
        return "E"
    mat_c4 = _rotation_matrix([0,0,1], np.pi/2)
    mat_c4_3 = _rotation_matrix([0,0,1], 3*np.pi/2)
    mat_c2 = _rotation_matrix([0,0,1], np.pi)
    mat_c2x = _rotation_matrix([1,0,0], np.pi)
    mat_c2y = _rotation_matrix([0,1,0], np.pi)
    if np.allclose(rot, mat_c4, atol=tol) or np.allclose(rot, mat_c4_3, atol=tol):
        return "2C4"
    if np.allclose(rot, mat_c2, atol=tol):
        return "C2"
    if np.allclose(rot, mat_c2x, atol=tol):
        return "2C2_1"
    if np.allclose(rot, mat_c2y, atol=tol):
        return "2C2_2"
    return "unknown"

def _assign_op_name_c4v(rot, tr, tol):
    if np.allclose(rot, np.eye(3), atol=tol):
        return "E"
    mat_c4 = _rotation_matrix([0,0,1], np.pi/2)
    mat_c4_3 = _rotation_matrix([0,0,1], 3*np.pi/2)
    mat_c2 = _rotation_matrix([0,0,1], np.pi)
    if np.allclose(rot, mat_c4, atol=tol) or np.allclose(rot, mat_c4_3, atol=tol):
        return "2C4"
    if np.allclose(rot, mat_c2, atol=tol):
        return "C2"
    if np.allclose(rot, np.diag([1,-1,1]), atol=tol) or np.allclose(rot, np.diag([-1,1,1]), atol=tol):
        return "2sigmav"
    mat_sigmad1 = np.array([[0,1,0],[1,0,0],[0,0,1]])
    mat_sigmad2 = np.array([[0,-1,0],[-1,0,0],[0,0,1]])
    if np.allclose(rot, mat_sigmad1, atol=tol) or np.allclose(rot, mat_sigmad2, atol=tol):
        return "2sigmad"
    return "unknown"

def _assign_op_name_d2d(rot, tr, tol):
    if np.allclose(rot, np.eye(3), atol=tol):
        return "E"
    mat_s4 = np.array([[0, 1, 0], [-1, 0, 0], [0, 0, -1]])
    mat_c2 = _rotation_matrix([0,0,1], np.pi)
    mat_c2d1 = _rotation_matrix([1,1,0], np.pi)
    mat_c2d2 = _rotation_matrix([-1,1,0], np.pi)
    if np.allclose(rot, mat_s4, atol=tol):
        return "2S4"
    if np.allclose(rot, mat_c2, atol=tol):
        return "C2"
    if np.allclose(rot, mat_c2d1, atol=tol) or np.allclose(rot, mat_c2d2, atol=tol):
        return "2C2_1"
    mat_sigmad1 = np.array([[0,1,0],[1,0,0],[0,0,1]])
    mat_sigmad2 = np.array([[0,-1,0],[-1,0,0],[0,0,1]])
    if np.allclose(rot, mat_sigmad1, atol=tol) or np.allclose(rot, mat_sigmad2, atol=tol):
        return "2sigmad"
    return "unknown"

def _assign_op_name_d4h(rot, tr, tol):
    if np.allclose(rot, np.eye(3), atol=tol):
        return "E"
    if np.allclose(rot, -np.eye(3), atol=tol):
        return "i"
    mat_c4 = _rotation_matrix([0,0,1], np.pi/2)
    mat_c4_3 = _rotation_matrix([0,0,1], 3*np.pi/2)
    mat_c2 = _rotation_matrix([0,0,1], np.pi)
    mat_c2x = _rotation_matrix([1,0,0], np.pi)
    mat_c2y = _rotation_matrix([0,1,0], np.pi)
    mat_s4 = np.array([[0, 1, 0], [-1, 0, 0], [0, 0, -1]])
    if np.allclose(rot, mat_c4, atol=tol) or np.allclose(rot, mat_c4_3, atol=tol):
        return "2C4"
    if np.allclose(rot, mat_c2, atol=tol):
        return "C2"
    if np.allclose(rot, mat_c2x, atol=tol):
        return "2C2_1"
    if np.allclose(rot, mat_c2y, atol=tol):
        return "2C2_2"
    if np.allclose(rot, mat_s4, atol=tol):
        return "2S4"
    # 2sigmad: 対角面鏡映 (x↔y, x↔-y, zそのまま or z反転)
    for z_sign in [1, -1]:
        mat1 = np.array([[0, 1, 0], [1, 0, 0], [0, 0, z_sign]])
        mat2 = np.array([[0, -1, 0], [-1, 0, 0], [0, 0, z_sign]])
        mat3 = np.array([[0, -1, 0], [1, 0, 0], [0, 0, z_sign]])
        # 丸めて比較（整数化）
        if (np.allclose(np.round(rot), mat1, atol=tol) or np.allclose(np.round(rot), mat2, atol=tol) or np.allclose(np.round(rot), mat3, atol=tol)):
            return "2sigmad"
    # 2sigmav: 垂直面鏡映 (yz, xz, zそのまま or z反転)
    for z_sign in [1, -1]:
        mat1 = np.diag([1, -1, z_sign])
        mat2 = np.diag([-1, 1, z_sign])
        if (np.allclose(np.round(rot), mat1, atol=tol) or np.allclose(np.round(rot), mat2, atol=tol)):
            return "2sigmav"
    if np.allclose(rot, np.diag([1, 1, -1]), atol=tol):
        return "sigmah"
    return "unknown"

def _assign_op_name_c3(rot, tr, tol):
    if np.allclose(rot, np.eye(3), atol=tol):
        return "E"
    mat_c3 = _rotation_matrix([0,0,1], 2*np.pi/3)
    mat_c3_2 = _rotation_matrix([0,0,1], 4*np.pi/3)
    if np.allclose(rot, mat_c3, atol=tol):
        return "C3"
    if np.allclose(rot, mat_c3_2, atol=tol):
        return "C3_2"
    return "unknown"

def _assign_op_name_c3i(rot, tr, tol):
    if np.allclose(rot, np.eye(3), atol=tol):
        return "E"
    if np.allclose(rot, -np.eye(3), atol=tol):
        return "i"
    mat_c3 = _rotation_matrix([0,0,1], 2*np.pi/3)
    mat_c3_2 = _rotation_matrix([0,0,1], 4*np.pi/3)
    if np.allclose(rot, mat_c3, atol=tol):
        return "C3"
    if np.allclose(rot, mat_c3_2, atol=tol):
        return "C3_2"
    return "unknown"

def _assign_op_name_d3(rot, tr, tol):
    if np.allclose(rot, np.eye(3), atol=tol):
        return "E"
    mat_c3 = _rotation_matrix([0,0,1], 2*np.pi/3)
    mat_c3_2 = _rotation_matrix([0,0,1], 4*np.pi/3)
    mat_c2x = _rotation_matrix([1,0,0], np.pi)
    mat_c2y = _rotation_matrix([0,1,0], np.pi)
    if np.allclose(rot, mat_c3, atol=tol) or np.allclose(rot, mat_c3_2, atol=tol):
        return "2C3"
    if np.allclose(rot, mat_c2x, atol=tol) or np.allclose(rot, mat_c2y, atol=tol):
        return "3C2"
    return "unknown"


def _assign_op_name_c3v(rot, tr, tol):
    if np.allclose(rot, np.eye(3), atol=tol):
        return "E"
    mat_c3 = _rotation_matrix([0,0,1], 2*np.pi/3)
    mat_c3_2 = _rotation_matrix([0,0,1], 4*np.pi/3)
    if np.allclose(rot, mat_c3, atol=tol) or np.allclose(rot, mat_c3_2, atol=tol):
        return "2C3"
    # 3sigmav: 3面鏡映
    sigmas = [np.diag([1,-1,1]), np.diag([-1,1,1]), np.array([[-0.5,0.866,0],[0.866,0.5,0],[0,0,1]]), np.array([[-0.5,-0.866,0],[-0.866,0.5,0],[0,0,1]])]
    for mat in sigmas:
        if np.allclose(rot, mat, atol=tol):
            return "3sigmav"
    return "unknown"

def _assign_op_name_d3d(rot, tr, tol):
    if np.allclose(rot, np.eye(3), atol=tol):
        return "E"
    if np.allclose(rot, -np.eye(3), atol=tol):
        return "i"
    mat_c3 = _rotation_matrix([0,0,1], 2*np.pi/3)
    mat_c3_2 = _rotation_matrix([0,0,1], 4*np.pi/3)
    mat_c2x = _rotation_matrix([1,0,0], np.pi)
    mat_c2y = _rotation_matrix([0,1,0], np.pi)
    if np.allclose(rot, mat_c3, atol=tol) or np.allclose(rot, mat_c3_2, atol=tol):
        return "2C3"
    if np.allclose(rot, mat_c2x, atol=tol) or np.allclose(rot, mat_c2y, atol=tol):
        return "3C2"
    # 2S6
    mat_s6 = np.array([[0.5,0.866,0],[0.866,-0.5,0],[0,0,-1]])
    if np.allclose(rot, mat_s6, atol=tol):
        return "2S6"
    # 3sigmad
    mat_sigmad1 = np.array([[0,1,0],[1,0,0],[0,0,1]])
    mat_sigmad2 = np.array([[0,-1,0],[-1,0,0],[0,0,1]])
    if np.allclose(rot, mat_sigmad1, atol=tol) or np.allclose(rot, mat_sigmad2, atol=tol):
        return "3sigmad"
    return "unknown"

def _assign_op_name_c6(rot, tr, tol):
    if np.allclose(rot, np.eye(3), atol=tol):
        return "E"
    mat_c6 = _rotation_matrix([0,0,1], np.pi/3)
    mat_c3 = _rotation_matrix([0,0,1], 2*np.pi/3)
    mat_c2 = _rotation_matrix([0,0,1], np.pi)
    if np.allclose(rot, mat_c6, atol=tol):
        return "C6"
    if np.allclose(rot, mat_c3, atol=tol):
        return "C3"
    if np.allclose(rot, mat_c2, atol=tol):
        return "C2"
    return "unknown"

def _assign_op_name_c3h(rot, tr, tol):
    if np.allclose(rot, np.eye(3), atol=tol):
        return "E"
    mat_c3 = _rotation_matrix([0,0,1], 2*np.pi/3)
    mat_c3_2 = _rotation_matrix([0,0,1], 4*np.pi/3)
    if np.allclose(rot, mat_c3, atol=tol) or np.allclose(rot, mat_c3_2, atol=tol):
        return "C3"
    # sigmah
    if np.allclose(rot, np.diag([1,1,-1]), atol=tol):
        return "sigmah"
    return "unknown"

def _assign_op_name_c6h(rot, tr, tol):
    if np.allclose(rot, np.eye(3), atol=tol):
        return "E"
    if np.allclose(rot, -np.eye(3), atol=tol):
        return "i"
    mat_c6 = _rotation_matrix([0,0,1], np.pi/3)
    mat_c3 = _rotation_matrix([0,0,1], 2*np.pi/3)
    mat_c2 = _rotation_matrix([0,0,1], np.pi)
    if np.allclose(rot, mat_c6, atol=tol):
        return "C6"
    if np.allclose(rot, mat_c3, atol=tol):
        return "C3"
    if np.allclose(rot, mat_c2, atol=tol):
        return "C2"
    # sigmah
    if np.allclose(rot, np.diag([1,1,-1]), atol=tol):
        return "sigmah"
    return "unknown"

def _assign_op_name_d6(rot, tr, tol):
    if np.allclose(rot, np.eye(3), atol=tol):
        return "E"
    mat_c6 = _rotation_matrix([0,0,1], np.pi/3)
    mat_c3 = _rotation_matrix([0,0,1], 2*np.pi/3)
    mat_c2 = _rotation_matrix([0,0,1], np.pi)
    if np.allclose(rot, mat_c6, atol=tol):
        return "C6"
    if np.allclose(rot, mat_c3, atol=tol):
        return "C3"
    if np.allclose(rot, mat_c2, atol=tol):
        return "C2"
    return "unknown"

def _assign_op_name_c6v(rot, tr, tol):
    if np.allclose(rot, np.eye(3), atol=tol):
        return "E"
    mat_c6 = _rotation_matrix([0,0,1], np.pi/3)
    mat_c3 = _rotation_matrix([0,0,1], 2*np.pi/3)
    mat_c2 = _rotation_matrix([0,0,1], np.pi)
    if np.allclose(rot, mat_c6, atol=tol):
        return "C6"
    if np.allclose(rot, mat_c3, atol=tol):
        return "C3"
    if np.allclose(rot, mat_c2, atol=tol):
        return "C2"
    return "unknown"


def _assign_op_name_d3h(rot, tr, tol):
    if np.allclose(rot, np.eye(3), atol=tol):
        return "E"
    mat_c3 = _rotation_matrix([0,0,1], 2*np.pi/3)
    mat_c3_2 = _rotation_matrix([0,0,1], 4*np.pi/3)
    mat_c2x = _rotation_matrix([1,0,0], np.pi)
    mat_c2y = _rotation_matrix([0,1,0], np.pi)
    if np.allclose(rot, mat_c3, atol=tol) or np.allclose(rot, mat_c3_2, atol=tol):
        return "2C3"
    if np.allclose(rot, mat_c2x, atol=tol) or np.allclose(rot, mat_c2y, atol=tol):
        return "3C2"
    if np.allclose(rot, np.diag([1,1,-1]), atol=tol):
        return "sigmah"
    return "unknown"

def _assign_op_name_d6h(rot, tr, tol):
    # E
    if np.allclose(rot, np.eye(3), atol=tol):
        return "E"
    # i
    if np.allclose(rot, -np.eye(3), atol=tol):
        return "i"
    # 2C6 (z軸)
    for angle in [np.pi/3, 5*np.pi/3]:
        if np.allclose(rot, _rotation_matrix([0, 0, 1], angle), atol=tol):
            return "2C6"
    # 2C3 (z軸)
    for angle in [2*np.pi/3, 4*np.pi/3]:
        if np.allclose(rot, _rotation_matrix([0, 0, 1], angle), atol=tol):
            return "2C3"
    # C2 (z軸)
    if np.allclose(rot, _rotation_matrix([0, 0, 1], np.pi), atol=tol):
        return "C2"
    # 3C2_1 (x, y, -x-y軸)
    c2_1_axes = [[1, 0, 0], [-1/2, np.sqrt(3)/2, 0], [-1/2, -np.sqrt(3)/2, 0]]
    for axis in c2_1_axes:
        if np.allclose(rot, _rotation_matrix(axis, np.pi), atol=tol):
            return "3C2_1"
    # 3C2_2 (y, -x-y, x軸)
    c2_2_axes = [[0, 1, 0], [np.sqrt(3)/2, -1/2, 0], [-np.sqrt(3)/2, -1/2, 0]]
    for axis in c2_2_axes:
        if np.allclose(rot, _rotation_matrix(axis, np.pi), atol=tol):
            return "3C2_2"
    # 2S3 (z軸)
    for angle in [2*np.pi/3, 4*np.pi/3]:
        s3 = _rotation_matrix([0, 0, 1], angle) @ np.diag([1, 1, -1])
        if np.allclose(rot, s3, atol=tol):
            return "2S3"
    # 2S6 (z軸)
    for angle in [np.pi/3, 5*np.pi/3]:
        s6 = _rotation_matrix([0, 0, 1], angle) @ np.diag([1, 1, -1])
        if np.allclose(rot, s6, atol=tol):
            return "2S6"
    # sigmah (xy)
    if np.allclose(rot, np.diag([1, 1, -1]), atol=tol):
        return "sigmah"
    # 3sigmad (d面)
    for angle in [np.pi/6, np.pi/2, 5*np.pi/6]:
        axis = [np.cos(angle), np.sin(angle), 0]
        n = np.array(axis).reshape(3, 1)
        sigma = np.eye(3) - 2 * n @ n.T
        if np.allclose(rot, sigma, atol=tol):
            return "3sigmad"
    # 3sigmav (v面)
    for angle in [0, 2*np.pi/3, 4*np.pi/3]:
        axis = [np.cos(angle), np.sin(angle), 0]
        n = np.array(axis).reshape(3, 1)
        sigma = np.eye(3) - 2 * n @ n.T
        if np.allclose(rot, sigma, atol=tol):
            return "3sigmav"
    return "unknown"

def _assign_op_name_t(rot, tr, tol):
    if np.allclose(rot, np.eye(3), atol=tol):
        return "E"
    mat_c2x = _rotation_matrix([1,0,0], np.pi)
    mat_c2y = _rotation_matrix([0,1,0], np.pi)
    mat_c2z = _rotation_matrix([0,0,1], np.pi)
    mat_c3_1 = _rotation_matrix([1,1,1], 2*np.pi/3)
    mat_c3_2 = _rotation_matrix([-1,1,1], 2*np.pi/3)
    mat_c3_3 = _rotation_matrix([1,-1,1], 2*np.pi/3)
    mat_c3_4 = _rotation_matrix([1,1,-1], 2*np.pi/3)
    if any(np.allclose(rot, m, atol=tol) for m in [mat_c2x, mat_c2y, mat_c2z]):
        return "3C2"
    if any(np.allclose(rot, m, atol=tol) for m in [mat_c3_1, mat_c3_2, mat_c3_3, mat_c3_4]):
        return "4C3"
    return "unknown"

def _assign_op_name_th(rot, tr, tol):
    # E
    if np.allclose(rot, np.eye(3), atol=tol):
        return "E"
    # i
    if np.allclose(rot, -np.eye(3), atol=tol):
        return "i"
    # 8C3 ([111]型)
    c3_axes = [[1, 1, 1], [-1, 1, 1], [1, -1, 1], [1, 1, -1]]
    for axis in c3_axes:
        axis = np.array(axis)/np.linalg.norm(axis)
        for angle in [2*np.pi/3, 4*np.pi/3]:
            if np.allclose(rot, _rotation_matrix(axis, angle), atol=tol):
                return "8C3"
    # 3C2 (x, y, z軸)
    c2_axes = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    for axis in c2_axes:
        if np.allclose(rot, _rotation_matrix(axis, np.pi), atol=tol):
            return "3C2"
    return "unknown"

def _assign_op_name_o(rot, tr, tol):
    # E
    if np.allclose(rot, np.eye(3), atol=tol):
        return "E"
    # 8C3 ([111]型)
    c3_axes = [[1, 1, 1], [-1, 1, 1], [1, -1, 1], [1, 1, -1]]
    for axis in c3_axes:
        axis = np.array(axis)/np.linalg.norm(axis)
        for angle in [2*np.pi/3, 4*np.pi/3]:
            if np.allclose(rot, _rotation_matrix(axis, angle), atol=tol):
                return "8C3"
    # 6C2 (x, y, z, face-diagonal axes)
    c2_axes = [[1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 0], [1, 0, 1], [0, 1, 1], [-1, 1, 0], [-1, 0, 1], [0, -1, 1]]
    for axis in c2_axes:
        axis = np.array(axis)/np.linalg.norm(axis)
        if np.allclose(rot, _rotation_matrix(axis, np.pi), atol=tol):
            return "6C2"
    # 6C4 (x, y, z軸)
    c4_axes = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    for axis in c4_axes:
        for angle in [np.pi/2, 3*np.pi/2]:
            if np.allclose(rot, _rotation_matrix(axis, angle), atol=tol):
                return "6C4"
    return "unknown"

def _assign_op_name_td(rot, tr, tol):
    # E
    if np.allclose(rot, np.eye(3), atol=tol):
        return "E"
    # 8C3 ([111]型)
    c3_axes = [[1, 1, 1], [-1, 1, 1], [1, -1, 1], [1, 1, -1]]
    for axis in c3_axes:
        axis = np.array(axis)/np.linalg.norm(axis)
        for angle in [2*np.pi/3, 4*np.pi/3]:
            if np.allclose(rot, _rotation_matrix(axis, angle), atol=tol):
                return "8C3"
    # 3C2 (x, y, z軸)
    c2_axes = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    for axis in c2_axes:
        if np.allclose(rot, _rotation_matrix(axis, np.pi), atol=tol):
            return "3C2"
    # 6S4 (x, y, z軸)
    for axis in c2_axes:
        for angle in [np.pi/2, 3*np.pi/2]:
            s4 = _rotation_matrix(axis, angle) @ np.diag([1, 1, -1])
            if np.allclose(rot, s4, atol=tol):
                return "6S4"
    # 6sigmad (立方体面心)
    c2p_axes = [[1, 1, 0], [1, 0, 1], [0, 1, 1], [-1, 1, 0], [-1, 0, 1], [0, -1, 1]]
    for axis in c2p_axes:
        axis = np.array(axis)/np.linalg.norm(axis)
        n = axis.reshape(3, 1)
        sigma = np.eye(3) - 2 * n @ n.T
        if np.allclose(rot, sigma, atol=tol):
            return "6sigmad"
    return "unknown"

def _assign_op_name_oh(rot, tr, tol):
    # E
    if np.allclose(rot, np.eye(3), atol=tol):
        return "E"
    # i
    if np.allclose(rot, -np.eye(3), atol=tol):
        return "i"
    # 6C4 (x, y, z軸)
    c4_axes = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    c4_angles = [np.pi/2, 3*np.pi/2]
    for axis in c4_axes:
        for angle in c4_angles:
            if np.allclose(rot, _rotation_matrix(axis, angle), atol=tol):
                return "6C4"
    # 8C3 ([111]型)
    c3_axes = [[1, 1, 1], [-1, 1, 1], [1, -1, 1], [1, 1, -1]]
    for axis in c3_axes:
        axis = np.array(axis)/np.linalg.norm(axis)
        for angle in [2*np.pi/3, 4*np.pi/3]:
            if np.allclose(rot, _rotation_matrix(axis, angle), atol=tol):
                return "8C3"
    # 3C2 (x, y, z軸)
    c2_axes = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    for axis in c2_axes:
        if np.allclose(rot, _rotation_matrix(axis, np.pi), atol=tol):
            return "3C2"
    # 6C2 (立方体面心軸)
    c2p_axes = [[1, 1, 0], [1, 0, 1], [0, 1, 1], [-1, 1, 0], [-1, 0, 1], [0, -1, 1]]
    for axis in c2p_axes:
        axis = np.array(axis)/np.linalg.norm(axis)
        if np.allclose(rot, _rotation_matrix(axis, np.pi), atol=tol):
            return "6C2"
    # 6S4 (x, y, z軸)
    for axis in c4_axes:
        for angle in c4_angles:
            s4 = _rotation_matrix(axis, angle) @ np.diag([1,1,-1])
            if np.allclose(rot, s4, atol=tol):
                return "6S4"
    # 8S6 ([111]型)
    for axis in c3_axes:
        axis = np.array(axis)/np.linalg.norm(axis)
        for angle in [np.pi/3, 5*np.pi/3]:
            s6 = _rotation_matrix(axis, angle) @ np.diag([1,1,-1])
            if np.allclose(rot, s6, atol=tol):
                return "8S6"
    # 3sigmah (xy, yz, zx)
    if np.allclose(rot, np.diag([1,1,-1]), atol=tol):
        return "3sigmah"
    if np.allclose(rot, np.diag([1,-1,1]), atol=tol):
        return "3sigmah"
    if np.allclose(rot, np.diag([-1,1,1]), atol=tol):
        return "3sigmah"
    # 6sigmad (立方体面心)
    for axis in c2p_axes:
        axis = np.array(axis)/np.linalg.norm(axis)
        n = axis.reshape(3,1)
        sigma = np.eye(3) - 2 * n @ n.T
        if np.allclose(rot, sigma, atol=tol):
            return "6sigmad"
    return "unknown"

def _rotation_matrix(axis, theta):
    axis = np.array(axis, dtype=float)
    axis = axis / np.linalg.norm(axis)
    a = np.cos(theta / 2.0)
    b, c, d = -axis * np.sin(theta / 2.0)
    return np.array([
        [a * a + b * b - c * c - d * d, 2 * (b * c + a * d), 2 * (b * d - a * c)],
        [2 * (b * c - a * d), a * a + c * c - b * b - d * d, 2 * (c * d + a * b)],
        [2 * (b * d + a * c), 2 * (c * d - a * b), a * a + d * d - b * b - c * c],
    ])

