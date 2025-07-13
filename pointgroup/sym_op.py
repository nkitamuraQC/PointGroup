import numpy as np
from itertools import product
import copy

def apply_sym_for_all(
    pos: np.ndarray, rot: np.ndarray, trans: np.ndarray
):  ## nsite, x, y, z
    return (np.dot(rot, pos.T) + trans).T


def apply_sym_single(
    pos: np.ndarray, rot: np.ndarray, trans: np.ndarray
):  ## nsite, x, y, z
    return np.dot(rot, pos) + trans


def apply_for_orb(amps_dict: dict, rot, trs):  ### dic[(x,y,z)] = amp
    transformed = {}
    for k, v in amps_dict.items():
        x, y, z = k
        amp = v
        coord = np.array([x, y, z])
        new_coord = apply_sym_single(coord, rot, trs)
        x, y, z = new_coord
        transformed[(x, y, z)] = amp
    return transformed

def get_diff(pos, cell):
    deltas = [-2, -1, 0, 1, 2]
    offsets = list(product(deltas, repeat=3))
    res = []
    pos_save = copy.deepcopy(pos)
    for dx, dy, dz in offsets:
        pos += np.array([dx, dy, dz])
        res.append(pos)
        pos = pos_save.copy()
    return res
        


def check_symmetry(amps_dict, transformed_dict, cell=None, tol=1e-6):
    flag = None
    for k1, amp1 in amps_dict.items():
        k1 = np.array(k1)
        for k2, amp2 in transformed_dict.items():
            k2 = np.array(k2)
            all_pos = get_diff(k2.copy(), cell)
            for pos in all_pos:
                # print(k1, pos)
                if np.linalg.norm(k1 - pos) < tol:
                    if abs(amp1 + amp2) < tol:
                        flag = -1
                        break
                    elif abs(amp1 - amp2) < tol:
                        flag = 1
                        break
                    else:
                        flag = 0
                        break

    return flag
        