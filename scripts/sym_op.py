import numpy as np

def apply_sym_for_all(pos: np.ndarray, rot: np.ndarray, trans: np.ndarray): ## nsite, x, y, z
    return (np.dot(rot, pos.T) + trans).T


def apply_sym_single(pos: np.ndarray, rot: np.ndarray, trans: np.ndarray): ## nsite, x, y, z
    return np.dot(rot, pos) + trans


def apply_for_orb(amps_dict: dict, rot, trs): ### dic[(x,y,z)] = amp
    transformed = {}
    for k, v in amps_dict.items():
        x, y, z = k
        amp = v
        coord = np.array([x, y, z])
        new_coord = apply_sym_single(coord, rot, trs)
        x, y, z = new_coord
        transformed[(x, y, z)] = amp
    return transformed

def check_symmetry(amps_dict, transformed_dict, tol=1e-6):
    results = []
    for k, amp1 in amps_dict.items():
        x, y, z = k
        amp2 = transformed_dict[(x, y, z)]
        if abs(amp1 - amp2) < tol:
            results.append(1)
        elif abs(amp1 + amp2) < tol:
            results.append(-1)
        else:
            results.append(0)
    unique = set(results)
    if unique == {1}:
        return 1
    elif unique == {-1}:
        return -1
    else:
        return 0

