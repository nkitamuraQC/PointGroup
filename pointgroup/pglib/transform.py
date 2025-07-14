### O(3)群に変換: L*R*L^-1
import numpy as np


def tr2o3(cell, R):
    L = cell
    L_inv = np.linalg.inv(L)
    return np.dot(np.dot(L, R), L_inv)


def from_o3(cell, R):
    L = cell
    L_inv = np.linalg.inv(L)
    return np.dot(np.dot(L_inv, R), L)


import re

def normalize_symmetry_name(name):
    # 大文字小文字無視
    name = name.lower()
    # 不要な文字除去
    name = re.sub(r'[^\w\d\(\)\-\s]', '', name)
    # 連続空白は1つに
    name = re.sub(r'\s+', ' ', name).strip()
    # 操作名と軸を分ける（例: c4 (rotation around (001)) → c4 (001)）
    m = re.match(r'(e|c\d+|s\d+|i|σ|sigma|sigmah|sigmad)(.*)', name)
    if not m:
        return name
    op = m.group(1)
    rest = m.group(2)

    # 軸・面を取り出す
    axis_match = re.search(r'\(([\d\-xyz\s]+)\)', rest)
    if axis_match:
        axis = axis_match.group(1).replace(' ', '')
    else:
        # 軸なしは空文字
        axis = ''

    # σやsigmaを統一
    if op in ['σ', 'sigma', 'sigmah', 'sigmad']:
        op = 'sigma'

    return f"{op}({axis})"


def is_similar_symmetry_name(name1, name2):
    n1 = normalize_symmetry_name(name1)
    n2 = normalize_symmetry_name(name2)
    return n1 == n2
