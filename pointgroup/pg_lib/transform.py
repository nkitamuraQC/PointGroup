### O(3)群に変換: L*R*L^-1
import numpy as np

def tr2o3(cell, R):
    L = cell
    L_inv = np.linalg.inv(L)
    return np.dot(np.dot(L, R), L_inv)

def fromo3(cell, R):
    L = cell
    L_inv = np.linalg.inv(L)
    return np.dot(np.dot(L_inv, R), L)

