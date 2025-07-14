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
    

def _assign_op_name_o(rot, tr, tol):
    dic = {}
    dic["E E"] = np.diag([1, 1, 1])

    dic["6C4 6C4_1"] = np.array([[0, 0, -1], [0, 1, 0], [1, 0, 0]])
    dic["6C2 6C2_1"] = np.diag([-1, 1, -1])
    dic["6C4 6C4_2"] = np.array([[0, 0, 1], [0, 1, 0], [-1, 0, 0]])
    dic["6C4 6C4_3"] = np.array([[1, 0, 0], [0, 0, -1], [0, 1, 0]])
    dic["6C2 6C2_2"] = np.diag([1, -1, -1])
    dic["6C4 6C4_4"] = np.array([[1, 0, 0], [0, 0, 1], [0, -1, 0]])
    dic["6C4 6C4_5"] = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
    dic["6C2 6C2_3"] = np.diag([-1, -1, 1])
    dic["6C4 6C4_6"] = np.array([[0, 1, 0], [-1, 0, 0], [0, 0, 1]])
    dic["6C2 6C2_3"] = np.array([[-1, 0, 0], [0, 0, 1], [0, 1, 0]])
    ###
    dic[""] = np.array([[-1, 0, 0], [0, 0, -1], [0, -1, 0]])
    dic["6C2 6C2_4"] = np.array([[0, -1, 0], [-1, 0, 0], [0, 0, -1]])
    dic[""] = np.array([[0, 1, 0], [1, 0, 0], [0, 0, -1]])
    dic[""] = np.array([[0, 0, 1], [0, -1, 0], [1, 0, 0]])
    dic[""] = np.array([[0, 0, -1], [0, -1, 0], [-1, 0, 0]])
    dic[""] = np.array([[0, -1, 0], [0, 0, 1], [-1, 0, 0]])
    dic[""] = np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]])
    dic[""] = np.array([[0, 1, 0], [0, 0, -1], [-1, 0, 0]])
    dic[""] = np.array([[0, -1, 0], [0, 0, -1], [1, 0, 0]])
    dic[""] = np.array([[0, 0, -1], [-1, 0, 0], [1, 0, 0]])
    dic[""] = np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]])
    dic[""] = np.array([[0, 0, -1], [1, 0, 0], [0, -1, 0]])
    dic[""] = np.array([[0, 0, 1], [-1, 0, 0], [0, -1, 0]])

    



def _assign_op_name_oh(rot, tr, tol):
    dic = {}
    dic["3sigmah sigmaz"] = np.diag([1, 1, -1])
    dic["3sigmah sigmax"] = np.diag([-1, 1, 1])
    dic["6C4 C4_1"] = np.array([[0, 0, 1], [0, 1, 0], [1, 0, 0]])
    #####
    dic["4S4 S4_3"] = np.array([[0, 0, -1], [0, 1, 0], [-1, 0, 0]])
    dic["6C4 C4_2"] = np.array([[1, 0, 0], [0, 0, 1], [0, 1, 0]])
    dic["3sigmah sigmay"] = np.diag([1, -1, 1])
    dic["6C2 C2_1"] = np.array([[1, 0, 0], [0, 0, -1], [0, -1, 0]])
    dic["6S4 S4_1"] = np.array([[0, -1, 0], [1, 0, 0], [0, 0, -1]])
    dic["i I"] = np.array([[-1, 0, 0], [0, -1, 0], [0, 0, -1]])
    dic["6S4 S4_2"] = np.array([[0, 1, 0], [-1, 0, 0], [0, 0, -1]])
    dic["6sigmad sigmad1"] = np.array([[-1, 0, 0], [0, 0, -1], [0, 1, 0]])
    dic["6sigmad sigmad2"] = np.array([[-1, 0, 0], [0, 0, 1], [0, -1, 0]])
    dic["6sigmad sigmad3"] = np.array([[0, -1, 0], [-1, 0, 0], [0, 0, 1]])
    dic["6sigmad sigmad4"] = np.array([[0, 1, 0], [1, 0, 0], [0, 0, 1]])

    ## S4
    dic["6S4 S4_3"] = np.array([[0, 0, -1], [0, -1, 0], [1, 0, 0]])
    dic["6S4 S4_4"] = np.array([[0, 0, 1], [0, -1, 0], [-1, 0, 0]])

    ### C3
    dic["8S6 S6_1"] = np.array([[0, -1, 0], [0, 0, -1], [-1, 0, 0]])
    dic["6S4 S4_5"] = np.array([[0, 1, 0], [0, 0, -1], [1, 0, 0]])
    dic["8S6 S6_2"] = np.array([[0, 1, 0], [0, 0, 1], [-1, 0, 0]])
    dic["6C4 C4_3"] = np.array([[0, -1, 0], [0, 0, 1], [1, 0, 0]])
    dic["8S6 S6_3"] = np.array([[0, 0, 1], [-1, 0, 0], [0, 1, 0]])
    dic["8S6 S6_4"] = np.array([[0, 0, -1], [1, 0, 0], [0, 1, 0]])
    dic["8S6 S6_5"] = np.array([[0, 0, 1], [1, 0, 0], [0, -1, 0]])
    dic["8S6 S6_6"] = np.array([[0, 0, -1], [-1, 0, 0], [0, -1, 0]])


