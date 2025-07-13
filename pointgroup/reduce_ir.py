from .get_character import GetCharacter
from .schonflies import point_group_map
from .pglib import point_group, assign_op_name
from .sym_op import apply_for_orb
import numpy as np


class GetIR:
    def __init__(self, getc: GetCharacter):
        self.getc = getc
        self.pg = point_group_map[getc.pg]
        self.ir_ch_all = point_group[self.pg]  ## 点群
        self.ir_ch = None
        self.h = 0

    def _get_h(self):
        h = 0
        for sym_op, ch in self.ir_ch_all.items():
            try:
                x = int(sym_op[0])
            except ValueError:
                x = 1
            h += x
        return h
    
    def _get_ir(self):
        return len(self.ir_ch_all["E"])

    def _get_ir_ch(self, ir_idx):
        op_names = self._sort_op()
        ir_ch = []
        # print(op_names)
        for sym_op in op_names:
            ch = self.ir_ch_all[sym_op]
            ir_ch.append(ch[ir_idx])
        return np.array(ir_ch)

    def reduce_ir(self, ch: np.ndarray, ir_idx=0):
        self.h = self._get_h()
        self.ir_ch = self._get_ir_ch(ir_idx)
        # print(ch, self.ir_ch)
        return np.dot(ch, self.ir_ch) / self.h
    
    def _search_op2(self, op):
        """
        与えられたO(3)行列opに最も近い点群操作名を検索する。
        Returns
        -------
        op_name : str or None
            一致する操作名（なければNone）
        """
        reps = get_rep()[self.pg]
        delta_list = []
        rep_names = []
        for rep_name, rep in reps.items():
            diff = np.linalg.norm(np.ravel(rep) - np.ravel(op))
            print(op, rep)
            delta_list.append(diff)
            rep_names.append(rep_name.split()[0])
        delta = np.array(delta_list)
        idx = np.argmin(delta)
        return rep_names[idx]
    
    def _sort_op(self):
        rot = self.getc.rot
        rot_o3 = self.getc.rot_o3
        # print(rot_o3)
        nops = rot.shape[0]
        op_names = []
        for n in range(nops):
            op_name = self._search_op2(rot[n])
            op_names.append(op_name)
                
        return op_names
