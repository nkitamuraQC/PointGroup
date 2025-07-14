from .get_character import GetCharacter
from .schonflies import point_group_map
from .pglib import point_group, assign_op_name, is_similar_symmetry_name
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
        ir_ch = []
        rot_o3 = self.getc.rot_o3
        for i in range(rot_o3.shape[0]):
            name = assign_op_name(rot_o3[i], self.pg)
            print(name, rot_o3[i])
            for k, v in self.ir_ch_all.items():
                # print(name, k, is_similar_symmetry_name(k, name))
                if is_similar_symmetry_name(k, name):
                    ch = self.ir_ch_all[name]
                    ir_ch.append(ch[ir_idx])
        return ir_ch

    def reduce_ir(self, ch: np.ndarray, ir_idx=0):
        self.h = self._get_h()
        self.ir_ch = self._get_ir_ch(ir_idx)
        # print(ch, self.ir_ch)
        return np.dot(ch, self.ir_ch) / self.h
