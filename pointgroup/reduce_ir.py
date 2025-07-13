from .get_character import GetCharacter
from .schonflies import point_group_map
from .pglib import point_group
from .sym_op import (
    apply_for_orb
)
import numpy as np

class GetIR:
    def __init__(self, getc: GetCharacter):
        self.getc = getc
        self.pg = point_group_map[getc.pg]
        self.ir_ch_all = point_group[self.pg] ## 点群
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
    
    def _get_ir_ch(self, ir_idx):
        ir_ch = []
        for sym_op, ch in self.ir_ch_all.items():
            try:
                nx = int(sym_op[0])
            except ValueError:
                nx = 1
            for _ in range(nx):
                ir_ch.append(ch[ir_idx])
        return np.array(ir_ch)


    def reduce_ir(self, ch: np.ndarray, ir_idx=0):
        self.h = self._get_h()
        self.ir_ch = self._get_ir_ch(ir_idx)
        
        return np.dot(ch, self.ir_ch) / self.h

        

    