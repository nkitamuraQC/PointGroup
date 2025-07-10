from .get_character import GetCharacter
from .schonflies import schonflies
from .point_groups import point_group
from .sym_op import (
    apply_for_orb
)
import numpy as np

class GetIR:
    def __init__(self, getc: GetCharacter):
        self.getc = getc
        self.pg = schonflies[getc.pg]
        self.ir_ch = point_group[self.pg] ## 点群
        self.h = len(self.ir_ch[0])

    def reduce_ir(self, ch, ir_idx):
        return np.dot(ch, self.ir_ch[ir_idx]) / self.h

        

    