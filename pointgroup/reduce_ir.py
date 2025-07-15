from .get_character import GetCharacter
from .schonflies import point_group_map
from .pglib import point_group, find_operation_type, output_sym_op_name
from .sym_op import apply_for_orb
import numpy as np
import copy

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

    def _get_ir_ch_old(self, ir_idx):
        ir_ch = []
        rot_o3 = self.getc.rot_o3
        for i in range(rot_o3.shape[0]):
            op_name, disc = find_operation_type(rot_o3[i])
            op_name_save = copy.copy(op_name)
            for k, v in self.ir_ch_all.items():
                try:
                    int(k[0])
                    k = k[1:]
                except ValueError:
                    pass
                if "(" in op_name and "(" not in k:
                    op_name_save = op_name_save.split("(")[0]
                if "rotation around principal" in disc:
                    op_name = op_name_save
                    op_name += "_1"
                elif "(perpendicular to principal axis)" in disc:
                    op_name = op_name_save
                    op_name += "_2"
                else:
                    op_name = op_name_save
                if k == op_name:
                    # print(op_name, k)
                    ch = v
                    ir_ch.append(ch[ir_idx])
                    break
        return ir_ch
    
    def _get_ir_ch_naiive(self, ir_idx):
        ir_ch = []
        for ir, ch in self.ir_ch_all.items():
            try:
                n = int(ir[0])
            except ValueError:
                n = 1
            for i in range(n):
                ir_ch.append(ch[ir_idx])
        
        return np.array(ir_ch)
    
    def _get_ir_ch(self, ir_idx):
        ir_ch = []
        rot_o3 = self.getc.rot_o3
        for i in range(rot_o3.shape[0]):
            op_name, disc = find_operation_type(rot_o3[i])
            op_name2 = output_sym_op_name(self.pg, op_name)
            for k, ch in self.ir_ch_all.items():
                if k == op_name2:
                    print(k, op_name, op_name2, disc, ch[ir_idx], ir_idx)
                    ir_ch.append(ch[ir_idx])

        return np.array(ir_ch)

    def reduce_ir(self, ch: np.ndarray, ir_idx=0):
        h = self._get_h()
        ir_ch = self._get_ir_ch(ir_idx)
        # print(ch, self.ir_ch)
        print(ch, ir_ch)

        return np.dot(ch, ir_ch) / h
