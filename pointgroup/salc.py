from .get_character import GetCharacter
from .schonflies import point_group_map
from .pglib import point_group, get_rep
from .sym_op import (
    apply_for_orb
)
import numpy as np

class GetSALC:
    def __init__(self, getc: GetCharacter):
        self.getc = getc
        self.pg = point_group_map[getc.pg]
        self.ir_ch = point_group[self.pg] ## 点群
        self.pos = getc.pos
        self.translate = getc.trans
        self.rot = getc.rot
        self.rot_o3 = getc.rot_o3

    def _gen_site_amps(self):
        results = []
        for i, coord in self.pos:
            amp_site = {}
            x, y, z = coord
            amp_site[(x, y, z)] = 1
            results.append(amp_site)
        return results

    def _search_op(self, op):
        reps = get_rep()
        for rep_name, rep in reps:
            diff = np.ravel(rep) - np.ravel(op)
            if np.linalg.norm(diff) < 1e-6:
                return rep_name
        return None
    
    def make_salc_site(self, sym_idx):
        ## s軌道的な対称性のみ対象;
        ## より精密には軌道のgrid情報必要
        rot = self.rot[sym_idx]
        rot_o3 = self.rot_o3[sym_idx]
        trs = self.translate[sym_idx]
        site_amps = self._gen_site_amps()
        weights = {}
        for samp in site_amps:
            for xyz in samp:
                x, y, z = xyz
                weights[(x, y, z)] = 0
        op_name = self._search_op(rot_o3)
        op_name = op_name.split(" ")[0]
        for samp in site_amps:
            applied = apply_for_orb(samp, rot, trs)
            for xyz, amp in applied.items():
                x, y, z = xyz
                weights[(x, y, z)] += amp * self.ir_ch[self.pg][op_name]
        return weights

    def reconst_umat(self):
        raise NotImplementedError

    def get_ham(self):
        raise NotImplementedError





