from .get_character import GetCharacter
from .schonflies import schonflies
from .point_groups import point_group
from .sym_op import (
    apply_for_orb
)
import numpy as np

class GetSALC:
    def __init__(self, getc: GetCharacter):
        self.getc = getc
        self.pg = schonflies[getc.pg]
        self.ir_ch = point_group[self.pg] ## 点群
        self.pos = getc.pos
        self.translate = getc.trans
        self.rot = getc.rot

    def gen_site_amps(self):
        results = []
        for i, coord in self.pos:
            amp_site = {}
            x, y, z = coord
            amp_site[(x, y, z)] = 1
            results.append(amp_site)
        return results
    
    def make_salc_site(self, sym_idx):
        ## s軌道的な対称性のみ対象;
        ## より精密には軌道のgrid情報必要
        rot = self.rot[sym_idx]
        trs = self.translate[sym_idx]
        site_amps = self.gen_site_amps()
        weights = {}
        for samp in site_amps:
            for xyz in samp:
                x, y, z = xyz
                weights[(x, y, z)] = 0
        for samp in site_amps:
            applied = apply_for_orb(samp, rot, trs)
            for xyz, amp in applied.items():
                x, y, z = xyz
                weights[(x, y, z)] += amp * self.ir_ch[sym_idx]
        return weights

    def get_ham(self):
        raise NotImplementedError





