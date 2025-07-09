from .tbmodel import TBModel
import spglib
from .sym_op import (
    apply_sym_for_all,
    apply_sym_single,
    apply_for_orb,
    check_symmetry,
)
import numpy as np

class GetCharacter:
    def __init__(self, tbmodel: TBModel):
        self.tbmodel = tbmodel
        self.spg = None
        self.pg = None
        self.symm_ops = None
        self.rot = None
        self.trans = None
        self.Umat = tbmodel.Umat ## site, k, band
        self.dataset = None
        self.pos = self.tbmodel.pos ## nsite, ndim
        self.site_species = self.tbmodel.site_species
        self.site_nlm = self.tbmodel.site_nlm

    def gen_space_group(self): # international
        cell = (self.tbmodel.cell, self.tbmodel.pos, self.tbmodel.site_species)
        if self.dataset is None:
            self.dataset = spglib.spglib.get_symmetry_dataset(cell)
        self.spg = self.dataset["international"]
        return
    
    def gen_point_group(self):
        cell = (self.tbmodel.cell, self.tbmodel.pos, self.tbmodel.site_species)
        if self.dataset is None:
            self.dataset = spglib.spglib.get_symmetry_dataset(cell)
        self.pg = self.dataset["pointgroup"]
        return
    
    def get_symm_ops(self):
        cell = (self.tbmodel.cell, self.tbmodel.pos, self.tbmodel.site_species)
        self.symm_ops = spglib.spglib.get_symmetry(cell)
        self.rot = self.symm_ops["rotations"]
        self.trans = self.symm_ops["translations"]
        return
    
    def get_site_amps(self, kidx, orb_idx):
        site_amps = {}
        umat = self.Umat[:, kidx, orb_idx] ## nsite
        pos = self.tbmodel.pos ## nsite, (x, y, z)
        nsites = pos.shape[0]
        for i in range(nsites):
            x, y, z = pos[i]
            site_amps[(x, y, z)] = umat[i]
        return site_amps


    def get_grid_amps(self, kidx, orb_idx):
        raise NotImplementedError

    def get_ch_occnum_vector(self):
        return
    
    def get_character(self, kidx, orb_idx, symm_idx, mode="site"):
        if mode == "site":
            amps = self.get_site_amps(kidx, orb_idx)
        elif mode == "grid":
            amps = self.get_grid_amps(kidx, orb_idx)
        else:
            raise NotImplementedError
        
        rot = self.rot[symm_idx]
        trs = self.trans[symm_idx]
        if mode == "site":
            result = apply_for_orb(amps, rot, trs)
            charcter = check_symmetry(amps, result)
            return charcter
            

        elif mode == "grid":
            raise NotImplementedError
        else:
            raise NotImplementedError
    
