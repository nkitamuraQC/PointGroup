from .tbmodel import TBModel

class GetCharacter:
    def __init__(self, tbmodel: TBModel):
        self.tbmodel = tbmodel
        self.spg = None
        self.pg = None
        self.symm_ops = None

    def gen_space_group(self):
        return
    
    def gen_point_group(self):
        return
    
    def get_symm_ops(self):
        return
    
    def get_character(self, kidx, orb_idx, symm_idx):
        return
    
