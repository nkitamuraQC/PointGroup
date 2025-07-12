from .pythtb_respack import w90, tb_model
from seekpath import get_path
import numpy as np

class TBModel:
    def __init__(
            self, 
            wannier_path: str = None, 
            wannier_prefix: str = None, 
            wannier_plot: str = None,
            nk: list = [4, 1, 1],
            site_species: list = [0],
            site_nlm: list = [(1,0,0)],
            pythtb_obj: tb_model = None
        ):
        self.path = wannier_path
        self.prefix = wannier_prefix
        self.wannier_plot = wannier_plot
        self.nk = nk
        self.site_species = site_species
        self.site_nlm = site_nlm
        self.pythtb_obj = pythtb_obj
        self.cell = None
        self.pos = None
        self.Umat = None
        if pythtb_obj is not None:
            assert(pythtb_obj._dim_r < 3)

        self.path_sym = None
        self.path_coord = None
        self.kpts = None

    def gen_pythtb(self):
        if self.pythtb_obj is None:
           self.pythtb_obj = w90(self.path, self.prefix).model()
        self.cell = self.pythtb_obj._lat
        self.pos = self.pythtb_obj._orb
        # self.pos_with_corner = self._fill_corners()
        self._get_uniform_kpts()
        eval = self._get_bands()
        return eval

    def _fill_corners(self):
        pos_arr = np.zeros((3))
        flag = False
        for pos in self.pos:
            if pos == [0, 0, 0]:
                flag = True
        if flag:
            self.pos.append((pos_arr+self.cell[0]).tolist())
            self.pos.append((pos_arr+self.cell[1]).tolist())
            self.pos.append((pos_arr+self.cell[2]).tolist())
            self.pos.append((pos_arr+self.cell[0]+self.cell[1]).tolist())
            self.pos.append((pos_arr+self.cell[0]+self.cell[1]+self.cell[2]).tolist())
            
        return
    
    def get_sym_kpts(self):
        cell = self.cell
        pos = self.pos
        n = [i for i in range(len(pos))]
        struct = (cell, pos, n)
        res = get_path(struct)
        self.coord = res["point_coords"]
        self.path = res['path']
        
        path_sym = []
        for i, p in enumerate(self.path):
            path_sym.append(p[0])
            if i == len(self.path) - 1:
                path_sym.append(p[1])

        path_coord = []
        for p in path_sym:
            path_coord.append(self.coord[p])

        self.path_coord = path_coord
        self.path_sym = path_sym
        return 
    
    def _get_uniform_kpts(self):
        self.kpts = self.pythtb_obj.k_uniform_mesh(self, self.nk)
        return
    
    def _get_bands(self):
        (eval, evec) = self.pythtb_obj.solve_all(self.kpts, eig_vectors=True)
        self.Umat = evec
        return eval
    

