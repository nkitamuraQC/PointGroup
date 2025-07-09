from .pythtb_respack import w90
from seekpath import get_path

class TBModel:
    def __init__(
            self, 
            wannier_path: str, 
            wannier_prefix: str, 
            wannier_plot: str,
            nk: list = [4, 1, 1],
            site_species: list = [0]
        ):
        self.path = wannier_path
        self.prefix = wannier_prefix
        self.wannier_plot = wannier_plot
        self.nk = nk
        self.site_species = site_species
        self.pythtb_obj = None
        self.cell = None
        self.pos = None
        self.Umat = None
        self.path_sym = None
        self.path_coord = None
        self.kpts = None

    def gen_pythtb(self):
        self.pythtb_obj = w90(self.path, self.prefix).model()
        self.cell = self.pythtb_obj._lat
        self.pos = self.pythtb_obj._orb
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
    
    def get_uniform_kpts(self):
        self.kpts = self.pythtb_obj.k_uniform_mesh(self, self.nk)
        return
    
    def get_bands(self):
        (eval, evec) = self.pythtb_obj.solve_all(self.kpts, eig_vectors=True)
        self.Umat = evec
        return eval
    

