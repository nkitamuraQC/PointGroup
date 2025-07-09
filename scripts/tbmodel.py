

class TBModel:
    def __init__(self, w90_ints_file: str, wannier_plot_file: str):
        self.w90_ints_file = w90_ints_file
        self.wannier_plot_file = wannier_plot_file
        self.kpts = [4, 1, 1]
        self.pythtb_obj = None
        self.cell = None
        self.pos = None
        self.Umat = None

    def gen_pythtb(self):
        return
    
    def get_bands(self):
        return
    
    
