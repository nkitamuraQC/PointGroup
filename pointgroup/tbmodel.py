from .pythtb_respack import w90, tb_model
from seekpath import get_path
import numpy as np

class TBModel:
    """
    タイトバインディングモデル（Wannier90/pythtb/respack形式）をラップし、
    格子・軌道・バンド情報や対称k点経路などの生成・管理を行うクラス。
    """
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
        """
        Parameters
        ----------
        wannier_path : str, optional
            wannier90出力ディレクトリのパス
        wannier_prefix : str, optional
            wannier90ファイルのプレフィックス
        wannier_plot : str, optional
            wannier90 plotファイル名
        nk : list, optional
            k点分割数（デフォルト: [4,1,1]）
        site_species : list, optional
            サイト種別リスト
        site_nlm : list, optional
            サイトごとの軌道情報
        pythtb_obj : tb_model, optional
            既存のpythtbオブジェクト
        """
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
            assert(pythtb_obj._dim_r == 3)

        self.path_sym = None
        self.path_coord = None
        self.kpts = None

    def gen_pythtb(self):
        """
        pythtbオブジェクトを生成し、格子・軌道・k点・バンド情報を初期化する。
        Returns
        -------
        eval : ndarray
            バンドエネルギー配列
        """
        if self.pythtb_obj is None:
           self.pythtb_obj = w90(self.path, self.prefix).model()
        self.cell = self.pythtb_obj._lat
        self.pos = self.pythtb_obj._orb
        # self.pos_with_corner = self._fill_corners()
        self._get_uniform_kpts()
        eval = self._get_bands()
        return eval

    def _fill_corners(self):
        """
        原点[0,0,0]が存在する場合、格子ベクトル平行移動分のコーナー座標もself.posに追加する。
        """
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
        """
        seekpathを用いて対称k点経路・座標を計算し、self.path, self.path_coord等に格納する。
        """
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
        """
        一様k点メッシュをpythtbから取得しself.kptsに格納する。
        """
        self.kpts = self.pythtb_obj.k_uniform_mesh(self.nk)
        return
    
    def _get_bands(self):
        """
        全k点でバンド計算を行い、固有ベクトル（波動関数）をself.Umatに格納する。
        Returns
        -------
        eval : ndarray
            バンドエネルギー配列
        """
        (eval, evec) = self.pythtb_obj.solve_all(self.kpts, eig_vectors=True)
        self.Umat = evec
        return eval
    

