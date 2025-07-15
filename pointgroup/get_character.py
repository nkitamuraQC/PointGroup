from .tbmodel import TBModel
import spglib
from .sym_op import (
    apply_sym_for_all,
    apply_sym_single,
    apply_for_orb,
    check_symmetry,
    check_symmetry2
)
import numpy as np
from .pglib import tr2o3
from .schonflies import point_group_map
from .utils import logger
from .pglib import find_operation_type

class GetCharacter:
    """
    Tight-bindingモデルから点群・空間群のキャラクター（指標）を計算するためのユーティリティクラス。
    spglibによる対称操作取得、軌道・サイトごとのアンプリチュード取得、
    およびキャラクター計算のための補助関数を提供する。
    """

    def __init__(self, tbmodel: TBModel):
        """
        Parameters
        ----------
        tbmodel : TBModel
            タイトバインディングモデルのインスタンス。
        """
        self.tbmodel = tbmodel
        self.spg = None
        self.pg = None
        self.symm_ops = None
        self.rot = None
        self.rot_o3 = None
        self.trans = None
        self.Umat = tbmodel.Umat  # site, k, band
        self.dataset = None
        self.nsymm_ops = 0
        self.pos = self.tbmodel.pos  # nsite, ndim
        self.tbcell = self.tbmodel.cell
        self.site_species = self.tbmodel.site_species
        self.site_nlm = self.tbmodel.site_nlm
        self.use_trace = True

    def _gen_space_group(self):
        """
        spglibを用いて空間群（international symbol）を取得し、self.spgに格納する。
        """
        cell = (self.tbmodel.cell, self.tbmodel.pos, self.tbmodel.site_species)
        if self.dataset is None:
            self.dataset = spglib.spglib.get_symmetry_dataset(cell)
        self.spg = self.dataset["international"]
        return

    def _gen_point_group(self):
        """
        spglibを用いて点群（point group symbol）を取得し、self.pgに格納する。
        """
        cell = (self.tbmodel.cell, self.tbmodel.pos, self.tbmodel.site_species)
        if self.dataset is None:
            self.dataset = spglib.spglib.get_symmetry_dataset(cell)
        self.pg = self.dataset["pointgroup"]
        return

    def _get_symm_ops(self):
        """
        spglibを用いて全対称操作（回転・並進）を取得し、
        O(3)回転行列へ変換してself.rot_o3に格納する。
        """
        cell = (self.tbmodel.cell, self.tbmodel.pos, np.array(self.site_species))
        self.symm_ops = spglib.spglib.get_symmetry(cell)
        self.rot = self.symm_ops["rotations"]
        self.trans = self.symm_ops["translations"]

        self.rot_o3 = np.zeros_like(self.rot)
        for i in range(self.rot_o3.shape[0]):
            R = self.rot[i, :, :]
            self.rot_o3[i, :, :] = tr2o3(self.tbcell, R)
        return self.rot.shape[0]

    def _get_site_amps(self, kidx, orb_idx):
        """
        指定したk点・バンドのサイトごとのアンプリチュード（波動関数成分）を取得する。
        原点サイトの場合は格子ベクトル平行移動分も登録する。

        Parameters
        ----------
        kidx : int
            k点インデックス
        orb_idx : int
            バンド（軌道）インデックス

        Returns
        -------
        site_amps : dict
            {(x, y, z): amp} の辞書
        """
        site_amps = {}
        umat = self.Umat[:, kidx, orb_idx]  # nsite
        pos = self.tbmodel.pos  # nsite, (x, y, z)
        nsites = pos.shape[0]
        for i in range(nsites):
            x, y, z = pos[i]
            site_amps[(x, y, z)] = umat[i]
            # current_pos = np.array([x, y, z])
            # if np.linalg.norm(current_pos) < 1e-6:
            #     x, y, z = current_pos+self.tbcell[0]
            #     site_amps[(x, y, z)] = umat[i]
            #     x, y, z = current_pos+self.tbcell[1]
            #     site_amps[(x, y, z)] = umat[i]
            #     x, y, z = current_pos+self.tbcell[2]
            #     site_amps[(x, y, z)] = umat[i]
            #     x, y, z = current_pos+self.tbcell[0]+self.tbcell[1]
            #     site_amps[(x, y, z)] = umat[i]
            #     # 3つの格子ベクトルの和
            #     x, y, z = current_pos + self.tbcell[0] + self.tbcell[1]
            #     site_amps[(x, y, z)] = umat[i]
            #     x, y, z = current_pos + self.tbcell[0] + self.tbcell[2]
            #     site_amps[(x, y, z)] = umat[i]
            #     x, y, z = current_pos + self.tbcell[1] + self.tbcell[2]
            #     site_amps[(x, y, z)] = umat[i]
            #     # 3つの格子ベクトルの和（1行が長くなるのを回避）
            #     vec = self.tbcell[0] + self.tbcell[1] + self.tbcell[2]
            #     x, y, z = current_pos + vec
            #     site_amps[(x, y, z)] = umat[i]
        return site_amps

    def _get_grid_amps(self, kidx, orb_idx):
        """
        グリッド上のアンプリチュード（波動関数成分）を取得する（未実装）。
        """
        raise NotImplementedError

    def get_ch_occnum_vector(self):
        """
        キャラクターの占有数ベクトルを返す（未実装）。
        """
        raise NotImplementedError

    def get_character(self, kidx=0, orb_idx=0, op_idx=0, mode="site"):
        """
        指定したk点・バンド・対称操作に対するキャラクター（指標）を計算する。

        Parameters
        ----------
        kidx : int
            k点インデックス
        orb_idx : int
            バンド（軌道）インデックス
        symm_idx : int
            対称操作インデックス
        mode : str
            'site'または'grid'（デフォルト: 'site'）

        Returns
        -------
        character : float or complex
            キャラクター値
        """
        self._gen_point_group()
        self._gen_space_group()
        self._get_symm_ops()
        if mode == "site":
            amps = self._get_site_amps(kidx, orb_idx)
        elif mode == "grid":
            amps = self._get_grid_amps(kidx, orb_idx)
        else:
            raise NotImplementedError

        rot = self.rot[op_idx]
        trs = self.trans[op_idx]
        try:
            logger.info(f"pg: {point_group_map[self.pg]}")
        except KeyError:
            raise KeyError("Point group is not Found.")
        logger.info(f"rot: {rot}")
        logger.info(f"trs: {trs}")
        if mode == "site":
            result = apply_for_orb(amps, rot, trs)
            logger.info(f"before: {amps}")
            logger.info(f"after: {result}")
            if self.use_trace:
                op_name, _ = find_operation_type(self.rot_o3[op_idx])
                charcter = check_symmetry2(amps, result)
                logger.info(f"op_name: {op_name}, character: {charcter}")
            else:
                charcter = check_symmetry(amps, result)
                logger.info(f"character: {charcter}")
            return charcter

        elif mode == "grid":
            raise NotImplementedError
        else:
            raise NotImplementedError
