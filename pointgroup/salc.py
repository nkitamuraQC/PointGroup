from .get_character import GetCharacter
from .schonflies import point_group_map
from .pglib import point_group
from .sym_op import apply_for_orb, get_diff
import numpy as np


class GetSALC:
    """
    SALC（Symmetry Adapted Linear Combination, 対称化線形結合）を計算するためのクラス。
    GetCharacterインスタンスから点群情報・座標・対称操作を受け取り、
    サイトごとのSALCや重み計算を行う。
    """

    def __init__(self, getc: GetCharacter):
        """
        Parameters
        ----------
        getc : GetCharacter
            キャラクター計算用インスタンス
        """
        self.getc = getc
        self.pg = point_group_map[getc.pg]
        self.ir_ch = point_group[self.pg]  # 点群
        self.pos = getc.pos
        self.translate = getc.trans
        self.rot = getc.rot
        self.rot_o3 = getc.rot_o3
        self.order = 4

    def _gen_site_amps(self):
        """
        各サイト座標ごとに1成分だけ1となるアンプリチュード辞書を生成する。
        Returns
        -------
        results : list of dict
            各サイトごとのアンプリチュード辞書リスト
        """
        results = []
        for coord in self.pos:
            amp_site = {}
            x, y, z = coord
            amp_site[(x, y, z)] = 1
            results.append(amp_site)
        return results

    def _search_op(self, op):
        """
        与えられたO(3)行列opに最も近い点群操作名を検索する。
        Returns
        -------
        op_name : str or None
            一致する操作名（なければNone）
        """
        reps = get_rep()[self.pg]
        for rep_name, rep in reps.items():
            diff = np.ravel(rep) - np.ravel(op)
            if np.linalg.norm(diff) < 1e-6:
                return rep_name
        return None
    
    def _search_op2(self, op):
        """
        与えられたO(3)行列opに最も近い点群操作名を検索する。
        Returns
        -------
        op_name : str or None
            一致する操作名（なければNone）
        """
        reps = get_rep()[self.pg]
        delta_list = []
        rep_names = []
        for rep_name, rep in reps.items():
            diff = np.linalg.norm(np.ravel(rep) - np.ravel(op))
            delta_list.append(diff)
            rep_names.append(rep_name)
        delta = np.array(delta_list)
        idx = np.argmin(delta)
        return rep_names[idx]

    def make_salc_site(self, sym_idx=0):
        """
        指定した対称操作インデックスに対するSALC重みを計算する。
        s軌道的な対称性のみ対応。
        Args
        -------
        sym_idx : int
            既約表現の番号
        Returns
        -------
        weights : dict
            各サイト座標ごとのSALC重み
        """
        # s軌道的な対称性のみ対象
        # より精密には軌道のgrid情報必要
        rot = self.rot
        rot_o3 = self.rot_o3
        trs = self.translate
        nops = self.rot.shape[0]
        site_amps = self._gen_site_amps()
        weights = {}
        for samp in site_amps:
            for xyz in samp:
                x, y, z = xyz
                x = round(x, self.order)
                y = round(y, self.order)
                z = round(z, self.order)
                weights[(x, y, z)] = 0
        for samp in site_amps:
            for n in range(nops):
                op_name = self._search_op2(rot_o3[n])
                op_name = op_name.split()[0]
                applied = apply_for_orb(samp, rot[n], trs[n])
                for xyz, amp in applied.items():
                    x, y, z = xyz
                    pos = np.array([x, y, z])
                    pos_diffs = get_diff(pos.copy(), None)
                    weights = self._get_new_weight(samp, amp, pos_diffs, weights, op_name, sym_idx)
                
        return weights
    
    def _get_new_weight(self, site_amp, amp, pos_diffs, weights, op_name, sym_idx):
        ok = False
        for xyz_ori in site_amp:
            x_, y_, z_ = xyz_ori
            pos_ori = np.array([x_, y_, z_])
            for pos_diff in pos_diffs:
                if np.linalg.norm(pos_ori - pos_diff) < 1e-6:
                    x_tr, y_tr, z_tr = pos_diff
                    x_tr = round(x_tr, self.order)
                    y_tr = round(y_tr, self.order)
                    z_tr = round(z_tr, self.order)
                    weights[(x_tr, y_tr, z_tr)] += amp * self.ir_ch[op_name][sym_idx]
                    ok = True
                    break
            if ok:
                break
        return weights
         

    def reconst_umat(self):
        """
        SALCから波動関数（Umat）を再構成する（未実装）。
        """
        raise NotImplementedError

    def get_ham(self):
        """
        SALC基底でのハミルトニアンを返す（未実装）。
        """
        raise NotImplementedError
