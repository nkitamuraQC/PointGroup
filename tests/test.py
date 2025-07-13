from pointgroup import (
    GetSALC,
    GetIR,
    GetCharacter,
    TBModel,
    tb_model,
)
import numpy as np

def test_overall():
    lat=[[1.0, 0, 0],[0, 1.0, 0],[0, 0, 1.0]]
    orb=[[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]
    site_species = [0, 1]
    my_model=tb_model(3, 3, lat, orb)
    my_model.set_hop(-1., 0, 0, [1, 0, 0])

    tb = TBModel(pythtb_obj=my_model, site_species=site_species)
    tb.gen_pythtb()

    gcclass = GetCharacter(tb)
    ch = gcclass.get_character(kidx=0, orb_idx=0, op_idx=1)
    print(ch)
    ir = GetIR(gcclass)
    h = ir._get_h()
    ch = np.ones((h))
    ir.reduce_ir(ch, ir_idx=0)

    salc = GetSALC(gcclass)
    w = salc.make_salc_site(sym_idx=0)
    print(w)
    return


if __name__ == "__main__":
    test_overall()
