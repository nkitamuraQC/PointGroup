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
    site_species = [1, 2]
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

def test_overall2():
    # define lattice vectors
    lat=[[1.0,0.0,0.0],[0.5,np.sqrt(3.0)/2.0,0.0],[0.0,0.0,10.0]]
    # define coordinates of orbitals
    orb=[[1./3.,1./3.,0.0],[2./3.,2./3.0,0.0]]
    site_species = [1, 1]
    my_model=tb_model(3, 3, lat, orb)
    delta=0.0
    t=-1.0
    
    # set on-site energies
    my_model.set_onsite([-delta,delta])
    # set hoppings (one for each connected pair of orbitals)
    # (amplitude, i, j, [lattice vector to cell containing j])
    my_model.set_hop(t, 0, 1, [ 0, 0, 0])
    my_model.set_hop(t, 1, 0, [ 1, 0, 0])
    my_model.set_hop(t, 1, 0, [ 0, 1, 0])

    nk = [4, 4, 1]

    tb = TBModel(pythtb_obj=my_model, site_species=site_species, nk=nk)
    tb.gen_pythtb()

    gcclass = GetCharacter(tb)
    n = gcclass._get_symm_ops()
    res = []
    for i in range(n):
        ch = gcclass.get_character(kidx=0, orb_idx=1, op_idx=i)
        res.append(ch)
    print(res.count(1))
    print(res.count(-1))
    ir = GetIR(gcclass)
    h = ir._get_h()
    ch = np.ones((h))
    w_ch = ir.reduce_ir(ch, ir_idx=0)
    print(w_ch)

    salc = GetSALC(gcclass)
    for i in range(12):
        w = salc.make_salc_site(sym_idx=i)
        print(w)
    return


if __name__ == "__main__":
    test_overall2()
