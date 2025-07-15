from pointgroup import (
    # GetSALC,
    GetIR,
    GetCharacter,
    TBModel,
    tb_model,
)
import numpy as np
import spglib

def test_overall():
    # Buckled layer
    # real space is 3D
    # define lattice vectors
    lat=[[1.0,0.0,0.0],[0.0,1.25,0.0],[0.0,0.0,3.0]]
    # define coordinates of orbitals
    orb=[[0.0,0.0,-0.15],[0.5,0.5,0.15]]
    site_species = [1, 2]
    
    # only first two lattice vectors repeat, so k-space is 2D
    my_model=tb_model(3,3,lat,orb)
    
    # set model parameters
    delta=1.1
    t=0.6
    
    # set on-site energies
    my_model.set_onsite([-delta,delta])
    # set hoppings (one for each connected pair of orbitals)
    # (amplitude, i, j, [lattice vector to cell containing j])
    my_model.set_hop(t, 1, 0, [0, 0, 0])
    my_model.set_hop(t, 1, 0, [1, 0, 0])
    my_model.set_hop(t, 1, 0, [0, 1, 0])
    my_model.set_hop(t, 1, 0, [1, 1, 0])

    tb = TBModel(pythtb_obj=my_model, site_species=site_species)
    tb.gen_pythtb()

    gcclass = GetCharacter(tb)
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
    nir = ir._get_ir()
    ch = np.array(res)
    w_chs = []
    for i in range(nir):
        w_ch = ir.reduce_ir(ch, ir_idx=i)
        w_chs.append(w_ch)
    print(w_chs)


    # salc = GetSALC(gcclass)
    # w = salc.make_salc_site(sym_idx=0)
    # print(w)
    return

def test_overall2():
    # graphene
    # define lattice vectors
    lat=[[1.0,0.0,0.0],[0.5,np.sqrt(3.0)/2.0,0.0],[0.0,0.0,1.0]]
    # define coordinates of orbitals
    orb=[[1./3.,1./3.,0.0],[2./3.,2./3.0,0.0]]
    site_species = [1, 1]
    cell = (lat, orb, site_species)
    ret = spglib.spglib.standardize_cell(cell)
    lat = ret[0]
    orb = ret[1]
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
    print("kpts: ", tb.kpts[15])

    gcclass = GetCharacter(tb)
    gcclass.use_trace = False
    n = gcclass._get_symm_ops()
    res1 = []
    for i in range(n):
        ch = gcclass.get_character(kidx=15, orb_idx=0, op_idx=i)
        res1.append(ch)
    
    res2 = []
    for i in range(n):
        ch = gcclass.get_character(kidx=15, orb_idx=1, op_idx=i)
        res2.append(ch)
    print(res1)
    print(res1.count(1))
    print(res1.count(-1))
    print(res1.count(2))
    print()
    print(res2)
    print(res2.count(1))
    print(res2.count(-1))
    print(res2.count(2))
    # ir = GetIR(gcclass)
    # h = ir._get_h()
    # nir = ir._get_ir()
    # ch = np.array(res1)
    # w_chs = []
    # for i in range(nir):
    #     w_ch = ir.reduce_ir(ch, ir_idx=i)
    #     w_chs.append(w_ch)
    # print(w_chs)

    # salc = GetSALC(gcclass)
    # for i in range(12):
    #     w = salc.make_salc_site(sym_idx=i)
    #     print(w)
    return

def test_overall3():
    # Checkerboard
    # define lattice vectors
    lat=[[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,10.0]]
    # define coordinates of orbitals
    orb=[[0.0,0.0,0.0],[0.5,0.5,0.0]]
    
    # make two dimensional tight-binding checkerboard model
    my_model=tb_model(3,3,lat,orb)
    
    # set model parameters
    delta=1.1
    t=0.6
    
    # set on-site energies
    my_model.set_onsite([delta,delta])
    # set hoppings (one for each connected pair of orbitals)
    # (amplitude, i, j, [lattice vector to cell containing j])
    my_model.set_hop(t, 1, 0, [0, 0, 0])
    my_model.set_hop(t, 1, 0, [1, 0, 0])
    my_model.set_hop(t, 1, 0, [0, 1, 0])
    my_model.set_hop(t, 1, 0, [1, 1, 0])

    site_species = [1, 2]

    tb = TBModel(pythtb_obj=my_model, site_species=site_species)
    tb.gen_pythtb()

    gcclass = GetCharacter(tb)
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
    nir = ir._get_ir()
    ch = np.array(res)
    w_chs = []
    for i in range(nir):
        w_ch = ir.reduce_ir(ch, ir_idx=i)
        w_chs.append(w_ch)
    print(w_chs)


    # salc = GetSALC(gcclass)
    # w = salc.make_salc_site(sym_idx=0)
    # print(w)
    return

def test_overall4():
    # Oh
    # real space is 3D
    # define lattice vectors
    lat=[[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]]
    # define coordinates of orbitals
    orb=[[0.0,0.0,0.0]]
    site_species = [1]
    
    # only first two lattice vectors repeat, so k-space is 2D
    my_model=tb_model(3,3,lat,orb)
    
    # set model parameters
    delta=1.1
    t=0.6
    
    # set on-site energies
    my_model.set_onsite([0])
    # set hoppings (one for each connected pair of orbitals)
    # (amplitude, i, j, [lattice vector to cell containing j])
    my_model.set_hop(t, 0, 0, [1, 0, 0])
    my_model.set_hop(t, 0, 0, [0, 1, 0])
    my_model.set_hop(t, 0, 0, [0, 0, 1])

    tb = TBModel(pythtb_obj=my_model, site_species=site_species)
    tb.gen_pythtb()

    gcclass = GetCharacter(tb)
    n = gcclass._get_symm_ops()
    res = []
    for i in range(n):
        ch = gcclass.get_character(kidx=0, orb_idx=0, op_idx=i)
        res.append(ch)
    print(res.count(1))
    print(res.count(-1))
    ir = GetIR(gcclass)
    h = ir._get_h()
    nir = ir._get_ir()
    ch = np.array(res)
    w_chs = []
    for i in range(nir):
        w_ch = ir.reduce_ir(ch, ir_idx=i)
        w_chs.append(w_ch)
    print(w_chs)


    # salc = GetSALC(gcclass)
    # w = salc.make_salc_site(sym_idx=0)
    # print(w)
    return


def test_overall5():
    site_species = [1, 2, 3]
    # define lattice vectors
    lat=[[1.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]]
    # define coordinates of orbitals
    orb=[[0.0, 0.0, 0.0],[1.0/3.0, 0.0, 0.0],[2.0/3.0, 0.0, 0.0]]
    
    # make one dimensional tight-binding model
    my_model=tb_model(3,3,lat,orb)
    
    # set model parameters
    delta=2.0
    t=-1.0

    # set hoppings (one for each connected pair of orbitals)
    # (amplitude, i, j, [lattice vector to cell containing j])
    my_model.set_hop(t, 0, 1, [0, 0, 0])
    my_model.set_hop(t, 1, 2, [0, 0, 0])
    my_model.set_hop(t, 2, 0, [1, 0, 0])

    tb = TBModel(pythtb_obj=my_model, site_species=site_species)
    tb.gen_pythtb()

    gcclass = GetCharacter(tb)
    n = gcclass._get_symm_ops()
    res = []
    for i in range(n):
        ch = gcclass.get_character(kidx=3, orb_idx=2, op_idx=i)
        res.append(ch)
    print(res.count(1))
    print(res.count(-1))
    ir = GetIR(gcclass)
    h = ir._get_h()
    nir = ir._get_ir()
    ch = np.array(res)
    w_chs = []
    for i in range(nir):
        w_ch = ir.reduce_ir(ch, ir_idx=i)
        w_chs.append(w_ch)
    print(w_chs)
    return



if __name__ == "__main__":
    # test_overall()
    test_overall2()
    # test_overall3()
    # test_overall4()
    # test_overall5()
