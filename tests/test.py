from pointgroup import (
    GetSALC,
    GetIR,
    GetCharacter,
    TBModel,
    tb_model,
)

def test_overall():
    lat=[[1.0, 0, 0],[0, 1.0, 0],[0, 0, 1.0]]
    orb=[[0.0]]
    my_model=tb_model(3, 3, lat, orb)
    my_model.set_hop(-1., 0, 0, [1, 0, 0])

    tb = TBModel(pythtb_obj=my_model)
    tb.gen_pythtb()

    gcclass = GetCharacter(tb)
    gcclass.get_character(0, 0, 0)
    ir = GetIR(gcclass)
    ir.reduce_ir(0)

    salc = GetSALC(gcclass)
    salc.make_salc_site(0)
    return


if __name__ == "__main__":
    test_overall()
