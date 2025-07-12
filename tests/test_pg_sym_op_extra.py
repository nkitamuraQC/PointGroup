import importlib
from pointgroup.pglib import get_rep

def test_import():
    importlib.import_module('pointgroup.pglib.pg_sym_op_extra')

def test_num():
    reps = get_rep()
    assert len(reps) == 32
