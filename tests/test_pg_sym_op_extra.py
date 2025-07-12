import importlib
from pointgroup.pg_lib import get_rep

def test_import():
    importlib.import_module('pointgroup.pg_lib.pg_sym_op_extra')

def test_num():
    reps = get_rep()
    assert len(reps) == 32
