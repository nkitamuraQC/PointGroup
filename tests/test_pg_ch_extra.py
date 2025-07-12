import importlib
from pointgroup.pg_lib import point_group

def test_import():
    importlib.import_module('pointgroup.pg_lib.pg_ch_extra')

def test_num():
    assert len(point_group) == 32

    

