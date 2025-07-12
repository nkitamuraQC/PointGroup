import importlib
from pointgroup.pglib import point_group

def test_import():
    importlib.import_module('pointgroup.pglib.pg_ch_extra')

def test_num():
    assert len(point_group) == 32



