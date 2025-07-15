import re


def output_sym_op_name(pg, sym_op_name):
    # pg_ch_extra.pyの全点群に対応したマッピング辞書
    mapping = {
        "C1": {"E": "E"},
        "Ci": {"E": "E", "i": "i"},
        "C2": {"E": "E", "C2": "C2"},
        "Cs": {"E": "E", "sigmah": "sigmah"},
        "C2h": {"E": "E", "C2": "C2", "i": "i", "sigmah": "sigmah"},
        "D2": {"E": "E", "C2(z)": "C2(z)", "C2(y)": "C2(y)", "C2(x)": "C2(x)"},
        "C2v": {"E": "E", "C2(z)": "C2(z)", "sigmav(xz)": "sigmav(xz)", "sigmav(yz)": "sigmav(yz)"},
        "D2h": {"E": "E", "sigmaxy": "sigmaxy", "sigmaxz": "sigmaxz", "sigmayz": "sigmayz", "i": "i", "C2(z)": "C2(z)", "C2(y)": "C2(y)", "C2(x)": "C2(x)"},
        "C4": {"E": "E", "C4": "2C4", "C2": "C2"},
        "S4": {"E": "E", "S4": "2S4", "C2": "C2"},
        "C4h": {"E": "E", "C4": "2C4", "C2": "C2", "sigmah": "sigmah", "S4": "2S4", "i": "i"},
        "D4": {"E": "E", "C4": "2C4", "C2": "C2", "C2_1": "2C2_1", "C2_2": "2C2_2"},
        "C4v": {"E": "E", "C4": "2C4", "C2": "C2", "sigmav": "2sigmav", "sigmad": "2sigmad"},
        "D2d": {"E": "E", "S4": "2S4", "C2": "C2", "C2(x)": "2C2_1", "sigmad": "2sigmad"},
        "D4h": {"E": "E", "C4": "2C4", "C2(z)": "C2", "C2_1": "2C2_1", "C2_2": "2C2_2", "i": "i", "S4": "2S4", "sigmah": "sigmah", "sigmav": "2sigmav", "sigmad": "2sigmad"},
        "C3": {"E": "E", "C3": "2C3"},
        "S6": {"E": "E", "C3": "2C3", "S6": "2S6", "i": "i"},
        "D3": {"E": "E", "C3": "2C3", "C2": "3C2"},
        "C3v": {"E": "E", "C3": "2C3", "sigmav": "3sigmav"},
        "D3d": {"E": "E", "S6": "2S6", "C3": "2C3", "i": "i", "C2": "3C2", "sigmad": "3sigmad"},
        "C6": {"E": "E", "C6": "2C6", "C3": "2C3", "C2": "C2"},
        "C3h": {"E": "E", "C3": "2C3", "sigmah": "sigmah", "S3": "2S3"},
        "C6h": {"E": "E", "C6": "2C6", "C3": "2C3", "C2": "C2", "sigmah": "sigmah", "S6": "2S6", "S3": "2S3", "i": "i"},
        "D6": {"E": "E", "C6": "2C6", "C3": "2C3", "C2": "C2", "C2_1": "3C2_1", "C2_2": "3C2_2"},
        "C6v": {"E": "E", "C6": "2C6", "C3": "2C3", "C2": "C2", "sigmav": "3sigmav", "sigmad": "3sigmad"},
        "D3h": {"E": "E", "C3": "2C3", "C2": "3C2", "sigmah": "sigmah", "S3": "2S3", "sigmav": "3sigmav"},
        "D6h": {"E": "E", "C6": "2C6", "C3": "2C3", "C2": "C2", "C2_1": "3C2_1", "C2_2": "3C2_2", "sigmah": "sigmah", "sigmav": "3sigmav", "sigmad": "3sigmad", "S6": "2S6", "S3": "2S3", "i": "i"},
        "T": {"E": "E", "C3": "8C3", "C2": "3C2"},
        "Th": {"E": "E", "C3": "8C3", "C2": "3C2", "sigmah": "3sigmah", "S6": "8S6", "i": "i"},
        "O": {"E": "E", "C3": "8C3", "C2_2": "6C2_2", "C4": "6C4", "C2_1": "6C2_2"},
        "Td": {"E": "E", "C3": "8C3", "sigmad": "6sigmad", "S4": "6S4", "C2": "3C2"},
        "Oh": {"E": "E", "C3": "8C3", "C2_2": "6C2", "C4": "6C4", "C2": "3C2", "i": "i", "S4": "6S4", "S6": "8S6", "sigmah": "3sigmah", "sigmad": "6sigmad", "sigmav(yz)": "3sigmah", "sigmav(xz)": "3sigmah"},
    }
    ret = ""
    if pg in mapping:
        for key, val in mapping[pg].items():
            if key in sym_op_name:
                if key == "i" and "sigma" in sym_op_name:
                    continue
                ret = val
                break
    return ret

        