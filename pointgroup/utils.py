from ase.io import write
from ase.data import chemical_symbols
from ase import Atoms
import logging


def visualize(cell, pos, species, name="vis.cif"):
    symbols = [chemical_symbols[Z] if isinstance(Z, int) else Z for Z in species]
    atoms = Atoms(symbols=symbols, positions=positions, cell=cell, pbc=True)
    write(name, atoms, format="cif")
    return


def get_logger(
    name: str = "app", log_file: str = None, level: int = logging.INFO
) -> logging.Logger:
    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.propagate = False  # 同じメッセージが複数回出力されるのを防ぐ

    if not logger.handlers:  # ハンドラがすでに設定されていなければ
        # コンソール出力
        console_handler = logging.StreamHandler()
        console_handler.setLevel(level)

        formatter = logging.Formatter(
            fmt="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )
        console_handler.setFormatter(formatter)
        logger.addHandler(console_handler)

        # ファイル出力が必要なら
        if log_file:
            file_handler = logging.FileHandler(log_file, encoding="utf-8")
            file_handler.setLevel(level)
            file_handler.setFormatter(formatter)
            logger.addHandler(file_handler)

    return logger


logger = get_logger()
