### 主要な32点群のうち、代表的な点群ごとに対称操作の行列表現を辞書で整理します。
# 各点群に対応する操作は簡潔な3x3行列（直交変換）として与えます。

from .rep_sym_op import *

# 代表的な点群の操作（例: C1, Ci, Cs, C2, C2v, D2, D2h, C3, C3v, D3, D3h, Td, Oh, Ih）
# 一部を例示し、残りは後で補完可能です

def get_rep(cell)
    E = identity()
    point_groups = {}
    
    # C1 (identity only)
    point_groups['C1'] = {'E': E}
    
    # Ci (identity, inversion)
    point_groups['Ci'] = {
        'E': E,
        'i': inversion(cell)
    }
    
    # Cs (mirror plane, e.g. xz)
    point_groups['Cs'] = {
        'E': E,
        'sigma(xz)': reflection_matrix('xz', cell)
    }
    
    # C2 (rotation 180° around z)
    point_groups['C2'] = {
        'E': E,
        'C2(z)': rotation_matrix('z', 180, cell)
    }
    
    # C2v (E, C2, σv(xz), σv'(yz))
    point_groups['C2v'] = {
        'E': E,
        'C2(z)': rotation_matrix('z', 180, cell),
        'sigma(xz)': reflection_matrix('xz', cell),
        'sigma(yz)': reflection_matrix('yz', cell)
    }
    
    # D2 (E, C2(x), C2(y), C2(z))
    point_groups['D2'] = {
        'E': E,
        'C2(x)': rotation_matrix('x', 180, cell),
        'C2(y)': rotation_matrix('y', 180, cell),
        'C2(z)': rotation_matrix('z', 180, cell)
    }
    
    # C3 (E, C3, C3^2)
    point_groups['C3'] = {
        'E': E,
        'C3(z)': rotation_matrix('z', 120, cell)
        'C3^2(z)': rotation_matrix('z', 240, cell)
    }
    
    # C3v (E, 2C3, 3σv)
    point_groups['C3v'] = {
        'E': E,
        'C3(z)': rotation_matrix('z', 120, cell),
        'C3^2(z)': rotation_matrix('z', 240, cell),
        'sigmav1': reflection_matrix('xz', cell),
        'sigmav2': gen_sigmav2(cell),
        'sigmav3': gen_sigmav3(cell)
    }
    
    # Td (tetrahedral)
    # Approximating with identity and representative operations
    point_groups['Td'] = {
        'E': E,
        'C3': rotation_matrix('z', 120, cell),
        'C2': rotation_matrix('x', 180, cell),
        'S4': rotation_matrix('z', 90, cell),  # improper
        'sigmad': reflection_matrix('xz', cell)
    }
    
    # Oh (octahedral)
    point_groups['Oh'] = {
        'E': E,
        'C4(z)': rotation_matrix('z', 90, cell),
        'C3': rotation_matrix('z', 120, cell),
        'C2(x)': rotation_matrix('x', 180, cell),
        'i': inversion(cell)
    }
    
    # point_groups.keys()
    
    # 残りの32点群の対称操作を追加して完成させる
    
    # 補助関数（再掲）
    
    # Monoclinic point groups
    point_groups['C2h'] = {
        'E': E,
        'C2(z)': rotation_matrix('z', 180, cell),
        'i': inversion(cell),
        'sigmah(xy)': reflection_matrix('xy', cell)
    }
    
    # Orthorhombic
    point_groups['D2h'] = {
        'E': E,
        'C2(x)': rotation_matrix('x', 180, cell),
        'C2(y)': rotation_matrix('y', 180, cell),
        'C2(z)': rotation_matrix('z', 180, cell),
        'i': inversion(cell),
        'sigma(xy)': reflection_matrix('xy', cell),
        'sigma(yz)': reflection_matrix('yz', cell),
        'sigma(xz)': reflection_matrix('xz', cell)
    }
    
    # Tetragonal point groups
    point_groups['C4'] = {
        'E': E,
        'C4(z)': rotation_matrix('z', 90, cell),
        'C4^2(z)': rotation_matrix('z', 180, cell),
        'C4^3(z)': rotation_matrix('z', 270, cell)
    }
    
    point_groups['S4'] = {
        'E': E,
        'S4': improper_rotation_matrix('z', 90, cell),
        'S4^2': rotation_matrix('z', 180, cell),
        'S4^3': improper_rotation_matrix('z', 270, cell)
    }
    
    point_groups['C4h'] = {
        'E': E,
        'C4(z)': rotation_matrix('z', 90, cell),
        'C2(z)': rotation_matrix('z', 180, cell),
        'C4^3(z)': rotation_matrix('z', 270, cell),
        'i': inversion(cell),
        'S4': improper_rotation_matrix('z', 90, cell),
        'S4^3': improper_rotation_matrix('z', 270, cell),
        'sigmah(xy)': reflection_matrix('xy', cell)
    }
    
    point_groups['D4'] = {
        'E': E,
        'C4(z)': rotation_matrix('z', 90),
        'C2(z)': rotation_matrix('z', 180),
        'C4^3(z)': rotation_matrix('z', 270),
        'C2(x)': rotation_matrix('x', 180),
        'C2(y)': rotation_matrix('y', 180),
        'C2(d1)': gen_c2_d1(cell),
        'C2(d2)': gen_c2_d2(cell)
    }
    
    point_groups['C4v'] = {
        'E': E,
        'C4(z)': rotation_matrix('z', 90),
        'C2(z)': rotation_matrix('z', 180),
        'C4^3(z)': rotation_matrix('z', 270),
        'sigmav(xz)': reflection_matrix('xz'),
        'sigmav(yz)': reflection_matrix('yz'),
        'sigmad1': gen_sigmad1(cell),
        'sigmad2': gen_sigmad2(cell)
    }
    
    point_groups['D2d'] = {
        'E': E,
        'C2(z)': rotation_matrix('z', 180),
        'C2(x)': rotation_matrix('x', 180),
        'C2(y)': rotation_matrix('y', 180),
        'S4': improper_rotation_matrix('z', 90),
        'S4^3': improper_rotation_matrix('z', 270),
        'sigmad1': gen_sigmad1(cell),
        'sigmad2': gen_sigmad2(cell)
    }
    
    point_groups['D4h'] = {
        'E': E,
        'C4(z)': rotation_matrix('z', 90, cell),
        'C2(z)': rotation_matrix('z', 180, cell),
        'C4^3(z)': rotation_matrix('z', 270, cell),
        'C2(x)': rotation_matrix('x', 180, cell),
        'C2(y)': rotation_matrix('y', 180, cell),
        'i': inversion(),
        'S4': improper_rotation_matrix('z', 90, cell),
        'S4^3': improper_rotation_matrix('z', 270, cell),
        'σh': reflection_matrix('xy', cell),
        'σv(xz)': reflection_matrix('xz', cell),
        'σv(yz)': reflection_matrix('yz', cell)
    }
    
    # Trigonal point groups
    point_groups['S6'] = {
        'E': E,
        'S6': improper_rotation_matrix('z', 60, cell),
        'S6^5': improper_rotation_matrix('z', 300, cell),
        'C3': rotation_matrix('z', 120, cell),
        'C3^2': rotation_matrix('z', 240, cell),
        'i': inversion(cell)
    }
    
    point_groups['D3'] = {
        'E': E,
        'C3': rotation_matrix('z', 120, cell),
        'C3^2': rotation_matrix('z', 240, cell),
        'C2(1)': rotation_matrix('x', 180, cell),
        'C2(2)': gen_c2_2(cell),
        'C2(3)': gen_c2_3(cell)
    }
    
    point_groups['D3d'] = {
        'E': E,
        'C3': rotation_matrix('z', 120),
        'C3^2': rotation_matrix('z', 240),
        'C2': rotation_matrix('x', 180),
        'i': inversion(),
        'S6': improper_rotation_matrix('z', 60),
        'sigmad1': gen_sigmad1_2(cell),
        'sigmad2': gen_sigmad2_2(cell),
        'sigmad3': gen_sigmad3(cell)
    }
    
    # Hexagonal point groups
    point_groups['C6'] = {
        'E': E,
        'C6': rotation_matrix('z', 60, cell),
        'C3': rotation_matrix('z', 120, cell),
        'C2': rotation_matrix('z', 180, cell),
        'C6^5': rotation_matrix('z', 300, cell)
    }
    
    point_groups['C3h'] = {
        'E': E,
        'C3': rotation_matrix('z', 120),
        'C3^2': rotation_matrix('z', 240),
        'sigmah': reflection_matrix('xy'),
        'S3': improper_rotation_matrix('z', 120),
        'S3^2': improper_rotation_matrix('z', 240)
    }
    
    point_groups['C6h'] = {
        'E': E,
        'C6': rotation_matrix('z', 60, cell),
        'C3': rotation_matrix('z', 120, cell),
        'C2': rotation_matrix('z', 180, cell),
        'C6^5': rotation_matrix('z', 300, cell),
        'i': inversion(cell),
        'S3': improper_rotation_matrix('z', 120, cell),
        'sigmah': reflection_matrix('xy', cell)
    }
    
    point_groups['D6'] = {
        'E': E,
        'C6': rotation_matrix('z', 60, cell),
        'C3': rotation_matrix('z', 120, cell),
        'C2': rotation_matrix('z', 180, cell),
        'C2(x)': rotation_matrix('x', 180, cell),
        'C2(y)': rotation_matrix('y', 180, cell)
    }
    
    point_groups['C6v'] = {
        'E': E,
        'C6': rotation_matrix('z', 60, cell),
        'C3': rotation_matrix('z', 120, cell),
        'C2': rotation_matrix('z', 180, cell),
        'sigmav': reflection_matrix('xz', cell),
        'sigmav2': gen_sigmav2_2(cell),
        'sigmav3': gen_sigmav3_2(cell)
    }
    
    point_groups['D3h'] = {
        'E': E,
        'C3': rotation_matrix('z', 120, cell),
        'C3^2': rotation_matrix('z', 240, cell),
        'C2': rotation_matrix('x', 180, cell),
        'sigmah': reflection_matrix('xy', cell),
        'sigmav1': reflection_matrix('xz', cell),
        'sigmav2': gen_sigmav2(cell),
        'sigmav3': gen_sigmav3(cell)
    }
    
    point_groups['D6h'] = {
        'E': E,
        'C6': rotation_matrix('z', 60, cell),
        'C3': rotation_matrix('z', 120, cell),
        'C2': rotation_matrix('z', 180, cell),
        'C2(x)': rotation_matrix('x', 180, cell),
        'i': inversion(),
        'sigmah': reflection_matrix('xy', cell),
        'sigmad1': reflection_matrix('xz', cell),
        'sigmad2': gen_sigmad2_3(cell)
    }
    
    # Cubic (すでに Td, Oh 追加済)
    point_groups['T'] = {
        'E': E,
        'C3': rotation_matrix('z', 120, cell),
        'C2': rotation_matrix('x', 180, cell)
    }
    
    point_groups['Th'] = {
        'E': E,
        'C3': rotation_matrix('z', 120, cell),
        'C2': rotation_matrix('x', 180, cell),
        'i': inversion()
    }
    
    point_groups['O'] = {
        'E': E,
        'C4(z)': rotation_matrix('z', 90, cell),
        'C3(1)': rotation_matrix('z', 120, cell),
        'C3(2)': rotation_matrix('z', 240, cell),
        'C2(x)': rotation_matrix('x', 180, cell),
        'C2(y)': rotation_matrix('y', 180, cell),
        'C2(z)': rotation_matrix('z', 180, cell)
    }


    # 完成した点群の個数を確認
    len(point_groups)  # Should be 32 now
    return point_groups


