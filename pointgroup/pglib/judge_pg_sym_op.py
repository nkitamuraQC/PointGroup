import numpy as np
from itertools import permutations, product
import math

def generate_all_representation_matrices(n=3):
    """n次元の対称操作表現行列を全て生成"""
    matrices = []
    for perm in permutations(range(n)):
        for signs in product([1, -1], repeat=n):
            matrix = np.zeros((n, n), dtype=int)
            for i in range(n):
                matrix[i, perm[i]] = signs[i]
            matrices.append(matrix)
    return matrices

def get_rotation_axis_and_angle(matrix):
    """回転軸と回転角を求める"""
    # 固有値から回転角を求める
    eigenvals = np.linalg.eigvals(matrix.astype(float))
    
    # 実固有値が1の場合、その固有ベクトルが回転軸
    eigenvals, eigenvecs = np.linalg.eig(matrix.astype(float))
    
    # 固有値が1に近いものを探す
    axis = None
    for i, val in enumerate(eigenvals):
        if abs(val - 1) < 1e-10:  # 実固有値が1
            axis = eigenvecs[:, i].real
            break
    
    if axis is not None:
        # 軸を正規化
        axis = axis / np.linalg.norm(axis)
        
        # 回転角を計算
        trace = np.trace(matrix)
        angle = math.acos((trace - 1) / 2) * 180 / math.pi
        
        return axis, angle
    
    return None, None

def get_mirror_plane_normal(matrix):
    """鏡映面の法線ベクトルを求める"""
    # 固有値が-1の固有ベクトルが鏡映面の法線
    eigenvals, eigenvecs = np.linalg.eig(matrix.astype(float))
    
    for i, val in enumerate(eigenvals):
        if abs(val + 1) < 1e-10:  # 固有値が-1
            normal = eigenvecs[:, i].real
            return normal / np.linalg.norm(normal)
    
    return None

def classify_symmetry_operation(matrix, principal_axis='z'):
    """対称操作を分類する
    
    Args:
        matrix: 対称操作行列
        principal_axis: 主回転軸 ('x', 'y', 'z')
    """
    det = np.linalg.det(matrix)
    trace = np.trace(matrix)
    
    # 恒等操作
    if np.allclose(matrix, np.eye(3)):
        return "E", "Identity"
    
    # 反転中心
    if np.allclose(matrix, -np.eye(3)):
        return "i", "Inversion center"
    
    # 行列式が正の場合（回転操作）
    if det > 0:
        axis, angle = get_rotation_axis_and_angle(matrix)
        
        if axis is not None:
            # 主軸方向を判定
            axis_abs = np.abs(axis)
            main_axis_idx = np.argmax(axis_abs)
            prime = 0
            if abs(axis[main_axis_idx] - 1) <= 1e-6 and main_axis_idx != principal_axis:
                prime = 1
            elif main_axis_idx != principal_axis:
                prime = 2
            axis_names = ['x', 'y', 'z']
            axis_name = axis_names[main_axis_idx]
            
            # 回転角から操作を判定
            if abs(angle - 180) < 1e-6:
                # 主回転軸の場合は区別
                if axis_name == principal_axis:
                    return f"C2({axis_name})", f"180° rotation around principal {axis_name}-axis; {axis}"
                else:
                    if prime == 1:
                        return f"C2_1({axis_name})", f"180° rotation around {axis_name}-axis (perpendicular to principal axis) {axis}"
                    elif prime == 2:
                        return f"C2_2({axis_name})", f"180° rotation around {axis_name}-axis (perpendicular to principal axis) {axis}"
            elif abs(angle - 120) < 1e-6:
                return f"C3({axis_name})", f"120° rotation around {axis_name}-axis"
            elif abs(angle - 90) < 1e-6:
                return f"C4({axis_name})", f"90° rotation around {axis_name}-axis"
            elif abs(angle - 60) < 1e-6:
                return f"C6({axis_name})", f"60° rotation around {axis_name}-axis"
            else:
                if trace == -1:
                    return f"C2({axis_name})", f"{angle:.1f}° rotation around {axis_name}-axis"
                if trace == 0:
                    return f"C3({axis_name})", f"{angle:.1f}° rotation around {axis_name}-axis"
                if trace == 1:
                    return f"C4({axis_name})", f"{angle:.1f}° rotation around {axis_name}-axis"
                if trace == 2:
                    return f"C6({axis_name})", f"{angle:.1f}° rotation around {axis_name}-axis"
    
    # 行列式が負の場合（反転を含む操作）
    else:
        # 鏡映操作の判定
        if trace == 1:  # σ操作
            normal = get_mirror_plane_normal(matrix)
            if normal is not None:
                normal_abs = np.abs(normal)
                main_idx = np.argmax(normal_abs)
                if abs(normal_abs[main_idx]) < 0.9:
                    main_idx = None
                axis_names = ['x', 'y', 'z']
                
                # 主回転軸との関係で分類
                if main_idx == 0:  # x方向が法線 → yz面
                    if principal_axis == 'x':
                        return "sigmah(yz)", "Horizontal mirror plane (yz-plane)"
                    else:
                        return "sigmav(yz)", "Vertical mirror plane (yz-plane)"
                elif main_idx == 1:  # y方向が法線 → xz面
                    if principal_axis == 'y':
                        return "sigmah(xz)", "Horizontal mirror plane (xz-plane)"
                    else:
                        return "sigmav(xz)", "Vertical mirror plane (xz-plane)"
                elif main_idx == 2:  # z方向が法線 → xy面
                    if principal_axis == 'z':
                        return "sigmah(xy)", "Horizontal mirror plane (xy-plane)"
                    else:
                        return "sigmav(xy)", "Vertical mirror plane (xy-plane)"
                else:
                    # 対角鏡映面の場合
                    if abs(normal[0]) > 0.5 and abs(normal[1]) > 0.5:
                        return "sigmad", "Dihedral mirror plane"
                    else:
                        return "sigmad", "Dihedral mirror plane"
                    # else:
                    #     return "sigmav", "Vertical mirror plane"
        
        # 回転反映操作の判定
        elif trace == -1:  # S4操作
            # 回転反映軸を求める
            # axis, _ = get_rotation_axis_and_angle(-matrix)  # -matrixで回転部分を取得
            # if axis is not None:
            #     axis_abs = np.abs(axis)
            #     main_axis_idx = np.argmax(axis_abs)
            #     axis_names = ['x', 'y', 'z']
            #     axis_name = axis_names[main_axis_idx]
            #     return f"S4({axis_name})", f"90° rotation-reflection around {axis_name}-axis"
            return "S4", "90° rotation-reflection"
            
        elif abs(trace) < 1e-6:
            return "S6", "60° rotation-reflection"

        elif abs(trace + 2) < 1e-6:
            return "S3", "120° rotation-reflection"
    
    # return "Unknown", "Unclassified operation"
    raise ValueError(f"不正な表現行列を検知しました: {matrix}")

def analyze_all_operations(principal_axis='z'):
    """全ての対称操作を分析
    
    Args:
        principal_axis: 主回転軸 ('x', 'y', 'z')
    """
    matrices = generate_all_representation_matrices(3)
    
    # 操作を分類
    operations = {}
    for matrix in matrices:
        op_name, description = classify_symmetry_operation(matrix, principal_axis)
        if op_name not in operations:
            operations[op_name] = []
        operations[op_name].append((matrix, description))
    
    # 結果を表示
    print(f"Classification of all symmetry operations (Principal axis: {principal_axis})")
    print("=" * 70)
    
    for op_name, op_list in sorted(operations.items()):
        print(f"\n{op_name}: {len(op_list)} operations")
        print("-" * 40)
        
        for i, (matrix, description) in enumerate(op_list):
            print(f"  {i+1}. {description}")
            print(f"     Matrix:")
            for row in matrix:
                print(f"     {row}")
            # print()
    
    return operations

def find_operation_type(target_matrix, principal_axis='z'):
    """特定の行列の操作タイプを判定"""
    op_name, description = classify_symmetry_operation(target_matrix, principal_axis)
    # print(f"Matrix:")
    # for row in target_matrix:
    #     print(f"  {row}")
    # print(f"Operation: {op_name}")
    # print(f"Description: {description}")
    return op_name, description

def compare_principal_axes():
    """異なる主軸設定での比較"""
    matrices = generate_all_representation_matrices(3)
    
    print("Comparison of mirror plane classifications:")
    print("=" * 60)
    
    # 鏡映操作のみを抽出
    mirror_matrices = []
    for matrix in matrices:
        if np.linalg.det(matrix) < 0 and np.trace(matrix) == 1:
            mirror_matrices.append(matrix)
    
    for i, matrix in enumerate(mirror_matrices):
        print(f"\nMirror operation {i+1}:")
        print("Matrix:")
        for row in matrix:
            print(f"  {row}")
        
        for axis in ['x', 'y', 'z']:
            op_name, description = classify_symmetry_operation(matrix, axis)
            print(f"  Principal axis {axis}: {op_name}")
        print("-" * 30)


if __name__ == "__main__":
    # 主軸をz軸として分析
    print("Analysis with z-axis as principal axis:")
    operations_z = analyze_all_operations('z')
    
    # 鏡映面の分類比較
    print("\n" + "="*70)
    compare_principal_axes()
    
    # 統計情報
    print("\nSummary (Principal axis: z):")
    print("=" * 30)
    for op_name, op_list in sorted(operations_z.items()):
        print(f"{op_name}: {len(op_list)} operations")
    
    # 特定の行列の例
    print("\nExample classifications (Principal axis: z):")
    print("=" * 40)
    
    examples = [
        np.array([[-1, 0, 0], [0, 0, 1], [0, 1, 0]]),
        np.array([[-1, 0, 0], [0, 0, -1], [0, -1, 0]]),
        np.array([[0, -1, 0], [-1, 0, 0], [0, 0, -1]]),
        np.array([[0, 1, 0], [1, 0, 0], [0, 0, -1]]),
        np.array([[0, 0, 1], [0, -1, 0], [1, 0, 0]]),
        np.array([[1, 0, 0], [0, -1, 0], [0, 0, -1]]),  # σ(xz)
        np.array([[-1, 0, 0], [0, 1, 0], [0, 0, -1]])   # σ(yz)
    ]
    
    for i, matrix in enumerate(examples):
        print(f"\nExample {i+1}:")
        find_operation_type(matrix, 'z')