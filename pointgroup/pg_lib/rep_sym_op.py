from scipy.spatial.transform import Rotation as R
import numpy as np

# 回転軸: z軸, 回転角: 90度（π/2ラジアン）
axis = np.array([0, 0, 1])
angle_rad = np.pi / 2
rotvec = axis * angle_rad  # 回転ベクトル（ロドリゲスベクトル）

# Rotation オブジェクトを生成
rotation = R.from_rotvec(rotvec)

# 回転行列を取得
rotation_matrix = rotation.as_matrix()

print("回転行列:")
print(rotation_matrix)
