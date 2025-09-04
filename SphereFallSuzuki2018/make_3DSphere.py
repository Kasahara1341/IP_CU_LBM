import numpy as np
import matplotlib.pyplot as plt

# 黄金比
phi = (1 + np.sqrt(5)) / 2

# 正20面体の頂点
vertices = np.array([
    [0,  1,  phi], [0, -1,  phi], [0,  1, -phi], [0, -1, -phi],
    [ 1,  phi, 0], [-1,  phi, 0], [ 1, -phi, 0], [-1, -phi, 0],
    [ phi, 0,  1], [-phi, 0,  1], [ phi, 0, -1], [-phi, 0, -1]
], dtype=float)

# 正規化して半径1の球面に乗せる
vertices /= np.linalg.norm(vertices[0])

# 正20面体の20枚の三角形面（頂点インデックス）
faces = np.array([
    [0, 1, 8], [0, 1, 9], [0, 4, 5], [0, 4, 8], [0, 5, 9],
    [1, 6, 7], [1, 6, 8], [1, 7, 9], [2, 3,10], [2, 3,11],
    [2, 4, 5], [2, 4,10], [2, 5,11], [3, 6, 7], [3, 6,10],
    [3, 7,11], [4, 8,10], [5, 9,11], [6, 8,10], [7, 9,11]
])

def subdivide_icosahedron(n):
    """
    正20面体をn分割して球面上の点群を返す
    """
    points = []
    for tri in faces:
        A, B, C = vertices[tri]
        # 三角形ABCをn分割
        for i in range(n+1):
            for j in range(n+1-i):
                k = n - i - j
                # 重心座標 (i,j,k)/n
                P = (i*A + j*B + k*C) / n
                P /= np.linalg.norm(P)  # 球面に射影
                points.append(tuple(P))
    # 重複を除去
    points = np.unique(np.array(points), axis=0)
    return points

# ===== 実行例 =====
if __name__ == "__main__":
    n = int(input("input number of division : "))  # 分割数を変えて試せる (1,2,3,...)
    pts = subdivide_icosahedron(n)

    print(f"n={n}, 点の数={len(pts)}")
    # === CSVに保存 ===
    filename = f"sphere_points_n{n}.csv"
    np.savetxt(filename, pts, delimiter=",", comments="")
    print(f"点群データを {filename} に保存しました")    