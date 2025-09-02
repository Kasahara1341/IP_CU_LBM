import numpy as np
import matplotlib.pyplot as plt
import csv

# パラメータ設定
a = 0.5   # 短軸 (x)
b = 1.0   # 長軸 (z)
# b = 1.0   # 
n_points = int(100*1.19)  # 点群の数

# 角度を細かく分割して近似的に弧長を計算
n_samples = 10000 
theta_samples = np.linspace(0, 2*np.pi, n_samples)
dx = -a * np.sin(theta_samples)
dy =  b * np.cos(theta_samples)
ds = np.sqrt(dx**2 + dy**2)  # 微小弧長
arc_length = np.cumsum(ds) * (2*np.pi / n_samples)  # 弧長累積
arc_length = arc_length - arc_length[0]
total_length = arc_length[-1]  # 楕円の全周長

# 等間隔の弧長位置を決定
target_s = np.linspace(0, total_length, n_points, endpoint=False)

# 弧長に対応するthetaを補間
theta_at_s = np.interp(target_s, arc_length, theta_samples)

# 座標を計算
x = a * np.cos(theta_at_s)
z = b * np.sin(theta_at_s)
y = np.zeros_like(x)

# CSVに出力
with open("ellipse_points_arclength.csv", "w", newline="") as f:
    writer = csv.writer(f)
    for xi, yi, zi in zip(x, y, z):
        writer.writerow([xi, yi, zi])

# Matplotlibで描画
plt.figure()
plt.scatter(x, z, c="red", marker="o", label="Arc-length points")
plt.gca().set_aspect("equal", adjustable="box")
plt.title("Ellipse point cloud (arc-length spacing)")
plt.xlabel("x")
plt.ylabel("y")
# plt.ylim(-0.05,0.05)
plt.grid(True)
plt.legend()
plt.savefig("./aaa.png")
