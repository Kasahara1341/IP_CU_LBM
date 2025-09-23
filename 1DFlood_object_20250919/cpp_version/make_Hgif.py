import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# ==== データ読み込み ====
# 例: CSV (各行が時間ステップ, 各列が空間点)
data1 = np.loadtxt("out/Hs.csv", delimiter=",")  # shape = (T, N)
data2 = np.loadtxt("out/zb.csv", delimiter=",")  # shape = (T, N)

# 時間の情報が別ファイルや配列にある場合
# Hs = np.loadtxt("out/Hs.csv", delimiter=",")  # shape = (T, 1) と仮定
Hs = np.arange(data1.shape[0]) 
zb = np.arange(data2.shape[0])
# 空間座標 (横軸)
x = np.linspace(0, 16000, data1.shape[1])
x = x[::-1]

# ==== 描画の準備 ====
fig, ax = plt.subplots(figsize=(11, 4))
line_zb, = ax.plot([], [], color="saddlebrown", linewidth=2, label="bed elevation")
line_H,  = ax.plot([], [], marker="o", color="blue", linewidth=2, label="water level")

ax.set_xlim(x.min(), x.max())
ax.set_ylim(-12,15)
ax.set_xlabel("L [m]")
ax.set_ylabel(r"H [m]")
ax.legend(loc="best")

# ==== update関数 ====
def update(frame):
    H = data1[frame, :]
    Z = data2[frame, :]
    line_H.set_data(x, H)
    line_zb.set_data(x, Z)
    ymin = min(zb.min(),H.min()) - 0.5
    ymax = H.max() + 0.5
    ax.set_title(f"t = {Hs[frame]:.1f} h")
    return line_H, line_zb

# ==== アニメーション ====
ani = animation.FuncAnimation(
    fig, update, frames=data1.shape[0], blit=True
)

# 保存（fpsを上げると速く見える）
ani.save("out/H_zb.gif", writer="pillow", fps=20) 
plt.close()
