import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def make_gif(file, outname):
    # ==== データ読み込み ====
    # 例: CSV (各行が時間ステップ, 各列が空間点)
    file_name = "out/"+file+".csv"
    data = np.loadtxt(file_name, delimiter=",")  # shape = (T, N)

    # 時間の情報が別ファイルや配列にある場合
    # Hs = np.loadtxt("out/Hs.csv", delimiter=",")  # shape = (T, 1) と仮定
    Hs = np.arange(data.shape[0]) 
    # 空間座標 (横軸)
    x = np.linspace(0, 16000, data.shape[1])
    x = x[::-1]

    # ==== 描画の準備 ====
    fig, ax = plt.subplots(figsize=(11, 4))
    line, = ax.plot([], [], marker="o", color="red", linewidth=2, label="Q")

    ax.set_xlim(x.min(), x.max())
    # ax.set_ylim(0, data.max() * 1.1)  # 固定にすると速い
    ax.set_ylim(0,3000)
    ax.set_xlabel("L [m]")
    ax.set_ylabel(r"Q [m$^3$/s]")
    ax.legend(loc="best")

    # ==== update関数 ====
    def update(frame):
        Q = data[frame, :]
        line.set_data(x, Q)
        ax.set_title(f"t = {Hs[frame]:.1f} h")
        return line,

    # ==== アニメーション ====
    ani = animation.FuncAnimation(
        fig, update, frames=data.shape[0], blit=True
    )

    # 保存（fpsを上げると速く見える）[]
    output_name= "out/"+outname+".gif"
    ani.save(output_name, writer="pillow", fps=20) 
    plt.close()

make_gif("Qs","discharge")
