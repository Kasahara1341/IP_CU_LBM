import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def update2(Qss,line_Q,x,ax,Hs,frame):
    Q = Qss[frame,:]
    line_Q.set_data(x, Q)
    ymin = Q.min() - 0.5
    ymax = Q.max() + 0.5
    ax.set_ylim(0, 2700)
    ax.set_title(f" t = {Hs[frame,0]:.1f} h")
    return line_Q,



def plot_results(Hs, Qs, x):

    Hss = Hs[:,1:]
    fig, ax = plt.subplots(figsize=(11,4))  # 縦横同じくらいに
    line_zb, = ax.plot([], [], color="saddlebrown", linewidth=2, label="bed elevation")
    line_H,  = ax.plot([], [], marker="o", color="blue", linewidth=2, label="water level")
    ax.set_xlim(x.min(), x.max())
    ax.set_xlabel(" L [m]")
    ax.set_ylabel("H [m]")
    ax.legend(loc="best")
    # 縦軸を毎フレーム更新して見やすくする
    def update(ax,line_zb,line_H,Hs,Hss,x,zb,frame):
        h = Hss[frame,:]
        H =  h
        line_zb.set_data(x, zb)
        line_H.set_data(x, H)
        ymin = min(zb.min(), H.min()) - 0.5
        ymax = H.max() + 0.5
        ax.set_ylim(-11, 15)
        ax.set_title(f" t = {Hs[frame,0]:.1f} h")
        return line_zb, line_H

    ani = animation.FuncAnimation(fig, update, frames=Hss.shape[0], blit=True)

    ani.save("./out/waterlevel.gif", writer="pillow", fps=20)
    plt.close()

    Qss = Qs[:,1:]   # 流量データ（時系列 × 空間）

    fig, ax = plt.subplots(figsize=(11,4))

    line_Q, = ax.plot([], [], marker="o", color="red", linewidth=2, label=" Q")

    ax.set_xlim(x.min(), x.max())
    ax.set_xlabel(" L [m]")
    ax.set_ylabel(r"Q [m^3/s]")   # 単位は適宜
    ax.legend(loc="best")

    ani = animation.FuncAnimation(fig, update, frames=Qss.shape[0], blit=True)

    ani.save("./out/discharge.gif", writer="pillow", fps=20)
    plt.close()

    np.savetxt("./out/Qs.csv",Qs, fmt="%.6f", delimiter=",")
    