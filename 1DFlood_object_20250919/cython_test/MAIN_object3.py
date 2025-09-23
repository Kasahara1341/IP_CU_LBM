from Element import Element
from Node import Node
from Calc_equation import calc_mainloop

import numpy as np
import pandas as pd
import time as timers
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# 以下，メイン文
# 河川の地形の作成
z = np.loadtxt("input/zb_river.txt")
w = np.loadtxt("input/width.txt")
sec = np.loadtxt("input/sec_map.txt")

#本川の地形配列の作成
zb = []
width = []
x = []  # x座標(今は雑)
count = 0
for i in range(sec.shape[0]):
    for j in range(sec.shape[1]):
        if sec[i,j] < 114 and sec[i,j] >= 1:
            zb.append(z[i,j])
            width.append(w[i,j])
            x.append(count*150)  # 150m間隔に格子を配置
            count = count +1
            #print(sec[i,j])
print("number of elements:",len(zb))

# 支川の地形配列の作成
zb2 = []
width2 = []
x2 = []  # x座標(今は雑)
count = 0
xjunc = 97*150   # secmapで言うところの97で合流していたので，このような処理としている

for i in range(sec.shape[0]):
    for j in range(sec.shape[1]):
        if sec[i,j] < 350 and sec[i,j] >= 329:
            zb2.append(z[i,j])
            width2.append(w[i,j])
            x2.append(xjunc + count*150)  # 合流点から見た河口からの距離
            count = count +1



# 地形配列
zb = np.array(zb)
x = np.array(x)
width = np.array(width)

zb2 = np.array(zb2)
x2 = np.array(x2)
width2 = np.array(width2)

# 直線流路化
# zb = np.linspace(zb[0], zb[-1], len(zb))

zb = zb[::-1]
x  = x[::-1]
width = width[::-1]

zb2 = zb2[::-1]
x2  = x2[::-1]
width2 = width2[::-1]


ib = (zb[-1]-zb[0])/(x[-1] - x[0] )

# 流量データの読み込み
init_dt = 2.5
dt = init_dt
df = pd.read_csv("input/Boundary2.csv")
Qb = list(df[df.keys()[1]])
Qb2 = list(df[df.keys()[2]])
Qb = np.array(Qb)
Qb2 = np.array(Qb2)
maxt = Qb.shape[0]
np.savetxt("input/zb.csv",zb, fmt="%.6f", delimiter=",")
np.savetxt("input/zb2.csv",zb2, fmt="%.6f", delimiter=",")
# print(x)
# print(zb)
# print(width)
# print(x2)
# print(zb2)
# print(width2)

# 河川のエレメント＆ノードオブジェクトへの属性の指定
elements = [Element(x[i],zb[i],0.03,width[i]) for i in range(zb.shape[0])]
nodes = [Node() for i in range(zb.shape[0])]
n_riv = len(elements)

elements2 = [Element(x2[i],zb2[i],0.03,width2[i]) for i in range(zb2.shape[0])]
nodes2 = [Node() for i in range(zb2.shape[0])]
n_riv2 = len(elements2)

# エレメント＆ノードの接続(本川)
for i in range(len(elements)):
    nodes[i].set_upelement(elements[i])
    
for i in range(1,len(elements)):
    elements[i].set_upnoad(nodes[i-1])

for i in range(len(elements)-1):
    nodes[i].set_downelement(elements[i+1])
    
for i in range(len(elements)):
    elements[i].set_dnnoad(nodes[i])

# エレメント＆ノードの接続(支川)
for i in range(len(elements2)):
    nodes2[i].set_upelement(elements2[i])
    
for i in range(1,len(elements2)):
    elements2[i].set_upnoad(nodes2[i-1])

for i in range(len(elements2)-1):
    nodes2[i].set_downelement(elements2[i+1])
    
for i in range(len(elements2)):
    elements2[i].set_dnnoad(nodes2[i])


# 境界条件(本川)
bc_upnode = Node()
elements[0].set_upnoad(bc_upnode)
bc_dnelement = Element(x[n_riv-1]-150,zb[n_riv-1]-ib*150,0.03,width[n_riv-1])
nodes[n_riv-1].set_downelement(bc_dnelement)
bc_dnelement.set_upnoad(nodes[n_riv-1])

# 境界条件(支川)
bc_upnode2 = Node()
elements2[0].set_upnoad(bc_upnode2)

#支川と本川の接続 # 96番目のエレメントが接続部
nodes2[-1].set_downelement(elements[16])   
elements[16].set_upnoad(nodes2[-1])

for taget_element in [elements, elements2]:
    for target in taget_element:
        # target.set_Euler()                ; number_of_stage = 1
        # target.set_AdamsBashforth4th(dt)  ; number_of_stage = 1
        target.set_Runge_Kutta_4th()      ; number_of_stage = 4
        # target.set_Runge_Kutta_6th()      ; number_of_stage = 6

# 初期条件
for target in [elements, elements2]:
    for element in target:
        element.set_depth(0)
for taget_node in [nodes, nodes2]:
    for node in taget_node:
        node.set_q(0)

bc_upnode.set_q(Qb[0])
bc_upnode2.set_q(Qb2[0])
bc_dnelement.set_depth(elements[n_riv-1].get_variable_depth())

Hs = []
Qs = []

time = 0
t = 0
maxt = 10
# maxt = 0.00125 *2

### check variables ###
total_Qin   = 0
total_water = 0
total_Qout  = 0
bc_upQin    = [0,0,0,0] ; bc_upQ2in   = [0,0,0,0]
bc_dnQout   = [0,0,0,0]
AB4thcoeff  = [55.0/24.0,-59.0/24.0,37.0/24.0,-9.0/24.0]
q_plot = []  # 流量の縦断分布(プロット用)
t_plot = []  # x座標(プロット用)# 計算開始時刻を記録
start_time = timers.time()
#####################
### 時間発展ループ ###
#####################
while time/3600 < maxt:
    H = []  # 水深の縦断分布
    Q = []  # 流量の縦断分布
    
    H.append(time/3600)
    Q.append(time/3600)
    total_water = 0

    # boundary contion
    bc_upnode.set_q(Qb[int(time // 3600)])
    bc_upnode2.set_q(Qb2[int(time // 3600)])

    """for i in range(3):
        bc_upQin[3-i] =  bc_upQin[2-i]
        bc_upQ2in[3-i] =  bc_upQ2in[2-i]
    bc_upQin[0] = Qb[int(time // 3600)]
    bc_upQ2in[0] = Qb2[int(time // 3600)]
    for i in range(4):
        total_Qin += bc_upQin[i] *AB4thcoeff[i]*dt
        total_Qin += bc_upQ2in[i]*AB4thcoeff[i]*dt
    bc_dnelement.set_depth(elements[n_riv-1].get_variable_depth())
    # Adams Bashforth採用時の本川下端での流出量計算
    for i in range(3):
        bc_dnQout[3-i] =  bc_dnQout[2-i]
        total_Qout     += bc_dnQout[3-i]*AB4thcoeff[3-i]*dt
    bc_dnQout[0]       =  elements[-1].dnnoads[0].get_variable_q()
    total_Qout += bc_dnQout[0]*AB4thcoeff[0]*dt """

    # main calculation
    calc_mainloop([elements, elements2], [nodes, nodes2], dt,H,Q,number_of_stage)

    # dtの判定
    H = np.array(H)
    Q = np.array(Q)


    if Qb[int(time // 3600)]>20000:
        dt = 1
        for taget_element in [elements, elements2]:
            for target in taget_element:
                target.set_Runge_Kutta_6th() ; number_of_stage = 6
    else:
        for taget_element in [elements, elements2]:
            for target in taget_element:
                target.set_Runge_Kutta_4th() ; number_of_stage = 4
        dt = init_dt

    t = t+1
    time = time+dt
    # print("total_Qin:",total_Qin," total_water:",total_water," total_Qout:",total_Qout,"sum/Qin:",(total_Qout+total_water-total_Qin)/total_Qin)
    # print(Qb[int(time // 3600)],Qb2[int(time // 3600)],dt,time/3600)
    if np.mod(time/3600,1)==0:
        # print("total_Qin:",total_Qin," total_water:",total_water," total_Qout:",total_Qout,"sum:",total_Qout+total_water-total_Qin)
        # print("numerical error = ",(total_Qout+total_water-total_Qin)/total_Qin)
        # print(Q[int(len(elements)/2)])
        Hs.append(H[:len(elements)+1])
        Qs.append(Q[:len(elements)+1])
        # q_plot.append(elements[len(elements)//2].dnnoads[0].get_variable_q())
        # t_plot.append(time/3600)
        print("time:",time/3600,"  [h]",r"Q(m^3/s):  dt=",dt)
# """
# 計算終了時刻を記録
end_time = timers.time()
# 経過時間を表示
print(f"計算時間: {end_time - start_time:.3f} 秒")
# save
Hs = np.array(Hs)
Qs = np.array(Qs)

figure = plt.figure()
plt.plot(t_plot,q_plot,label="Q")
plt.xlabel("t [h]")
plt.ylabel("Q [m^3/s]")
plt.grid()
plt.legend()
plt.savefig("./out/hydroGraph.png")

# 水深⇒zbの変換
for i in range(zb.shape[0]):
    Hs[:,i+1] = Hs[:,i+1]+zb[i]

np.savetxt("./out/Hs.csv",Hs, fmt="%.6f", delimiter=",")
np.savetxt("./out/Qs.csv",Qs, fmt="%.6f", delimiter=",")
print("output Hs.csv, Qs.csv")

Hss = Hs[:,1:]

fig, ax = plt.subplots(figsize=(11,4))  # 縦横同じくらいに

line_zb, = ax.plot([], [], color="saddlebrown", linewidth=2, label="bed elevation")
line_H,  = ax.plot([], [], marker="o", color="blue", linewidth=2, label="water level")

ax.set_xlim(x.min(), x.max())
ax.set_xlabel(" L [m]")
ax.set_ylabel("H [m]")
ax.legend(loc="best")

# 縦軸を毎フレーム更新して見やすくする
def update(frame):
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
print("output waterlevel.gif")

Qss = Qs[:,1:]   # 流量データ（時系列 × 空間）

fig, ax = plt.subplots(figsize=(11,4))

line_Q, = ax.plot([], [], marker="o", color="red", linewidth=2, label=" Q")

ax.set_xlim(x.min(), x.max())
ax.set_xlabel(" L [m]")
ax.set_ylabel(r"Q [m^3/s]")   # 単位は適宜
ax.legend(loc="best")

def update(frame):
    Q = Qss[frame,:]
    line_Q.set_data(x, Q)
    ymin = Q.min() - 0.5
    ymax = Q.max() + 0.5
    ax.set_ylim(0, 2700)
    ax.set_title(f" t = {Hs[frame,0]:.1f} h")
    return line_Q,

ani = animation.FuncAnimation(fig, update, frames=Qss.shape[0], blit=True)

ani.save("./out/discharge.gif", writer="pillow", fps=20)
plt.close()

np.savetxt("./out/Qs.csv",Qs, fmt="%.6f", delimiter=",")

# """