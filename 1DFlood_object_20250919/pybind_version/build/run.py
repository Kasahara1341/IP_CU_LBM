import my_module
import numpy as np
import pandas as pd
import time as timers
import matplotlib.pyplot as plt
import matplotlib.animation as animation

z = np.loadtxt("../input/zb_river.txt")
w = np.loadtxt("../input/width.txt")
sec = np.loadtxt("../input/sec_map.txt")

#本川の地形配列の作成
# 条件に合うマスクを作る
mask = (sec >= 1) & (sec < 114)
# 条件に合う要素の数
n_points = np.count_nonzero(mask)
# x座標 (合流点からの距離)
x = np.arange(n_points) * 150
# zb2
zb = z[mask]
# width2 (今は固定値100)
width = w[mask]

# 支川の地形配列の作成
xjunc = 97*150   # secmapで言うところの97で合流していたので，このような処理としている
mask = (sec >= 329) & (sec < 350)
n_points = np.count_nonzero(mask)
x2 = xjunc + np.arange(n_points) * 150
zb2 = z[mask]
width2 = w[mask]

# 直線流路化
#zb = np.linspace(zb[0], zb[-1], len(zb))

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
df = pd.read_csv("../input/Boundary2.csv")
Qb = list(df[df.keys()[1]])
Qb2 = list(df[df.keys()[2]])
print(Qb)
Qb = np.array(Qb)
Qb2 = np.array(Qb2)
maxt = Qb.shape[0]
np.savetxt("../input/zb.csv",zb, fmt="%.6f", delimiter=",")
np.savetxt("../input/zb2.csv",zb2, fmt="%.6f", delimiter=",")


# 河川のエレメント＆ノードオブジェクトへの属性の指定
elements = [my_module.Element(x[i],150,zb[i],0.03,width[i]) for i in range(zb.shape[0])]
nodes = [my_module.Node() for i in range(zb.shape[0])]
n_riv = len(elements)

elements2 = [my_module.Element(x2[i],150,zb2[i],0.03,width2[i]) for i in range(zb2.shape[0])]
nodes2 = [my_module.Node() for i in range(zb2.shape[0])]
n_riv2 = len(elements2)

# エレメント＆ノードの接続(本川)
for i in range(len(elements)):
    nodes[i].set_up_element(elements[i])
    elements[i].set_dn_node(nodes[i])
    
for i in range(1,len(elements)):
    elements[i].set_up_node(nodes[i-1])
    nodes[i-1].set_dn_element(elements[i])

# エレメント＆ノードの接続(支川)
for i in range(len(elements2)):
    nodes2[i].set_up_element(elements2[i])
    elements2[i].set_dn_node(nodes2[i])
    
for i in range(1,len(elements2)):
    elements2[i].set_up_node(nodes2[i-1])
    nodes2[i-1].set_dn_element(elements2[i])


# 境界条件(本川)
bc_upnode = my_module.Node()
elements[0].set_up_node(bc_upnode)
bc_dnelement = my_module.Element(x[n_riv-1]-150,150,zb[n_riv-1]-ib*150,0.03,width[n_riv-1])
nodes[n_riv-1].set_dn_element(bc_dnelement)

# 境界条件(支川)
bc_upnode2 = my_module.Node()
elements2[0].set_up_node(bc_upnode2)

#支川と本川の接続 # 96番目のエレメントが接続部
nodes2[-1].set_dn_element(elements[16])   
elements[16].set_up_node(nodes2[-1])

for elements_list in [elements,elements2]:
    for target_element in elements_list:
        # solver = my_module.Euler() ; num_stage=1
        solver = my_module.Runge_Kutta(my_module.RK4()) ; num_stage=4
        target_element.set_time_solver(solver)

for i in range(n_riv):
    elements[i].set_depth(0)
    nodes[i].set_flux(0)
for i in range(n_riv2):
    elements2[i].set_depth(0)
    nodes[i].set_flux(0)

time=0
maxt = 10

# 計算開始時刻を記録
start_time = timers.time()
#####################
### 時間発展ループ ###
#####################
count = 1
while time/3600 < maxt :
    H=[] ; Q=[]
    # boundary contion
    bc_upnode.set_flux(Qb[int(time // 3600)])
    bc_upnode2.set_flux(Qb2[int(time // 3600)])
    bc_dnelement.set_depth(elements[n_riv-1].get_depth())
    for stage in range(num_stage):
        for target_element in elements:
            target_element.solve_mass_equation(dt)# 質量保存側
        for target_element in elements2:
            target_element.solve_mass_equation(dt)# 質量保存側
        for target_node in nodes:
            target_node.solve_momentum_equation()
        tmp=0
        for target_node in nodes2:
            target_node.solve_momentum_equation()  
            tmp+=1
    for target_element in elements:
        H.append(target_element.get_elevation()+target_element.get_depth())
    for target_node in nodes:
        Q.append(target_node.get_flux)
    
    time += dt

    if np.mod(time/3600,1)==0:
        print("time:",time/3600,"  [h]",r"Q(m^3/s):","dt:",dt,"Qup:",Qb[int(time // 3600)],"Q[50]:",nodes[50].get_flux())
        fig = plt.figure(figsize=(11,4))
        plt.plot(x,H,label="water level")
        plt.plot(x,zb,label="elevatino")
        plt.legend() ; plt.grid()
        plt.savefig("figure/fig%04d.png"%(time/3600))
        plt.close()

# 計算終了時刻を記録
end_time = timers.time()
# 経過時間を表示
print(f"計算時間: {end_time - start_time:.3f} 秒")        