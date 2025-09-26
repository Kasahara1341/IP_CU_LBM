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
# 条件に合うマスクの座標
zb = z[mask]
# 条件に合うマスクの川幅
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

zb  = zb[::-1]  ; x  = x[::-1]  ; width  = width[::-1]
zb2 = zb2[::-1] ; x2 = x2[::-1] ; width2 = width2[::-1]

ib = (zb[-1]-zb[0])/(x[-1] - x[0] )

# 流量データの読み込み
init_dt = 2.5
dt = init_dt
df = pd.read_csv("../input/Boundary2.csv")
Qb = list(df[df.keys()[1]])
Qb2 = list(df[df.keys()[2]])
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
bc_dnelement.set_up_node(nodes[n_riv-1])
nodes[n_riv-1].set_dn_element(bc_dnelement)

# 境界条件(支川)
bc_upnode2 = my_module.Node()
elements2[0].set_up_node(bc_upnode2)

#支川と本川の接続 # 96番目のエレメントが接続部
nodes2[-1].set_dn_element(elements[16])   
elements[16].set_up_node(nodes2[-1])

# 水系内で配列結合
# domein_hogeの変更は元のhogesやhoges2にも反映される
domein_elements = elements + elements2
domein_nodes    = nodes    + nodes2

# 時間積分方法の設定
for target_element in domein_elements:
    # solver = my_module.Euler() ; num_stage=1
    solver = my_module.Runge_Kutta(my_module.RK4()) ; num_stage=4
    # solver = my_module.Runge_Kutta(my_module.RK6()) ; num_stage=6
    target_element.set_time_solver(solver)

# 初期化
for i in range(n_riv):
    elements[i].set_depth(0)
    nodes[i].set_flux(0)
for i in range(n_riv2):
    elements2[i].set_depth(0)
    nodes[i].set_flux(0)


time=0.0
maxt = 100

# 計算開始時刻を記録
start_time = timers.time()
# 質量(水量)保存確認用
total_Qin=0 ; total_water=0 ; total_Qout=0
#####################
### 時間発展ループ ###
#####################
count = 1
Hall=[] ; Qall=[]
while time/3600 < maxt :
    H=[] ; Q=[]

    # boundary condition
    bc_upnode.set_flux(Qb[int(time//3600)])
    bc_upnode2.set_flux(Qb2[int(time//3600)])
    bc_dn_depth=elements[n_riv-1].get_depth()
    bc_dnelement.set_depth(bc_dn_depth)

    # solv_mass_equationなどはc++だがforで回す分多少遅くなる
    # ここではpythonコードで簡単に質量保存を確認するために使用する
    """qout_stage = []
    for stage in range(num_stage):
        qout_stage.append(nodes[-1].get_flux())
        for target_element in domein_elements:
            target_element.solve_mass_equation(dt)# 質量保存側
        for target_node in domein_nodes:
            target_node.solve_momentum_equation()
    # 4次Runge Kutta用の係数
    total_Qout+= (qout_stage[0]+2.0*qout_stage[1]+2.0*qout_stage[2]+qout_stage[3])*dt/6.0
    total_Qin += Qb[int(time // 3600)]  * dt
    total_Qin += Qb2[int(time // 3600)] * dt
    total_water=0 
    for target_element in domein_elements:
        tmp_depth = target_element.get_depth()
        tmp_width = target_element.get_width()
        tmp_length= target_element.get_length()
        total_water += tmp_depth*tmp_length*tmp_width
    # """
    # c++に渡して計算する．
    my_module.compute_all(domein_elements,domein_nodes,dt,num_stage)

    if Qb[int(time//3600)]>200:
        dt = 1
    else:
        dt = init_dt

    if 3600*count<time+dt:
        dt = 3600*count-time    
    time += dt

    if np.mod(time/3600,1)==0:
        for target_element in elements:
            H.append(target_element.get_elevation()+target_element.get_depth())
        for target_node in nodes:
            Q.append(target_node.get_flux())
        print("element[50]",elements[50].get_depth(),"domein_element[50]",domein_elements[50].get_depth())
        print("time:",time/3600,"  [h]",r"Q(m^3/s):","dt:",dt,"Qup:",Qb[int((time-dt) // 3600)],"Q[50]:",nodes[50].get_flux())
        print("total_Qin:",total_Qin," total_water:",total_water," total_Qout:",total_Qout,"sum:",total_Qout+total_water-total_Qin)
        print("numerical error = ",(total_Qout+total_water-total_Qin)/total_Qin)        
        count += 1
        Hall.append(H) 
        Qall.append(Q)
        dt = init_dt
        """fig = plt.figure(figsize=(11,4))
        plt.plot(x,H,label="water level")
        plt.plot(x,zb,label="elevatino")
        plt.legend() ; plt.grid()
        plt.savefig("figure/fig%04d.png"%(time/3600))
        plt.close()"""


# 計算終了時刻を記録
end_time = timers.time()
# 経過時間を表示
print(f"計算時間: {end_time - start_time:.3f} 秒")        

np.savetxt("../out/Hs.csv",Hall, fmt="%.6f", delimiter=",")
np.savetxt("../out/Qs.csv",Qall, fmt="%.6f", delimiter=",")