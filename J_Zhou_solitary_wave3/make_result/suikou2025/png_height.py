#from PIL import Image
# 0514
"""
ここではx,y,z面に壁関数を適用した結果をr=3,12で可視化する．r=3の場合
1200,60,58,0.005,0.0025,3.0,12.0,-9.81,1,600,100,
nx,ny,nz,dx,dt,ratiox,ratioy,gz,nu,save_interval,total_count,

"""
from os import sep
import numpy as np #,sys]
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import japanize_matplotlib
import csv 

def read_csv(filename,time,height,number):
    csv_file=open(filename+".csv","r",encoding="utf-8",errors="",newline="")
    f = csv.reader(csv_file, delimiter=",", doublequote=True, lineterminator="\r\n", quotechar='"', skipinitialspace=True)
    count=0
    for row in f:   #行ごとに回す
        if row[0]=="" :
            break
        if row[0]!='X' :  #float検知時のみ処理を行う
            time[number][count]=row[0] ; height[number][count]=row[1] 
            count+=1
    csv_file.close()

start=0
end=601
var=[]
x=[] ; y=[] ; z=[]
hei=[] ; ti=[] 
time=np.zeros((10,605)) ; height=np.zeros((10,605))

read_csv("make_result/suikou2025/Horn_2001_B",time,height,0)
read_csv("make_result/suikou2025/1_3_B",time,height,1)
read_csv("make_result/suikou2025/1_12_B",time,height,2)

for k in range(7):
    for i in range(len(time[0])):
        # height[k][i]/=100
        height[k][i]*=1
for k in range(7):
    for i in range(len(time[0])):
        if time[k][i]==0.0 and i!=0:
            for j in range(i,len(time[0])):
                time[k][j]=time[k][j-1] ; height[k][j]=height[k][j-1]
########### plot ######################################
figure(figsize=(8,4))
plt.xlim(0,400)
plt.ylabel("height (m)",fontsize=18)
plt.xticks(fontsize=18)
plt.ylim(-5.0,5.0)
plt.yticks(fontsize=18)
plt.xlabel("time (s)",fontsize=18)
plt.plot(time[1],height[1],"-o",markersize=4,color='gray',label="ケース1")
plt.plot(time[2],height[2],color='gray',label="ケース2")
plt.plot(time[0],height[0],color="black",label="Horn et al. (2001)")
plt.legend(loc="upper left",fontsize="15")
plt.grid()
plt.savefig("make_result/suikou2025/interface_HornB.png",bbox_inches="tight",pad_inches=0.3)# ,bbox_inches="tight",pad_inches=0.3)   # dpi=〇(解像度)bbox.inches="tight"(完成後の画像の余白を調整),pad-igches=0.3(横のカラーバーの位置調整)
print("done")
