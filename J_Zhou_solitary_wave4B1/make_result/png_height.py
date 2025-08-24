#from PIL import Image
from os import sep
import numpy as np #,sys]
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
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

def make_fig(time,height,idx1,idx2,ratio,num):
    figure(figsize=(8,4))
    plt.xlim(0,16.25)
    plt.ylabel(r"$\eta/d_0$",fontsize=18)
    plt.xticks(fontsize=18)
    plt.ylim(-0.2,0.6)
    plt.yticks(fontsize=18)
    plt.xlabel(r"$t\sqrt{g/d_0}$",fontsize=18)
    plt.grid()
    # plt.plot(time[3],height[3],color='gray',label="ratio=%d"%ratio)
    plt.plot(time[idx1],height[idx1],color="black",label="LBM %d"%num)
    plt.plot(time[idx2],height[idx2],color="gray",label="exp.")
    plt.legend(loc="upper left",fontsize="15")
    plt.savefig("make_result/gauge%d_%dB.png"%(num,ratio),bbox_inches="tight",pad_inches=0.3)# ,bbox_inches="tight",pad_inches=0.3)   # dpi=〇(解像度)bbox.inches="tight"(完成後の画像の余白を調整),pad-igches=0.3(横のカラーバーの位置調整)


start=0
end=601
var=[]
x=[] ; y=[] ; z=[]
hei=[] ; ti=[] 
time=np.zeros((10,605)) ; height=np.zeros((10,605))

ratio= int(input("ratio: "))
# ratio=6
read_csv("make_result/gauge1_%d_B"%ratio,time,height,0)
read_csv("make_result/gauge2_%d_B"%ratio,time,height,1)
read_csv("make_result/gauge3_%d_B"%ratio,time,height,2)
read_csv("make_result/case1_gauge1",time,height,3)
read_csv("make_result/case1_gauge2",time,height,4)
read_csv("make_result/case1_gauge3",time,height,5)

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
make_fig(time,height,0,3,ratio,1)
make_fig(time,height,1,4,ratio,2)
make_fig(time,height,2,5,ratio,3)
print("done")
