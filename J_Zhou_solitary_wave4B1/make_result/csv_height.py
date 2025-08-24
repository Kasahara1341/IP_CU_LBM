#from PIL import Image
import glob
import sys
from os import sep
import numpy as np #,sys]
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import csv 
import matplotlib.ticker as ticker

def get_hei(nxyz,pointx,dx,salt,salinity,x,z,hei):
    for count2 in range(nxyz):
        if(x[count2]>=pointx and x[count2]<=pointx+dx*12 ):
            if(salinity[count2]>salt and salinity[count2+1]<salt):
                hei.append((z[count2]+np.abs(salt-salinity[count2])/(np.abs(salt-salinity[count2+1])+np.abs(salt-salinity[count2]))*dx-0.34)/0.34)
                break
        else:
            continue

start=0
end=200
type=1 # AMR適用時は2 その他は1
salt=0.5
hei1=[] ; hei2=[] ; hei3=[] ; ti=[] 
#csvファイルから値を取得
csv_file=open("./input/initial.csv","r",encoding="utf-8",errors="",newline="")
f = csv.reader(csv_file, delimiter=",", doublequote=True, lineterminator="\r\n", quotechar='"', skipinitialspace=True)
ratio= int(input("ratio: "))
for row in f:   #行ごとに回す
    if row[0]=='nx':  #2行目が来たら(値が出きったら)ループから抜ける
        break
    for i in range(18):
        if row[i]=="":
            break
        row[i]=float(row[i])
    row[0]=int(row[0]/ratio) ; row[1]=int(row[1]/ratio) ; row[2]=int(row[2]) 
    nx   =int(row[0 ] )
    ny   =int(row[1 ] )
    nz   =int(row[2 ])
    dx   =row[3 ]
    dt   =row[4 ]
    loop1=row[8]
if(ny<1):
    ny=1 
print("dx=%f dt=%f"%(dx,dt))
print("ratio=%f type=%f"%(ratio,type))
print("nx=%d ny=%d nz=%d"%(nx,ny,nz))
tmp=0 
end = glob.glob("./data/Fakhari_LBM*.csv")
end = len(end)
# end = 35
for count1 in range(start,end): 
    # print(count1)
    csv_file=open("./data/Fakhari_LBM%04d.csv"%(count1),"r",encoding="utf-8",errors="",newline="")
    f = csv.reader(csv_file, delimiter=",", doublequote=True, lineterminator="\r\n", quotechar='"', skipinitialspace=True)
    tmp=0 
    salinity  = [] ; x = [] ; y= [] ; z= [] 
    for row in f:   #行ごとに回す
        if row[0]=='':  #空白が来たら(値が出きったら)ループから抜ける
            break
        if row[0]!="num_calc" and row[0]!="x       ":
            for i in range(8):
                row[i] = float(row[i])
            salinity.append(row[7]) ; x.append(row[0]) ; y.append(row[1]) ; z.append(row[2])
            tmp+=1 
    csv_file.close()
    get_hei(len(salinity),4.0712,dx,0.5,salinity,x,z,hei1)
    get_hei(len(salinity),4.6688,dx,0.5,salinity,x,z,hei2)
    get_hei(len(salinity),5.1224,dx,0.5,salinity,x,z,hei3)
    ti.append(0.1*count1*np.sqrt(9.81/0.34) - 2.75)
    # print("lin(hei)=%d  lin(ti)=%d"%(len(hei1),len(ti)),"hei2 hei3",len(hei2),len(hei3))
    # print("loop=%d"%(count1))
    # print("heiBght=%f, time=%f"%(heiB[int(count1)],ti[int(count1)]))
data1=np.zeros((len(ti),2)) ; data2=np.zeros((len(ti),2)) ; data3=np.zeros((len(ti),2))
for j in range(len(ti)):
    data1[j][0]=ti[j] ; data1[j][1]=hei1[j] 
    data2[j][0]=ti[j] ; data2[j][1]=hei2[j] 
    data3[j][0]=ti[j] ; data3[j][1]=hei3[j] 
np.savetxt("./make_result/gauge1_%d_B.csv"%ratio,data1, delimiter=",")
np.savetxt("./make_result/gauge2_%d_B.csv"%ratio,data2, delimiter=",")
np.savetxt("./make_result/gauge3_%d_B.csv"%ratio,data3, delimiter=",")
print("done")
