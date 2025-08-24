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
start=0
end=200
type=1 # AMR適用時は2 その他は1
salt=0.5
heiB=[] ; ti=[] 
#csvファイルから値を取得
csv_file=open("./input/initial.csv","r",encoding="utf-8",errors="",newline="")
f = csv.reader(csv_file, delimiter=",", doublequote=True, lineterminator="\r\n", quotechar='"', skipinitialspace=True)
for row in f:   #行ごとに回す
    if row[0]=='nx':  #2行目が来たら(値が出きったら)ループから抜ける
        break
    for i in range(18):
        if row[i]=="":
            break
        row[i]=float(row[i])
    ratio=row[5 ]
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
for count1 in range(start,end): 
    # print(count1)
    csv_file=open("./data/Fakhari_LBM%04d.csv"%(count1),"r",encoding="utf-8",errors="",newline="")
    f = csv.reader(csv_file, delimiter=",", doublequote=True, lineterminator="\r\n", quotechar='"', skipinitialspace=True)
    tmp=0 
    phi  = [] ; x = [] ; y= [] ; z= [] 
    for row in f:   #行ごとに回す
        if row[0]=='':  #空白が来たら(値が出きったら)ループから抜ける
            break
        if row[0]!="num_calc" and row[0]!="x       ":
            phi.append(float(row[7])) ; x.append(float(row[0])) ; y.append(float(row[1])) ; z.append(float(row[2]))
            tmp+=1 
    csv_file.close()
    for count2 in range(len(phi)):
        if(x[count2]>=3.02 and x[count2]<=3.02+dx*ratio ):
        # if(x[count2]>=3.02 and x[count2]<=3.02+dx*ratio and y[count2]==0.165):  # center
            if(phi[count2]>salt and phi[count2+type]<salt):
                heiB.append((z[count2]+np.abs(salt-phi[count2])/(np.abs(salt-phi[count2+type])+np.abs(salt-phi[count2]))*dx/type-0.087)*100)
                break
        else:
            continue
    ti.append(2*count1)
    print("lin(heiB)=%d  lin(ti)=%d"%(len(heiB),len(ti)))
    # print("loop=%d"%(count1))
    # print("heiBght=%f, time=%f"%(heiB[int(count1)],ti[int(count1)]))
data=np.zeros((len(ti),2))
for j in range(len(ti)):
    data[j][0]=ti[j] ; data[j][1]=heiB[j] 
np.savetxt("./make_result/1_%d_B.csv"%ratio,data, delimiter=",")
np.savetxt("./make_result/1_0.csv",data, delimiter=",")
print("done")
