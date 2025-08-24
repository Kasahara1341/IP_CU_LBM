import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
import csv
import glob
from multiprocessing import Pool, cpu_count

def process_file(counter):
    # print("loop %d"%counter)
    csv_file=open("./data/Fakhari_LBM%04d.csv"%(counter),"r",encoding="utf-8",errors="",newline="")
    f = csv.reader(csv_file, delimiter=",", doublequote=True, lineterminator="\r\n", quotechar='"', skipinitialspace=True)
    count=0
    for row in f:
        if row[0]=="num_calc":
            num_calc=0
            sal = np.zeros(num_calc) ; posx=np.zeros(num_calc) ; posy=np.zeros(num_calc) ; posz=np.zeros(num_calc)
            pressure = np.zeros(num_calc) ; velx = np.zeros(num_calc) ; vely = np.zeros(num_calc) ; velz = np.zeros(num_calc) 
            Fx = np.zeros(num_calc) ; Fy = np.zeros(num_calc) ; Fz = np.zeros(num_calc) ; 
        elif row[0]=="x       ":
            pass
        else:
            posx     = np.append(posx     ,float(row[0]))      
            posy     = np.append(posy     ,float(row[1]))      
            posz     = np.append(posz     ,float(row[2]))      
            sal      = np.append(sal      ,float(row[6]))     
            velx     = np.append(velx     ,float(row[9])) 
            count+=1
    xz = np.stack([posx,posz],axis=1)
    xz_unique = np.unique(xz, axis=0)

    x_unique = np.unique(posx)
    z_unique = np.unique(posz)

    sal_avg = np.full((len(z_unique), len(x_unique)),np.nan)
    for i, x in enumerate(x_unique):
        for j, z in enumerate(z_unique):
            # (x,z)が一致するインデックスを取得
            mask = (posx == x) & (posz == z)
            if np.any(mask):
                sal_avg[j, i] = np.mean(sal[mask])  # y方向の平均
    X, Z = np.meshgrid(x_unique, z_unique)

    plt.figure(figsize=(8,1))
    plt.pcolormesh(X, Z, sal_avg, shading='auto', cmap='jet',vmin=0,vmax=1.22)
    # plt.tricontourf(posx,posz,sal,levels=100,cmap="jet",vmin=0,vmax=24.26)
    plt.xlabel("time step %04d"%counter)
    plt.colorbar()
    plt.savefig("figure/sal_%04d.png"%counter,bbox_inches='tight', dpi=300)

    """plt.figure(figsize=(8,1))
    plt.tricontourf(posx,posz,pressure,levels=1,cmap="jet",vmin=0,vmax=25)
    plt.xlabel("time step %04d"%counter)
    # plt.colorbar()
    plt.savefig("figure/pressure_%04d.png"%counter,bbox_inches='tight', dpi=300)"""
def main():
    files = sorted(glob.glob('./data/Fakhari_LBM*.csv'))  
    total_count=len(files)
    with Pool(cpu_count()) as pool:
        pool.map(process_file, range(total_count))
        # pool.map(process_file, range(0,10))
    #make gif and mp4
    files = sorted(glob.glob('./figure/sal*.png'))  
    titlegif="./Lock_exchange"  
    images = list(map(lambda file : Image.open(file) , files))
    images[0].save( titlegif+'.gif' , save_all = True , append_images = images[1:] , duration = 100 , loop = 0)
    """files = sorted(glob.glob('./figure/pressure*.png'))  
    titlegif="./pressure"  
    images = list(map(lambda file : Image.open(file) , files))
    images[0].save( titlegif+'.gif' , save_all = True , append_images = images[1:] , duration = 100 , loop = 0)"""

if __name__ == "__main__":
    main()



