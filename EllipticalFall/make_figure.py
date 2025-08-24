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
            vel_abs= np.zeros(num_calc)
            Fx = np.zeros(num_calc) ; Fy = np.zeros(num_calc) ; Fz = np.zeros(num_calc) ; 
        elif row[0]=="x       ":
            pass
        else:
            posx     = np.append(posx     ,float(row[0]))      
            posy     = np.append(posy     ,float(row[1]))      
            posz     = np.append(posz     ,float(row[1]))      
            sal      = np.append(sal      ,float(row[6]))     
            velx     = np.append(velx     ,float(row[9])) 
            vely     = np.append(vely     ,float(row[10])) 
            velz     = np.append(velz     ,float(row[11])) 
            vel_abs  = np.append(vel_abs  , np.sqrt(velx[count]**2+vely[count]**2+velz[count]**2))
            count+=1
    max_vel = np.max(vel_abs)
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
                sal_avg[j, i] = np.mean(vel_abs[mask])  # y方向の平均
    X, Z = np.meshgrid(x_unique, z_unique)

    plt.figure(figsize=(4,2))
    plt.pcolormesh(X, Z, sal_avg, shading='auto', cmap='jet',vmin=0,vmax=0.03)
    # plt.tricontourf(posx,posz,sal,levels=100,cmap="jet",vmin=0,vmax=24.26)
    plt.xlabel("time step %04d"%counter)
    plt.colorbar()
    plt.savefig("figure/velx_%04d.png"%counter,bbox_inches='tight', dpi=300)

    """plt.figure(figsize=(8,1))
    plt.tricontourf(posx,posz,pressure,levels=1,cmap="jet",vmin=0,vmax=25)
    plt.xlabel("time step %04d"%counter)
    # plt.colorbar()
    plt.savefig("figure/pressure_%04d.png"%counter,bbox_inches='tight', dpi=300)"""
def main():
    files = sorted(glob.glob('./data/Fakhari_LBM*.csv'))  
    total_count=len(files)
    with Pool(cpu_count()) as pool:
        pool.map(process_file, range(50,total_count))
        # pool.map(process_file, range(0,10))
    #make gif and mp4
    files = sorted(glob.glob('./figure/velx*.png'))  
    titlegif="./karman_vortex"  
    images = list(map(lambda file : Image.open(file) , files))
    images[0].save( titlegif+'.gif' , save_all = True , append_images = images[1:] , duration = 100 , loop = 0)
    """files = sorted(glob.glob('./figure/pressure*.png'))  
    titlegif="./pressure"  
    images = list(map(lambda file : Image.open(file) , files))
    images[0].save( titlegif+'.gif' , save_all = True , append_images = images[1:] , duration = 100 , loop = 0)"""

if __name__ == "__main__":
    main()



