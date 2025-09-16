import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
import csv
import glob
import pandas as pd
from multiprocessing import Pool, cpu_count

def process_file(counter):
    # print("loop %d"%counter)
    center_x = 1.4
    csv_file=open("./data/Fakhari_LBM%04d.csv"%(counter),"r",encoding="utf-8",errors="",newline="")
    f = csv.reader(csv_file, delimiter=",", doublequote=True, lineterminator="\r\n", quotechar='"', skipinitialspace=True)
    for row in f:
        if row[0]=="num_calc":
            num_calc=0
            sal = np.zeros(num_calc) ; posx=np.zeros(num_calc) ; posy=np.zeros(num_calc) ; posz=np.zeros(num_calc)
        elif row[0]=="x       ":
            pass
        # else:
        elif float(row[0]) > center_x-0.51 and float(row[0]) <center_x+0.51:
            posx     = np.append(posx     ,float(row[0]))      
            posy     = np.append(posy     ,float(row[1]))      
            posz     = np.append(posz     ,float(row[2]))      
            sal      = np.append(sal      ,float(row[6]))      
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

    length = 1.1
    plt.figure(figsize=(length*5,length))
    plt.pcolormesh(X, Z, sal_avg, shading='auto', cmap='jet',vmin=0,vmax=24.26)
    plt.xlim(center_x-0.5,center_x+0.5)
    plt.xlabel("x(m)")
    plt.ylabel("z(m)")
    # plt.xlabel("time step %04d"%counter)
    plt.colorbar(label="salinity")
    plt.savefig("figure/ave_sal_%04d.png"%counter,bbox_inches='tight', dpi=300) # """

def main():
    files = sorted(glob.glob('./data/Fakhari_LBM*.csv'))  
    total_count=len(files)
    total_count= [125, 129, 134, 138, 142, 145]
    with Pool(cpu_count()) as pool:
        pool.map(process_file, total_count)

if __name__ == "__main__":
    main()



