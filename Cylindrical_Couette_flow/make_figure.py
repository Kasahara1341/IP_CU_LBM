from os import sep
import numpy as np #,sys]
import csv
import pyvista as pv
import matplotlib.pyplot as plt
import glob

from multiprocessing import Pool, cpu_count

def process_file(counter):
    # print(f"Processing loop {counter}")
    count=0 
    csv_file=open("input/initial.csv","r",encoding="utf-8",errors="",newline="")
    f = csv.reader(csv_file, delimiter=",", doublequote=True, lineterminator="\r\n", quotechar='"', skipinitialspace=True)
    for row in f:
        nx   =int(row[0]) ; nz = int(row[2])
        dz   =float(row[3])
        ratiox=float(row[5])
        ratioy=float(row[6])
        dx=dz*ratiox
        if(int(row[1])/ratioy<1):   ratioy==1
        dy=dz*ratioy
        break
    with open(f"./data/Fakhari_LBM{counter:04d}.csv", "r", encoding="utf-8", errors="", newline="") as csv_file:
        f = csv.reader(csv_file, delimiter=",", doublequote=True, lineterminator="\r\n",
                       quotechar='"', skipinitialspace=True)
        for row in f:
            if row[0]=="num_calc":
                num_calc=int(row[1])
                sal = np.zeros(num_calc) ; phi = np.zeros(num_calc) ; rho = np.zeros(num_calc) ; posx=np.zeros(num_calc) ; posy=np.zeros(num_calc) ; posz=np.zeros(num_calc)
                pressure = np.zeros(num_calc) ; velx = np.zeros(num_calc) ; vely = np.zeros(num_calc) ; velz = np.zeros(num_calc) 
                Fx = np.zeros(num_calc) ; Fy = np.zeros(num_calc) ; Fz = np.zeros(num_calc) ; delX = np.zeros(num_calc) ; delY = np.zeros(num_calc)
            elif row[0]=="x       ":
                pass
            else:
                posx[count]      = float(row[0])
                posy[count]      = float(row[1])
                posz[count]      = float(row[2])
                delX[count]      = float(row[3])
                delY[count]      = float(row[4])
                pressure [count] = float(row[5])
                sal [count]      = float(row[6])
                phi [count]      = float(row[7])
                rho [count]      = float(row[8])
                velx[count]      = float(row[9])
                vely[count]      = float(row[10])
                velz[count]      = float(row[11])
                Fx  [count]      = float(row[12])
                Fy  [count]      = float(row[13])
                Fz  [count]      = float(row[14])
                count+=1
    velU = np.zeros(nx) ; x = np.zeros(nx)
    exact_vel = np.zeros(nx) 
    Re = 6 ; nu = 1.0*10**(-4)
    L = 1.0 ; R1 = L/4.8 ; R2 = L/2.4
    vel_theta = nu*Re/(L/4.8)**2  ; vel_U = vel_theta*R1
    x_line0 = [0.5*L-R1,0.5*L-R1] ; x_line1 = [0.5*L+R1,0.5*L+R1]
    x_line3 = [0.5*L-R2,0.5*L-R2] ; x_line4 = [0.5*L+R2,0.5*L+R2]

    y_line = [-1,1]
    error = 0
    for i in range(nx):
        velU[i] = velz[int(nz/2+1 + i*nz)]
        x[i] = posx[int(nz/2+1 + i*nz)]
        R0 = posx[int(nz/2+1 + i*nz)] - L/2.0
        # exact_vel[i] = vel_theta * (R0/R2-R2/R0)/(R1/R2-R2/R1) * R0 # seta
        exact_vel[i] = -vel_theta*R1 * (R1*R2/(R2**2-R1**2))*(R2/R0 - R0/R2) # Zhang et al. 2019
        if x[i] < 0.5*L-R2 or x[i] > 0.5*L+R2 or (x[i]>0.5*L-R1 and x[i]<0.5*L+R1):
            velU[i] = 0
            exact_vel[i] = 0
        error += (exact_vel[i] - velU[i])/(exact_vel[i]+10**(-10))
    print("error = ",error/float(len(exact_vel)))
    plt.figure()
    plt.plot(x,velU,label="my LBM")
    plt.plot(x,exact_vel,label="exact")
    plt.plot(x_line0,y_line,color="black")
    plt.plot(x_line1,y_line,color="black")
    plt.plot(x_line3,y_line,color="black")
    plt.plot(x_line4,y_line,color="black")
    plt.ylim(-vel_U,vel_U)
    plt.xlabel(r"$x$(m)")
    plt.ylabel(r"$u$(m/s)")
    plt.legend()
    plt.savefig("figure/cylindrical%04d.png"%counter)


def main():
    files = sorted(glob.glob('./data/Fakhari_LBM*.csv'))  
    total_count=len(files)
    with Pool(cpu_count()) as pool:
        # pool.map(process_file, range(total_count))
        pool.map(process_file, [100])

if __name__ == "__main__":
    main()

