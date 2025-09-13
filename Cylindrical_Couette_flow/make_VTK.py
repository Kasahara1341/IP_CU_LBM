from os import sep
import numpy as np #,sys]
import csv
import pyvista as pv
import glob

from multiprocessing import Pool, cpu_count

def process_file(counter):
    # print(f"Processing loop {counter}")
    count=0 
    csv_file=open("input/initial.csv","r",encoding="utf-8",errors="",newline="")
    f = csv.reader(csv_file, delimiter=",", doublequote=True, lineterminator="\r\n", quotechar='"', skipinitialspace=True)
    for row in f:
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
    # 以下はVTK形式のデータ作成
    points = [] ; cells = [] ; cell_types = [] ; 
    sal_data = [] ; phi_data = [] ;rho_data= [] ; pressure_data = [] ; vel_data = [] ; force_data = []
    point_lookup = {}
    point_id = 0
    for i in range(count):
        p0 = (posx[i] - delX[i]/2.0, posy[i] - delY[i]/2.0, posz[i] - dz/2.0)
        p0 = (posx[i] - delX[i]/2.0, posy[i] - delY[i]/2.0, posz[i] - dz/2.0)
        p1 = (posx[i] + delX[i]/2.0, posy[i] - delY[i]/2.0, posz[i] - dz/2.0)
        p2 = (posx[i] + delX[i]/2.0, posy[i] + delY[i]/2.0, posz[i] - dz/2.0)
        p3 = (posx[i] - delX[i]/2.0, posy[i] + delY[i]/2.0, posz[i] - dz/2.0)
        p4 = (posx[i] - delX[i]/2.0, posy[i] - delY[i]/2.0, posz[i] + dz/2.0)
        p5 = (posx[i] + delX[i]/2.0, posy[i] - delY[i]/2.0, posz[i] + dz/2.0)
        p6 = (posx[i] + delX[i]/2.0, posy[i] + delY[i]/2.0, posz[i] + dz/2.0)
        p7 = (posx[i] - delX[i]/2.0, posy[i] + delY[i]/2.0, posz[i] + dz/2.0)
        corners = [p0, p1, p2, p3, p4, p5, p6, p7]
        corner_ids = []
        for p in corners:
            if p not in point_lookup:
                point_lookup[p] = point_id
                points.append(p)
                point_id += 1
            corner_ids.append(point_lookup[p])
        # セル定義（VTK_HEXAHEDRON = 12）
        cells.append([8] + corner_ids)
        cell_types.append(pv.CellType.HEXAHEDRON)
        sal_data.append(sal[i]) ; phi_data.append(phi[i]) ; rho_data.append(rho[i])
        pressure_data.append(pressure[i])
        vel_data.append([velx[i],vely[i],velz[i]])
        force_data.append([Fx[i],Fy[i],Fz[i]])   
    # numpy形式に変換
    points = np.array(points)
    cells_flat = np.hstack(cells).astype(np.int64)
    cell_types = np.array(cell_types)
    # UnstructuredGrid 作成
    grid = pv.UnstructuredGrid(cells_flat, cell_types, points)
    grid.cell_data["salinity"] = np.array(sal_data)
    grid.cell_data["phi"]      = np.array(phi_data)
    grid.cell_data["rho"]      = np.array(rho_data)
    grid.cell_data["pressure"] = np.array(pressure_data)
    grid.cell_data["velocity"] = np.array(vel_data)
    grid.cell_data["force"] = np.array(force_data)
    # 保存
    grid.save("vtu_CC/LBM%04d.vtu"%counter)

def main():
    files = sorted(glob.glob('./data/Fakhari_LBM*.csv'))  
    total_count=len(files)
    with Pool(cpu_count()) as pool:
        pool.map(process_file, range(total_count))
        # pool.map(process_file, range(300,301))

if __name__ == "__main__":
    main()

