from os import sep
import numpy as np #,sys]
import csv
import pyvista as pv
import pandas as pd
import glob

from multiprocessing import Pool, cpu_count

def process_file(counter):
    # print(f"Processing loop {counter}")
    count=0 
    df = pd.read_csv(f"./IB_point/ibm_{counter:04d}.csv")
    # 座標を points に変換
    points = df[["x", "y", "z"]].to_numpy()

    # PolyData 作成
    pdata = pv.PolyData(points)

    # ベクトルを追加
    vectors = df[["ax", "ay", "az"]].to_numpy()
    pdata["acceleration"] = vectors

    vectors = df[["velx", "vely", "velz"]].to_numpy()
    pdata["velocity"] = vectors

    # VTP形式で保存
    pdata.save(f"vtu_KM/IB_point{counter:04d}.vtp")    

def main():
    files = sorted(glob.glob('./IB_point/ibm*.csv'))  
    total_count=len(files)
    with Pool(cpu_count()) as pool:
        # pool.map(process_file, range(0, 200, 10))
        pool.map(process_file, range(0,total_count,1))
        # pool.map(process_file, range(300,301))

if __name__ == "__main__":
    main()

