import matplotlib.pyplot as plt
import csv
import numpy as np

def read_csv(csv_file,plot_flag):
    xarray = [] ; yarray = []
    f = csv.reader(open(csv_file+".csv", 'r'), delimiter=',', doublequote=True, lineterminator='\r\n',
                   quotechar='"', skipinitialspace=True)
    for row in f:
        if row[0] != 'x_H':
            yarray.append(float(row[1]))
            xarray.append(float(row[0]))
    if plot_flag=="plot" :
        plt.plot(xarray, yarray, label=csv_file)
    elif plot_flag=="scatter" :
        plt.scatter(xarray, yarray, label=csv_file)


plt.figure()
plt.title(" Velocity profile ")
plt.ylabel(r'$u$[m/s]')
plt.xlabel(r'$t$[s]')
plt.grid()
plt.ylim(-0.14, 0.0)

read_csv("my_LBM_Re1.5","plot")
read_csv("my_LBM_Re4.1","plot")
read_csv("my_LBM_Re11.6","plot")
read_csv("my_LBM_Re31.9","plot")
read_csv("Sphere_SuzukiRe1.5","scatter") 
read_csv("Sphere_SuzukiRe4.1","scatter") 
read_csv("Sphere_SuzukiRe11","scatter") 
read_csv("Sphere_SuzukiRe32","scatter") 
plt.legend()

plt.savefig("x_H.png")

# x_H    = [] ; y_H    = []
# x_Href = [] ; y_Href = []
# read_csv(x_H,y_H,"x_theta")
# read_csv(x_Href,y_Href,"Suzuki_1.1theta")
# make_figure(x_H,y_H,x_Href,y_Href, "x_theta", -30, 50) #

