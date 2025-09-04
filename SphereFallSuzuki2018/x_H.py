import matplotlib.pyplot as plt
import csv
import numpy as np

def read_csv(xarray,yarray,csv_file):
    f = csv.reader(open(csv_file+".csv", 'r'), delimiter=',', doublequote=True, lineterminator='\r\n',
                   quotechar='"', skipinitialspace=True)
    for row in f:
        if row[0] != 'x_H':
            yarray.append(float(row[1]))
            xarray.append(float(row[0]))

def make_figure(x_H,y_H,x_Href,y_Href, filename, ymin,ymax):
    plt.figure()
    plt.plot(x_H,y_H, marker='o', linestyle='-', color='b', label="my LBM")
    plt.plot(x_Href,y_Href, label="Suzuki & Yoshino(2018)")
    plt.title(filename)
    plt.ylabel(r'$u$[m/s]')
    plt.xlabel(r'$t$[s]')
    plt.grid()
    plt.legend()
    plt.ylim(ymin,ymax)
    plt.savefig(filename+'.png')    

x_H    = [] ; y_H    = []
x_Href = [] ; y_Href = []
read_csv(y_H,x_H,"x_H")
read_csv(x_Href,y_Href,"Sphere_SuzukiRe32") ; 
for i in range(len(x_Href)):
    y_Href[i] *= -1
    pass
make_figure(x_H,y_H,x_Href,y_Href, "x_H", -0.14, 0.0) #


# x_H    = [] ; y_H    = []
# x_Href = [] ; y_Href = []
# read_csv(x_H,y_H,"x_theta")
# read_csv(x_Href,y_Href,"Suzuki_1.1theta")
# make_figure(x_H,y_H,x_Href,y_Href, "x_theta", -30, 50) #

